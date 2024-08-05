/**
* \file       AVM.cpp
* \date        September 24, 2021
* \version     v0.8
* \copyright   <2009-2021> Forschungszentrum Jülich GmbH. All rights reserved.
*
* \section License
* This file is part of JuPedSim.
*
* JuPedSim is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* any later version.
*
* JuPedSim is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with JuPedSim. If not, see <http://www.gnu.org/licenses/>.
*
* \section Description
* Implementation of first-order model
* Anticipation velocity model: Qiancheng (8)
*
*
**/

# define NOMINMAX
#include "../pedestrian/Pedestrian.h"
#include "../mpi/LCGrid.h"
#include "../geometry/Wall.h"
#include "../geometry/SubRoom.h"

#include "AVMModel.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads()  1
#endif

using std::vector;
using std::string;

AVMModel::AVMModel(std::shared_ptr<DirectionStrategy> dir, int model,
    double aped, double Dped, double awall, double Dwall,
    double Ts, double Td,
    double AntiT, bool calpha,
    double lb, double rb, double ub, double db, double co)
{
    _direction = dir;
    _Model = model;
    // Force_rep_PED Parameter
    _aPed = aped;
    _DPed = Dped;

    // Force_rep_WALL Parameter
    _aWall = awall;
    _DWall = Dwall;

    // GCVM Parameter
    _Ts = Ts; // Speed module
    _Td = Td; // Direction module

    // AVM Parameter
    _AntiTime = AntiT;
    _ConstantAlpha = calpha;

    // Boundary Case
    _LeftBoundary = lb;
    _RightBoundary = rb;
    _UpBoundary = ub;
    _DownBoundary = db;
    _CutOff = co;
}


AVMModel::~AVMModel()
{

}

string AVMModel::GetDescription()
{
    string rueck;
    char tmp[CLENGTH];

    sprintf(tmp, "\t\ta: \t\tPed: %f \tWall: %f\n", _aPed, _aWall);
    rueck.append(tmp);
    sprintf(tmp, "\t\tD: \t\tPed: %f \tWall: %f\n", _DPed, _DWall);
    rueck.append(tmp);
    return rueck;
}

// Initialize the direction of pedestrians according to the target 
bool AVMModel::Init(Building* building)
{

    if (auto dirff = dynamic_cast<DirectionFloorfield*>(_direction.get())) {
        Log->Write("INFO:\t Init DirectionFloorfield starting ...");
        double _deltaH = building->GetConfig()->get_deltaH();
        double _wallAvoidDistance = building->GetConfig()->get_wall_avoid_distance();
        bool _useWallAvoidance = building->GetConfig()->get_use_wall_avoidance();
        dirff->Init(building, _deltaH, _wallAvoidDistance, _useWallAvoidance);
        Log->Write("INFO:\t Init DirectionFloorfield done");
    }

    if (auto dirlocff = dynamic_cast<DirectionLocalFloorfield*>(_direction.get())) {
        Log->Write("INFO:\t Init DirectionLOCALFloorfield starting ...");
        double _deltaH = building->GetConfig()->get_deltaH();
        double _wallAvoidDistance = building->GetConfig()->get_wall_avoid_distance();
        bool _useWallAvoidance = building->GetConfig()->get_use_wall_avoidance();
        dirlocff->Init(building, _deltaH, _wallAvoidDistance, _useWallAvoidance);
        Log->Write("INFO:\t Init DirectionLOCALFloorfield done");
    }

    if (auto dirsublocff = dynamic_cast<DirectionSubLocalFloorfield*>(_direction.get())) {
        Log->Write("INFO:\t Init DirectionSubLOCALFloorfield starting ...");
        double _deltaH = building->GetConfig()->get_deltaH();
        double _wallAvoidDistance = building->GetConfig()->get_wall_avoid_distance();
        bool _useWallAvoidance = building->GetConfig()->get_use_wall_avoidance();
        dirsublocff->Init(building, _deltaH, _wallAvoidDistance, _useWallAvoidance);
        Log->Write("INFO:\t Init DirectionSubLOCALFloorfield done");
    }

    const vector< Pedestrian* >& allPeds = building->GetAllPedestrians();
    size_t peds_size = allPeds.size();
    for (unsigned int p = 0; p < peds_size; p++)
    {
        Pedestrian* ped = allPeds[p];
        double cosPhi, sinPhi;
        //a destination could not be found for that pedestrian
        if (ped->FindRoute() == -1) {
            Log->Write(
                "ERROR:\tAVMModel::Init() can not initialise route. ped %d is deleted in Room %d %d.\n", ped->GetID(), ped->GetRoomID(), ped->GetSubRoomID());
            building->DeletePedestrian(ped);
            p--;
            peds_size--;
            continue;
        }
        Point target = ped->GetExitLine()->ShortestPoint(ped->GetPos());
        Point d = target - ped->GetPos();
        double dist = d.Norm();
        if (dist != 0.0) {
            cosPhi = d._x / dist;
            sinPhi = d._y / dist;
        }
        else {
            Log->Write(
                "ERROR: \tallPeds::Init() cannot initialise phi! "
                "dist to target is 0\n");
            return false;
        }

        ped->InitV0(target);

        JEllipse E = ped->GetEllipse();
        E.SetCosPhi(cosPhi);
        E.SetSinPhi(sinPhi);
        ped->SetEllipse(E);
    }
    return true;
}

// Compute the status of pedestrians at next timestep
void AVMModel::ComputeNextTimeStep(double current, double deltaT, Building* building, int periodic)
{
    const vector< Pedestrian* >& allPeds = building->GetAllPedestrians();
    vector<Pedestrian*> pedsToRemove;
    pedsToRemove.reserve(500);
    unsigned long nSize;
    nSize = allPeds.size();
    //---------------------------------------------------------------
    vector< Point > resultAcc = vector<Point >();
    resultAcc.reserve(nSize);

    vector< Point > resultDir = vector<Point >();
    resultDir.reserve(nSize);

    vector<Point> resultForce = vector<Point>();
    resultForce.reserve(nSize);

    int start = 0;
    int end = nSize - 1;

    // calculate direction of each agent
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped1 = allPeds[p];
        Point p1 = ped1->GetPos();
        Room* room1 = building->GetRoom(ped1->GetRoomID());
        SubRoom* subroom1 = room1->GetSubRoom(ped1->GetSubRoomID());
        // Effect from neighbours
        Point IniDirection = DesireDirection(ped1, room1);//desired moving direction, direction3, and using core_size here.
        vector<Pedestrian*> neighbours;
        building->GetGrid()->GetNeighbourhood(ped1, neighbours);
        int size = (int)neighbours.size();

        Point repPed = Point(0, 0); //CSM effect
        // Calculate the influence from neighbours
        for (int i = 0; i < size; i++)
        {
            Pedestrian* ped2 = neighbours[i];
            Point p2 = ped2->GetPos();
            //Check if two pedestrians can see each other
            SubRoom* subroom2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());
            vector<SubRoom*> emptyVector;
            emptyVector.push_back(subroom1);
            emptyVector.push_back(subroom2);
            bool isVisible = building->IsVisible(p1, p2, emptyVector, false);
            if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom1->IsDirectlyConnectedWith(subroom2))
            {
                // Different model with different effect
                switch (GetModel())
                {
                case 0:
                    repPed += ForceRepPedCSM(ped1, ped2, building, periodic);
                    break;
                case 1:
                    repPed += ForceRepPedGCVM(ped1, ped2, building, periodic);
                    break;
                case 2:
                    repPed += ForceRepPedAVM(ped1, ped2, building, periodic);
                    break;
                default:
                    break;
                }
            }
        } //for i

        // Calculate the influence from walls
        Point repWall = ForceRepRoom(ped1, subroom1);

        // Turning process
        Point direction = Point(0, 0);
        Point dDirection = IniDirection + repPed + repWall;
        switch (GetModel())
        {
        case 0:
            direction = dDirection;
            break;
        case 1:
        case 2:
            Point aDirection = ped1->GetMoveDirection();
            Point AccTu = Point(0, 0);
            double angleTau = GetTd();
            AccTu = (dDirection.Normalized() - aDirection) / angleTau;
            direction = aDirection + AccTu * deltaT;
            if (IniDirection.ScalarProduct(dDirection)*IniDirection.ScalarProduct(aDirection) < 0)
            {
                direction = dDirection;
            }
            break;
        }
        direction = direction.Normalized();
        resultDir.push_back(direction);
    }

    // Update direction of each pedestrian
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped = allPeds[p];
        Point direction = resultDir[p];
        ped->SetMoveDirection(direction);
    }

    // Calculate speed
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped1 = allPeds[p];
        Room* room1 = building->GetRoom(ped1->GetRoomID());
        SubRoom* subroom1 = room1->GetSubRoom(ped1->GetSubRoomID());

        vector<Pedestrian*> neighbours;
        building->GetGrid()->GetNeighbourhood(ped1, neighbours);
        int size = (int)neighbours.size();

        vector< my_pair > spacings = vector<my_pair >();
        spacings.reserve(size);

        //Calculating spacing in front -------------------------------------------------------------------------------------------------------
        for (int i = 0; i < size; i++)
        {
            Pedestrian* ped2 = neighbours[i];
            SubRoom* subroom2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());
            if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom1->IsDirectlyConnectedWith(subroom2))
            {
                spacings.push_back(GetSpacing(ped1, ped2, periodic));
                //spacings.push_back(GetSpacingEllipse(ped1, ped2, periodic));
            }
        }

        // Calculate min spacing
        std::sort(spacings.begin(), spacings.end(), sort_pred_agcvm());
        double spacing = spacings.size() == 0 ? 100 : spacings[0].first;
        double spacingWall = GetSpacingRoom(ped1, subroom1);
        // neighbour subroom needs to be considered
        for (const auto & subr : subroom1->GetNeighbors())
        {
            double spacingWallNext = GetSpacingRoom(ped1, subr);
            spacingWall = spacingWall > spacingWallNext ? spacingWallNext : spacingWall;
        }
        spacing = spacing < spacingWall ? spacing : spacingWall;
        // Optimal speed function
        Point speed;
        Point ei = ped1->GetMoveDirection();
        //speed = ei * OptimalSpeed(ped1, spacing);
        // Here this speed model can create space
        speed = ei * PushSpeed(ped1, spacing);
        resultAcc.push_back(speed);
    }

    // Calculate the pushing force, from both neighbors and walls
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped1 = allPeds[p];
        Point p1 = ped1->GetPos();
        Room* room1 = building->GetRoom(ped1->GetRoomID());
        SubRoom* subroom1 = room1->GetSubRoom(ped1->GetSubRoomID());
        vector<Pedestrian*> neighbours;
        building->GetGrid()->GetNeighbourhood(ped1, neighbours);
        int size = (int)neighbours.size();

        Point contactForce = Point(0, 0); //Contacting results in force.
        // Calculate the influence from neighbours
        for (int i = 0; i < size; i++)
        {
            Pedestrian* ped2 = neighbours[i];
            Point p2 = ped2->GetPos();
            //Check if two pedestrians can see each other
            SubRoom* subroom2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());
            vector<SubRoom*> emptyVector;
            emptyVector.push_back(subroom1);
            emptyVector.push_back(subroom2);
            bool isVisible = building->IsVisible(p1, p2, emptyVector, false);
            if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom1->IsDirectlyConnectedWith(subroom2))
            {
                // Different model with different effect
                contactForce += ForceConPed(ped1, ped2, building, periodic);
            }
            contactForce += ForceConRoom(ped1, subroom1);
        } //for i
         resultForce.push_back(contactForce);
    }


    //Update everything
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped = allPeds[p];
        Point vNeu = resultAcc[p];
        Point dirNeu = resultDir[p];
        Point force = resultForce[p];
        // The final velocity considers acceleration
        // TODO: Pmax should be setted here to avoid cross wall and pedestrians
        Point FinalV = vNeu + force*deltaT;
        vNeu = FinalV;
        dirNeu=FinalV.Normalized();
        UpdatePed(ped, vNeu, dirNeu, deltaT, periodic);
    }

    // remove the pedestrians that have left the building
    for (unsigned int p = 0; p < pedsToRemove.size(); p++)
    {
        building->DeletePedestrian(pedsToRemove[p]);
    }
    pedsToRemove.clear();
}

/*----------Functions important----------*/
Point AVMModel::DesireDirection(Pedestrian* ped, Room* room) const
{

    Point target = this->GetDirection()->GetTarget(room, ped);
    Point desiredDirection;
    const Point pos = ped->GetPos();
    double dist = (target - pos).Norm();
    if (dist > J_EPS_GOAL)
    {
        Point NewE0;
        NewE0 = (target - pos).Normalized();
        ped->SetLastE0(NewE0);
    }
    desiredDirection = ped->GetLastE0();
    return desiredDirection;
}

Point AVMModel::ForceRepPedCSM(Pedestrian * ped1, Pedestrian * ped2, Building * building, int periodic) const
{
    Point FRep(0.0, 0.0);
    Point p1 = ped1->GetPos();
    Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();
    //Direction
    Point distp12 = p2 - p1;
    Point ep12 = distp12.Normalized();

    // Distance
    double Distance = distp12.Norm();
    double dist = Distance - (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());

    double R_ij = _aPed * exp((-dist) / _DPed);
    FRep = ep12 * (-R_ij);
    return FRep;
}

Point AVMModel::ForceRepPedGCVM(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const
{
    Point FRep(0.0, 0.0);
    Point p1 = ped1->GetPos();
    Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();

    //Direction
    Point distp12 = p2 - p1;
    Point ep12 = distp12.Normalized();

    // Distance
    double Distance = distp12.Norm();
    double dist = Distance - (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());

    Point e1 = ped1->GetMoveDirection();
    Room* room1 = building->GetRoom(ped1->GetRoomID());
    Point d1 = DesireDirection(ped1, room1);

    Point e2 = ped2->GetMoveDirection();
    Room* room2 = building->GetRoom(ped2->GetRoomID());
    Point d2 = DesireDirection(ped2, room2);

    double condition1 = d1.ScalarProduct(ep12);
    double condition2 = e1.ScalarProduct(ep12);

    if (condition1 >= 0 || condition2 >= 0)
    {
        Point infd = GetInfDirection(d1, ep12);
        double R_ij = _aPed * exp((-dist) / _DPed);
        FRep = infd * R_ij;
    }
    return FRep;
}


Point AVMModel::ForceRepPedAVM(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const
{
    Point FRep(0.0, 0.0);
    Point p1 = ped1->GetPos();
    Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();

    //Direction
    Point distp12 = p2 - p1;
    Point ep12 = distp12.Normalized();

    // Distance
    double Distance = distp12.Norm();
    double dist = Distance - (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());

    Point e1 = ped1->GetMoveDirection();
    Room* room1 = building->GetRoom(ped1->GetRoomID());
    Point d1 = DesireDirection(ped1, room1);

    Point e2 = ped2->GetMoveDirection();
    Room* room2 = building->GetRoom(ped2->GetRoomID());
    Point d2 = DesireDirection(ped2, room2);

    double condition1 = d1.ScalarProduct(ep12);
    double condition2 = e1.ScalarProduct(ep12);

    if (condition1 >= 0 || condition2 >= 0)
    {
        double S_Gap = (ped1->GetV() - ped2->GetV()).ScalarProduct(ep12);
        double Dis_Gap = S_Gap * GetAntiT();
        double R_dist = dist - Dis_Gap;
        R_dist = R_dist < 0 ? 0 : R_dist;

        double multi_d = d1.ScalarProduct(e2);
        double alpha = GetConstantAlpha() ? 1 : (1 + 0.5*(1 - multi_d));
        double R_ij = _aPed * exp((-R_dist) / _DPed) * alpha;

        Point newep12 = distp12 + ped2->GetV() * GetAntiT();
        Point infd = GetInfDirection(d1, newep12);
        FRep = infd * R_ij;
    }
    return FRep;
}

Point AVMModel::ForceConPed(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const
{
    Point FCon(0.0, 0.0);
    Point p1 = ped1->GetPos();
    Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();
    //Direction
    Point distp12 = p2 - p1;
    Point ep12 = distp12.Normalized();

    // Distance
    double Distance = distp12.Norm();
    double dist = Distance - (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());
    if (dist >= 0)
    {
        return FCon;
    }
    //TODO: the value of aForce and DForce should be set from inifile
    double aForce = 0.2;
    double DForce = 0.5;
    double R_ij = aForce * exp((-dist) / DForce);
    FCon = ep12 * (-R_ij);
    return FCon;
}

Point AVMModel::ForceRepRoom(Pedestrian* ped, SubRoom* subroom) const
{
    Point f(0., 0.);
    const Point& centroid = subroom->GetCentroid();
    bool inside = subroom->IsInSubRoom(centroid);
    //first the walls
    for (const auto & wall : subroom->GetAllWalls())
    {
        f += ForceRepWall(ped, wall, centroid, inside);
    }
    //then the obstacles
    for (const auto & obst : subroom->GetAllObstacles())
    {
        if (obst->Contains(ped->GetPos()))
        {
            Log->Write("ERROR:\t Agent [%d] is trapped in obstacle in room/subroom [%d/%d]",
                ped->GetID(), subroom->GetRoomID(), subroom->GetSubRoomID());
            exit(EXIT_FAILURE);
        }
        for (const auto & wall : obst->GetAllWalls())
        {
            f += ForceRepWall(ped, wall, centroid, inside);
        }
    }
    return f;
}

Point AVMModel::ForceRepWall(Pedestrian* ped, const Line& w, const Point& centroid, bool inside) const
{
    Point F_wrep = Point(0.0, 0.0);
    Point pt = w.ShortestPoint(ped->GetPos());
    Point dist = pt - ped->GetPos(); // ped ---> wall
    double Distance = dist.Norm(); //Distance between the center of pedestrian and walls 
    Point e_iw = Point(0.0, 0.0);
    if (Distance > J_EPS)
    {
        e_iw = dist.Normalized();
    }
    else
    {
        e_iw = (centroid - pt).Normalized();
        printf("ERROR: \tIn ID=%d, AGCVMModel::ForceRepWall() eiw can not be calculated!!!\n", ped->GetID());
        //exit(EXIT_FAILURE);
    }

    //rule: wall in behind has no influence 
    Point pt1 = w.GetPoint1();
    Point pt2 = w.GetPoint2();
    Point dist_pt1 = pt1 - ped->GetPos();
    Point dist_pt2 = pt2 - ped->GetPos();
    Point e_iw1 = dist_pt1.Normalized();
    Point e_iw2 = dist_pt2.Normalized();
    Point ei = ped->GetMoveDirection();
    //vision area----------------------------------
    Point e0 = ped->GetLastE0();
    double result_e01 = e0.ScalarProduct(e_iw1);
    double result_e02 = e0.ScalarProduct(e_iw2);
    double result_ei1 = ei.ScalarProduct(e_iw1);
    double result_ei2 = ei.ScalarProduct(e_iw2);
    //----------------------------------------------
    if (GetModel() == 0)//all
    {
        result_e01 = 1;
        result_e02 = 1;
        result_ei1 = 1;
        result_ei2 = 1;
    }
    if (result_e01 < 0 && result_e02 < 0 && result_ei1 < 0 && result_ei1 < 0)
    {
        return F_wrep;
    }

    // The distance between the real body(which is r2 > r1)
    double effdis = Distance - ped->GetEllipse().GetBmax(); //Using circle now.
    double R_iw = _aWall * exp((-effdis) / _DWall);
    if (GetModel() == 0)
    {
        Point inf_direction = GetInfDirection(e0, e_iw);
        F_wrep = inf_direction * R_iw;//new method
    }
    else
    {
        F_wrep = e_iw * R_iw*-1;//original method
    }
    return F_wrep;
}

Point AVMModel::ForceConRoom(Pedestrian* ped, SubRoom* subroom) const
{
    Point f(0., 0.);
    const Point& centroid = subroom->GetCentroid();
    bool inside = subroom->IsInSubRoom(centroid);
    //first the walls
    for (const auto& wall : subroom->GetAllWalls())
    {
        f += ForceConWall(ped, wall, centroid, inside);
    }
    //then the obstacles
    for (const auto& obst : subroom->GetAllObstacles())
    {
        if (obst->Contains(ped->GetPos()))
        {
            Log->Write("ERROR:\t Agent [%d] is trapped in obstacle in room/subroom [%d/%d]",
                ped->GetID(), subroom->GetRoomID(), subroom->GetSubRoomID());
            exit(EXIT_FAILURE);
        }
        for (const auto& wall : obst->GetAllWalls())
        {
            f += ForceConWall(ped, wall, centroid, inside);
        }
    }
    return f;
}

Point AVMModel::ForceConWall(Pedestrian* ped, const Line& w, const Point& centroid, bool inside) const
{
    Point F_wrep = Point(0.0, 0.0);
    Point pt = w.ShortestPoint(ped->GetPos());
    Point dist = pt - ped->GetPos(); // ped ---> wall
    double Distance = dist.Norm(); //Distance between the center of pedestrian and walls 
    Point e_iw = Point(0.0, 0.0);
    if (Distance > J_EPS)
    {
        e_iw = dist.Normalized();
    }
    else
    {
        e_iw = (centroid - pt).Normalized();
        printf("ERROR: \tIn ID=%d, AGCVMModel::ForceRepWall() eiw can not be calculated!!!\n", ped->GetID());
        //exit(EXIT_FAILURE);
    }

    // The distance between the real body(which is r2 > r1)
    double effdis = Distance - ped->GetEllipse().GetBmax(); //Using circle now.
    if (effdis>=0)
    {
        return F_wrep;
    }
    //TODO: the value of aForce and DForce should be set from inifile
    double aForce = 0.2;
    double DForce = 0.5;
    double R_iw = aForce * exp((-effdis) / DForce);
    F_wrep = e_iw * (-R_iw);
    return F_wrep;
}

my_pair AVMModel::GetSpacing(Pedestrian* ped1, Pedestrian* ped2, int periodic) const
{
    Point e1 = ped1->GetMoveDirection();
    Point p1 = ped1->GetPos();
    Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();
    Point distp12 = p2 - p1; //ped1 ---> ped2
    double Distance = distp12.Norm();
    Point ep12;
    if (Distance >= J_EPS)
    {
        ep12 = distp12.Normalized();
    }
    else
    {
        printf("ERROR: \tIn AGCVMModel::GetSpacing() ep12 can not be calculated!!!\n");
        exit(EXIT_FAILURE);
    }

    //Judge conllision
    double condition1 = e1.ScalarProduct(ep12); // < e_i , e_ij > should be positive
    double l = ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax();
    double condition2 = e1.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
    condition2 = (condition2 > 0) ? condition2 : -condition2; // abs
    if ((condition1 > 0) && (condition2 < l / Distance))
    {
        return  my_pair((Distance - l), ped2->GetID());
    }
    else
    {
        return  my_pair(FLT_MAX, -1);
    }
}

double AVMModel::GetSpacingRoom(Pedestrian* ped, SubRoom* subroom) const
{
    double spacing = FLT_MAX;
    //first the walls
    for (const auto & wall : subroom->GetAllWalls())
    {
        double distance = GetSpacingWall(ped, wall);
        spacing = spacing > distance ? distance : spacing;
    }
    //then the obstacles
    for (const auto & obst : subroom->GetAllObstacles())
    {
        for (const auto & wall : obst->GetAllWalls())
        {
            double distance = GetSpacingWall(ped, wall);
            spacing = spacing > distance ? distance : spacing;
        }
    }
    //and finally the closed doors
    for (const auto & goal : subroom->GetAllTransitions())
    {
        if (!goal->IsOpen())
        {
            double distance = GetSpacingWall(ped, *(static_cast<Line*>(goal)));
            spacing = spacing > distance ? distance : spacing;
        }
        /*
        //door is open, bur not my door
        int uid1 = goal->GetUniqueID();
        int uid2 = ped->GetExitIndex();
        if ((uid1 != uid2) && (goal->IsOpen() == true))
        {
            double distance = GetSpacingWall(ped, *(static_cast<Line*>(goal)));
            spacing = spacing > distance ? distance : spacing;
        }
        */
    }
    return spacing;
}

double AVMModel::GetSpacingWall(Pedestrian* ped, const Line& l) const
{
    double spacing = FLT_MAX;
    Point pp = ped->GetPos();
    Point pt = l.ShortestPoint(ped->GetPos());
    Point p1 = l.GetPoint1();
    Point p2 = l.GetPoint2();
    Point dist = pt - pp;
    Point ei_vertical;
    Point ei = ped->GetMoveDirection();
    ei_vertical._x = -ei.Normalized()._y;
    ei_vertical._y = ei.Normalized()._x;
    //--------------------------------------------------------
    double b = ped->GetEllipse().GetBmax();
    //--------------------------------------------------------
    Point A1 = pp + ei_vertical * b;
    Point A2 = pp - ei_vertical * b;
    Point p1_A1 = p1 - A1;
    Point p2_A1 = p2 - A1;
    Point p12 = p1 - p2;
    double A1_result1 = ei.CrossProduct(p1_A1);
    double A1_result2 = ei.CrossProduct(p2_A1);
    double result3 = ei.ScalarProduct(dist);
    Point p1_A2 = p1 - A2;
    Point p2_A2 = p2 - A2;
    double A2_result1 = ei.CrossProduct(p1_A2);
    double A2_result2 = ei.CrossProduct(p2_A2);
    if (result3 <= 0)
    {
        return spacing;
    }
    if (A1_result1*A1_result2 > 0 && A2_result1*A2_result2 > 0 && A1_result1*A2_result1 > 0 && A1_result2*A2_result2 > 0)
    {
        return spacing;
    }
    double effdis = dist.Norm() - b;
    double cosangle = dist.ScalarProduct(ei) / (dist.Norm()*ei.Norm());
    if (cosangle < 0.00001)
    {
        return spacing;
    }
    spacing = effdis / cosangle;
    return spacing;
}

double AVMModel::OptimalSpeed(Pedestrian* ped, double spacing) const
{
    double v0 = ped->GetV0Norm();
    v0 = v0 < 0.1 ? 0.1 : v0; //To avoid pedestrians who's desired speed is zero
    double T = GetTs();
    double speed = (spacing) / T;
    speed = (speed > 0) ? speed : 0;
    speed = (speed < v0) ? speed : v0;
    return speed;
}


double AVMModel::PushSpeed(Pedestrian* ped, double spacing) const
{
    int plevel = 0;
    // TODO: the plevel of pedestrians should be generated by the RF classifier
    // Plevel=RFClassifier(pedestrian,environment,classifier)
    int random = rand() % 10000;
    if (random > 5000)
    {
        plevel = 1;
    }
    //
    double v0 = ped->GetV0Norm();
    v0 = v0 < 0.1 ? 0.1 : v0; //To avoid pedestrians who's desired speed is zero
    // TODO: the value of d and T are set from the inifile
    double d = 0;
    double T = 0.5;
    if (plevel == 1)
    {
        d = 0.2;
        T = 0.1;
    }
    double speed = (spacing+d) / T;
    speed = (speed > 0) ? speed : 0;
    speed = (speed < v0) ? speed : v0;
    return speed;
}

/*----------Functions helpful----------*/
Point AVMModel::GetPosPeriodic(Pedestrian* ped1, Pedestrian* ped2) const
{
    double x1 = ped1->GetPos()._x;
    double y1 = ped1->GetPos()._y;
    double x2 = ped2->GetPos()._x;
    double y2 = ped2->GetPos()._y;
    double xL = GetLeftBoundary();
    double xR = GetRightBoundary();
    double yU = GetUpBoundary();
    double yD = GetDownBoundary();
    double cutoff = GetCutoff();
    if ((xR - x1) + (x2 - xL) <= cutoff) {
        double x2_periodic = x2 + xR - xL;
        return Point(x2_periodic, y2);
    }
    if ((x1 - xL) + (xR - x2) <= cutoff) {
        double x2_periodic = xL - xR + x2;
        return Point(x2_periodic, y2);
    }
    if ((y1 - yD) + (yU - y2) <= cutoff) {
        double y2_periodic = yD - yU + y2;
        return Point(x2, y2_periodic);
    }
    if ((y2 - yD) + (yU - y1) <= cutoff) {
        double y2_periodic = yU + y2 - yD;
        return Point(x2, y2_periodic);
    }
    return Point(x2, y2);
}

Point AVMModel::GetInfDirection(Point e0, Point ep12) const
{
    Point InfDirection = Point(-e0._y, e0._x).Normalized();
    double result = InfDirection.ScalarProduct(ep12);
    Point zero = Point(0, 0);
    if (fabs(result) < J_EPS)//if neighbour is in front or behind pedestrian
    {
        int random = rand() % 10000;//choose one direciton bu random
        if (random < 5000)//when random is larger than 50, influence's direction is right, otherwise is left
        {
            InfDirection = zero - InfDirection;
        }
    }
    else if (result > 0)
    {
        InfDirection = zero - InfDirection;
    }
    return InfDirection;
}

void AVMModel::UpdatePed(Pedestrian* ped, Point speed, Point direction, double deltaT, int periodic)
{
    Point e0 = ped->GetLastE0();
    Point MD = direction;
    Point MS = speed;

    Point pos_neu = ped->GetPos() + MS * deltaT;
    ped->SetMoveDirection(MD);
    ped->SetV(MS);

    if (periodic)
    {
        double xLeft_gcvm = GetLeftBoundary();
        double xRight_gcvm = GetRightBoundary();
        double yUp_gcvm = GetUpBoundary();
        double yDown_gcvm = GetDownBoundary();
        if ((ped->GetPos()._x < xRight_gcvm) && (pos_neu._x >= xRight_gcvm))
        {
            ped->SetmoveManually(true);
            ped->SetPos(Point(pos_neu._x - (xRight_gcvm - xLeft_gcvm), pos_neu._y));
        }
        else if ((ped->GetPos()._x > xLeft_gcvm) && (pos_neu._x <= xLeft_gcvm))
        {
            ped->SetmoveManually(true);
            ped->SetPos(Point(pos_neu._x + (xRight_gcvm - xLeft_gcvm), pos_neu._y));
        }
        else if ((ped->GetPos()._y < yUp_gcvm) && (pos_neu._y >= yUp_gcvm))
        {
            ped->SetmoveManually(true);
            ped->SetPos(Point(pos_neu._x, pos_neu._y - (yUp_gcvm - yDown_gcvm)));
        }
        else if ((ped->GetPos()._y > yDown_gcvm) && (pos_neu._y <= yDown_gcvm))
        {
            ped->SetmoveManually(true);
            ped->SetPos(Point(pos_neu._x, pos_neu._y + (yUp_gcvm - yDown_gcvm)));
        }
        else
        {
            ped->SetPos(pos_neu);
        }
    }
    else
    {
        ped->SetPos(pos_neu);
    }
}

