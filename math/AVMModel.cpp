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
    double lb, double rb, double ub, double db, double co,
    double apush, double Dpush, double Tpush, double Spush, double Snorm)
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

    // PVM Parameter
    _aPush = apush;
    _DPush = Dpush;
    _TPush = Tpush;
    _Spush = Spush;
    _Snorm = Snorm;
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
        // check if pedestrian cross wall----------------------------------------------------------------------------------
        Point posNew = ped->GetPos();
        Room* room1= building->GetRoom(ped->GetRoomID());
        SubRoom* subroom1 = room1->GetSubRoom(ped->GetSubRoomID());
        bool outside = 1;
        for (auto it = room1->GetAllSubRooms().begin(); it != room1->GetAllSubRooms().end(); ++it)
        {
            if (it->second->IsInSubRoom(posNew))
            {
                outside = 0;
                break;
            }
        }
        if (outside == 1 && subroom1->GetAllTransitions().size() == 0)
        {
            printf("ERROR: ped %d cross walls.\n", ped->GetID());
        }
        //----------------------------------------------------------------------------------------------------------------------------
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
    // the value of aForce and DForce should be set from inifile
    double aForce = GetaPush();
    double DForce = GetDPush();
    // printf("Test: aForce is %0.2f, DForce is %0.2f.\n", aForce, DForce);
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
    // the value of aForce and DForce should be set from inifile
    double aForce = GetaPush();
    double DForce = GetDPush();
   //  printf("Test: aForce is %0.2f, DForce is %0.2f.\n", aForce, DForce);
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
    double meanPlevel = ped->GetP0();
    double plevel = ped->GetPlevel();
    // The plevel of pedestrians should be generated by the RF classifier
    // Here the RF classifier code is generated by m2cgen
    // The input should be the surrouding neighbors
    double input[9] = { 2.434, 0.092, 0.043, 0.206, 0.19, 3.588, 2.42, 2.167, 2.0 };
    double output=-1;
    double *po=&output;
    score(input, po);
    // pushing: score<0.5, nonpushing: score>0.5
    if (output < 0.5)
    {
        plevel = 1;
        ped->SetPlevel(plevel);
    }
    //-------------------------------------------------------------------------------------------
    double v0 = ped->GetV0Norm();
    v0 = v0 < 0.1 ? 0.1 : v0; //To avoid pedestrians who's desired speed is zero
    // the value of d and T are set from the inifile
    double extraSpace = GetSNorm();
    double T = GetTs();
    if (plevel == 1)
    {
        extraSpace = GetSPush();
        T = GetTPush();
    }
    //the value of T and extraSpace should be normal distribution
    static std::default_random_engine gen;
    std::normal_distribution<double> spaceDist = std::normal_distribution<double>(0, 0.1);
    std::normal_distribution<double> TDist = std::normal_distribution<double>(0, 0.1);
    extraSpace = extraSpace + spaceDist(gen);
    T = T + TDist(gen);
    //printf("id=%d, p0=%f, plevel=%f, extraSpace=%f, T=%f\n", ped->GetID(), meanPlevel, plevel, extraSpace, T);
    double speed = (spacing+extraSpace) / T;
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

void AVMModel::add_vectors(double* v1, double* v2, int size, double* result) const
{
    for (int i = 0; i < size; ++i)
        result[i] = v1[i] + v2[i];
}
void AVMModel::mul_vector_number(double* v1, double num, int size, double* result) const
{
    for (int i = 0; i < size; ++i)
        result[i] = v1[i] * num;
}
void AVMModel::score(double* input, double* output) const
{
    double var0[2];
    double var1[2];
    double var2[2];
    double var3[2];
    double var4[2];
    double var5[2];
    double var6[2];
    double var7[2];
    double var8[2];
    double var9[2];
    double var10[2];
    double var11[2];
    double var12[2];
    double var13[2];
    double var14[2];
    double var15[2];
    double var16[2];
    double var17[2];
    double var18[2];
    double var19[2];
    double var20[2];
    double var21[2];
    if (input[0] <= 2.4404999017715454) {
        if (input[2] <= 0.06149999983608723) {
            if (input[0] <= 2.149500012397766) {
                double tempArray[2] = { 0.9820643781921844, 0.017935621807815653 };
                memcpy(var21, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.6955100351748397, 0.30448996482516033 };
                memcpy(var21, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.0649999380111694) {
                double tempArray[2] = { 0.9967614300398485, 0.003238569960151509 };
                memcpy(var21, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.7909264699025504, 0.20907353009744953 };
                memcpy(var21, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[0] <= 2.7415000200271606) {
            if (input[1] <= 0.02049999963492155) {
                double tempArray[2] = { 0.3837880674742443, 0.6162119325257557 };
                memcpy(var21, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.5444735120994114, 0.45552648790058864 };
                memcpy(var21, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[2] <= 0.07349999994039536) {
                double tempArray[2] = { 0.04654462497217895, 0.9534553750278211 };
                memcpy(var21, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.12502986382097636, 0.8749701361790236 };
                memcpy(var21, tempArray, 2 * sizeof(double));
            }
        }
    }
    double var22[2];
    if (input[0] <= 2.438499927520752) {
        if (input[1] <= 0.021499999798834324) {
            if (input[2] <= 0.11649999767541885) {
                double tempArray[2] = { 0.8251261352169525, 0.17487386478304742 };
                memcpy(var22, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9316432775011317, 0.06835672249886826 };
                memcpy(var22, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.070499897003174) {
                double tempArray[2] = { 0.9966471812640999, 0.0033528187359000883 };
                memcpy(var22, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8079449894884373, 0.1920550105115627 };
                memcpy(var22, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[8] <= 2.072000026702881) {
            if (input[7] <= 2.598499894142151) {
                double tempArray[2] = { 0.39850136239782014, 0.6014986376021798 };
                memcpy(var22, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.13147502903600464, 0.8685249709639954 };
                memcpy(var22, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.8289999961853027) {
                double tempArray[2] = { 0.36508972267536705, 0.6349102773246329 };
                memcpy(var22, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.04207065171570694, 0.9579293482842931 };
                memcpy(var22, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var21, var22, 2, var20);
    double var23[2];
    if (input[0] <= 2.4404999017715454) {
        if (input[7] <= 2.4954999685287476) {
            if (input[1] <= 0.060499999672174454) {
                double tempArray[2] = { 0.90564116188452, 0.09435883811547999 };
                memcpy(var23, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9701938084342279, 0.029806191565772126 };
                memcpy(var23, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[1] <= 0.025500000454485416) {
                double tempArray[2] = { 0.8086973180076629, 0.19130268199233716 };
                memcpy(var23, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9225866050808315, 0.07741339491916858 };
                memcpy(var23, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[6] <= 3.448499917984009) {
            if (input[0] <= 2.9809999465942383) {
                double tempArray[2] = { 0.495396044783928, 0.504603955216072 };
                memcpy(var23, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.05414118954974986, 0.9458588104502501 };
                memcpy(var23, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[8] <= 1.9419999718666077) {
                double tempArray[2] = { 0.31069223269418267, 0.6893077673058173 };
                memcpy(var23, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.11969464999574722, 0.8803053500042528 };
                memcpy(var23, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var20, var23, 2, var19);
    double var24[2];
    if (input[5] <= 4.888499975204468) {
        if (input[0] <= 2.3554999828338623) {
            if (input[7] <= 2.5325000286102295) {
                double tempArray[2] = { 0.9517899651020361, 0.0482100348979639 };
                memcpy(var24, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8500606428138265, 0.14993935718617343 };
                memcpy(var24, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[8] <= 2.1134999990463257) {
                double tempArray[2] = { 0.4332552693208431, 0.5667447306791569 };
                memcpy(var24, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.23803598774885146, 0.7619640122511485 };
                memcpy(var24, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[0] <= 2.5689998865127563) {
            if (input[0] <= 2.277999997138977) {
                double tempArray[2] = { 0.9571764055808814, 0.04282359441911866 };
                memcpy(var24, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.586925828313253, 0.413074171686747 };
                memcpy(var24, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.853999972343445) {
                double tempArray[2] = { 0.30638051044083525, 0.6936194895591647 };
                memcpy(var24, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.02346882131622772, 0.9765311786837723 };
                memcpy(var24, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var19, var24, 2, var18);
    double var25[2];
    if (input[6] <= 3.549499988555908) {
        if (input[5] <= 6.019500017166138) {
            if (input[7] <= 2.4774999618530273) {
                double tempArray[2] = { 0.8907649216931635, 0.10923507830683656 };
                memcpy(var25, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.6936197317668708, 0.30638026823312925 };
                memcpy(var25, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[6] <= 2.612499952316284) {
                double tempArray[2] = { 0.7083251714005877, 0.29167482859941235 };
                memcpy(var25, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.49050632911392406, 0.509493670886076 };
                memcpy(var25, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[0] <= 2.385499954223633) {
            if (input[8] <= 2.3174999952316284) {
                double tempArray[2] = { 0.9238722219767596, 0.0761277780232404 };
                memcpy(var25, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.848224709841357, 0.15177529015864294 };
                memcpy(var25, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.7730000019073486) {
                double tempArray[2] = { 0.35219319335363575, 0.6478068066463643 };
                memcpy(var25, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.017969765157037364, 0.9820302348429626 };
                memcpy(var25, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var18, var25, 2, var17);
    double var26[2];
    if (input[6] <= 3.6085000038146973) {
        if (input[0] <= 2.503499984741211) {
            if (input[1] <= 0.025500000454485416) {
                double tempArray[2] = { 0.8612535230170447, 0.1387464769829553 };
                memcpy(var26, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9556618051733626, 0.04433819482663731 };
                memcpy(var26, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[7] <= 2.6200000047683716) {
                double tempArray[2] = { 0.3532453143754039, 0.646754685624596 };
                memcpy(var26, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.1336494295451178, 0.8663505704548822 };
                memcpy(var26, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[5] <= 3.4299999475479126) {
            if (input[5] <= 1.4305000305175781) {
                double tempArray[2] = { 0.44675925925925924, 0.5532407407407407 };
                memcpy(var26, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8083629893238434, 0.19163701067615657 };
                memcpy(var26, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[8] <= 2.013000011444092) {
                double tempArray[2] = { 0.6547239576692592, 0.3452760423307408 };
                memcpy(var26, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.37182625389968804, 0.628173746100312 };
                memcpy(var26, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var17, var26, 2, var16);
    double var27[2];
    if (input[6] <= 3.54449999332428) {
        if (input[0] <= 2.503499984741211) {
            if (input[7] <= 2.4854999780654907) {
                double tempArray[2] = { 0.9467159286271076, 0.05328407137289239 };
                memcpy(var27, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8634917341080446, 0.13650826589195544 };
                memcpy(var27, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[1] <= 0.13749999552965164) {
                double tempArray[2] = { 0.2714176726652605, 0.7285823273347395 };
                memcpy(var27, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.48583162217659137, 0.5141683778234086 };
                memcpy(var27, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[8] <= 2.0234999656677246) {
            if (input[1] <= 0.015500000212341547) {
                double tempArray[2] = { 0.5693946310183675, 0.43060536898163243 };
                memcpy(var27, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.7860749973138498, 0.21392500268615022 };
                memcpy(var27, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.4609999656677246) {
                double tempArray[2] = { 0.82843370222011, 0.17156629777989002 };
                memcpy(var27, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.08792434652536572, 0.9120756534746343 };
                memcpy(var27, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var16, var27, 2, var15);
    double var28[2];
    if (input[8] <= 2.188499927520752) {
        if (input[0] <= 2.3615000247955322) {
            if (input[0] <= 2.1200000047683716) {
                double tempArray[2] = { 0.9925235028499518, 0.007476497150048116 };
                memcpy(var28, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.7845370216995197, 0.21546297830048036 };
                memcpy(var28, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[7] <= 2.6080000400543213) {
                double tempArray[2] = { 0.37926523492859476, 0.6207347650714052 };
                memcpy(var28, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.2062552831783601, 0.7937447168216399 };
                memcpy(var28, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[0] <= 2.4609999656677246) {
            if (input[0] <= 2.1320000886917114) {
                double tempArray[2] = { 0.9792362645834404, 0.020763735416559594 };
                memcpy(var28, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.6734477761676513, 0.32655222383234866 };
                memcpy(var28, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.8289999961853027) {
                double tempArray[2] = { 0.3369087898639407, 0.6630912101360593 };
                memcpy(var28, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.041770346040277365, 0.9582296539597226 };
                memcpy(var28, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var15, var28, 2, var14);
    double var29[2];
    if (input[0] <= 2.430500030517578) {
        if (input[0] <= 2.1290000677108765) {
            if (input[0] <= 2.037000060081482) {
                double tempArray[2] = { 0.9976754772535764, 0.002324522746423625 };
                memcpy(var29, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9200034668053388, 0.07999653319466112 };
                memcpy(var29, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[1] <= -0.019499999471008778) {
                double tempArray[2] = { 0.6397853361375009, 0.3602146638624991 };
                memcpy(var29, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.7811963464047799, 0.2188036535952201 };
                memcpy(var29, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[8] <= 2.072000026702881) {
            if (input[8] <= 1.9419999718666077) {
                double tempArray[2] = { 0.4391833722789325, 0.5608166277210676 };
                memcpy(var29, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.29992273653157264, 0.7000772634684274 };
                memcpy(var29, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.8289999961853027) {
                double tempArray[2] = { 0.36785680991889624, 0.6321431900811038 };
                memcpy(var29, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.04226696997781335, 0.9577330300221867 };
                memcpy(var29, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var14, var29, 2, var13);
    double var30[2];
    if (input[0] <= 2.430500030517578) {
        if (input[0] <= 2.1290000677108765) {
            if (input[6] <= 3.8244999647140503) {
                double tempArray[2] = { 0.9936065591180875, 0.006393440881912509 };
                memcpy(var30, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9741464867010718, 0.025853513298928148 };
                memcpy(var30, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[1] <= 0.022499999962747097) {
                double tempArray[2] = { 0.6771554162858511, 0.3228445837141489 };
                memcpy(var30, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.7947425102800079, 0.20525748971999216 };
                memcpy(var30, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[8] <= 2.072000026702881) {
            if (input[6] <= 3.4614999294281006) {
                double tempArray[2] = { 0.44712182061579653, 0.5528781793842035 };
                memcpy(var30, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.22831727205337288, 0.7716827279466272 };
                memcpy(var30, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[1] <= 0.0925000011920929) {
                double tempArray[2] = { 0.1417514778475227, 0.8582485221524773 };
                memcpy(var30, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.29502487562189056, 0.7049751243781095 };
                memcpy(var30, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var13, var30, 2, var12);
    double var31[2];
    if (input[1] <= 0.039499999955296516) {
        if (input[0] <= 2.393499970436096) {
            if (input[0] <= 2.1290000677108765) {
                double tempArray[2] = { 0.9820884828945011, 0.017911517105498837 };
                memcpy(var31, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.7036364370182023, 0.2963635629817976 };
                memcpy(var31, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[2] <= 0.1314999982714653) {
                double tempArray[2] = { 0.15898594276264041, 0.8410140572373596 };
                memcpy(var31, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.31569365214225026, 0.6843063478577497 };
                memcpy(var31, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[6] <= 3.875499963760376) {
            if (input[8] <= 2.2885000705718994) {
                double tempArray[2] = { 0.9091898140278108, 0.09081018597218916 };
                memcpy(var31, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.6276035432128322, 0.3723964567871678 };
                memcpy(var31, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.4644999504089355) {
                double tempArray[2] = { 0.9095064838137886, 0.0904935161862114 };
                memcpy(var31, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.1061109061846797, 0.8938890938153203 };
                memcpy(var31, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var12, var31, 2, var11);
    double var32[2];
    if (input[1] <= 0.025500000454485416) {
        if (input[0] <= 2.393499970436096) {
            if (input[2] <= 0.10249999910593033) {
                double tempArray[2] = { 0.8364963881624938, 0.16350361183750614 };
                memcpy(var32, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9270798058471611, 0.07292019415283892 };
                memcpy(var32, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[6] <= 3.1614999771118164) {
                double tempArray[2] = { 0.32024550952252007, 0.6797544904774799 };
                memcpy(var32, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.1359332948865821, 0.8640667051134179 };
                memcpy(var32, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[7] <= 2.489500045776367) {
            if (input[6] <= 3.862499952316284) {
                double tempArray[2] = { 0.9099225897255454, 0.09007741027445461 };
                memcpy(var32, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.6652804642166344, 0.33471953578336555 };
                memcpy(var32, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[8] <= 2.0304999351501465) {
                double tempArray[2] = { 0.8620161290322581, 0.13798387096774192 };
                memcpy(var32, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.5363481749596616, 0.4636518250403384 };
                memcpy(var32, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var11, var32, 2, var10);
    double var33[2];
    if (input[6] <= 3.551500082015991) {
        if (input[7] <= 2.3914999961853027) {
            if (input[1] <= 0.0025000000605359674) {
                double tempArray[2] = { 0.7801990852838311, 0.21980091471616894 };
                memcpy(var33, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9125629690552072, 0.08743703094479284 };
                memcpy(var33, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[4] <= 0.7004999816417694) {
                double tempArray[2] = { 0.7014114879435405, 0.29858851205645953 };
                memcpy(var33, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.35262943334692215, 0.6473705666530779 };
                memcpy(var33, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[0] <= 2.385499954223633) {
            if (input[8] <= 2.3174999952316284) {
                double tempArray[2] = { 0.9201303275801338, 0.07986967241986614 };
                memcpy(var33, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8521661147902869, 0.14783388520971302 };
                memcpy(var33, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[8] <= 1.9419999718666077) {
                double tempArray[2] = { 0.33482913358646876, 0.6651708664135313 };
                memcpy(var33, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.12560437696157434, 0.8743956230384257 };
                memcpy(var33, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var10, var33, 2, var9);
    double var34[2];
    if (input[0] <= 2.4469999074935913) {
        if (input[1] <= 0.016500000841915607) {
            if (input[8] <= 1.8149999976158142) {
                double tempArray[2] = { 0.9596559673701011, 0.04034403262989892 };
                memcpy(var34, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8219700651180872, 0.1780299348819127 };
                memcpy(var34, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.070499897003174) {
                double tempArray[2] = { 0.9964360254360376, 0.0035639745639623587 };
                memcpy(var34, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8036902221664366, 0.19630977783356346 };
                memcpy(var34, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[6] <= 3.436500072479248) {
            if (input[7] <= 2.634999990463257) {
                double tempArray[2] = { 0.40487264673311185, 0.5951273532668881 };
                memcpy(var34, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.1617867867867868, 0.8382132132132132 };
                memcpy(var34, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.7730000019073486) {
                double tempArray[2] = { 0.3452929741738406, 0.6547070258261594 };
                memcpy(var34, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.01954417164881187, 0.9804558283511882 };
                memcpy(var34, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var9, var34, 2, var8);
    double var35[2];
    if (input[8] <= 2.188499927520752) {
        if (input[0] <= 2.3615000247955322) {
            if (input[6] <= 1.6665000319480896) {
                double tempArray[2] = { 0.9899439019322668, 0.010056098067733223 };
                memcpy(var35, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9409528046093882, 0.05904719539061176 };
                memcpy(var35, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[8] <= 2.072000026702881) {
                double tempArray[2] = { 0.3933800863131936, 0.6066199136868065 };
                memcpy(var35, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.22016548665017488, 0.7798345133498251 };
                memcpy(var35, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[0] <= 2.4609999656677246) {
            if (input[0] <= 2.1320000886917114) {
                double tempArray[2] = { 0.9811899407074218, 0.018810059292578205 };
                memcpy(var35, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.6714213226588487, 0.32857867734115126 };
                memcpy(var35, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[6] <= 3.1424999237060547) {
                double tempArray[2] = { 0.2762792762792763, 0.7237207237207237 };
                memcpy(var35, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.09827173317175913, 0.9017282668282408 };
                memcpy(var35, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var8, var35, 2, var7);
    double var36[2];
    if (input[7] <= 2.367500066757202) {
        if (input[1] <= 0.0005000000237487257) {
            if (input[0] <= 2.503499984741211) {
                double tempArray[2] = { 0.8836240855374227, 0.11637591446257738 };
                memcpy(var36, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.208778840742825, 0.7912211592571751 };
                memcpy(var36, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[6] <= 3.8545000553131104) {
                double tempArray[2] = { 0.914172058712455, 0.085827941287545 };
                memcpy(var36, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.6571834992887624, 0.34281650071123754 };
                memcpy(var36, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[8] <= 2.022499918937683) {
            if (input[0] <= 2.2419999837875366) {
                double tempArray[2] = { 0.9739815404201145, 0.02601845957988542 };
                memcpy(var36, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.3946357266631508, 0.6053642733368493 };
                memcpy(var36, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.4609999656677246) {
                double tempArray[2] = { 0.8187876951887858, 0.18121230481121423 };
                memcpy(var36, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.13122992549216683, 0.8687700745078332 };
                memcpy(var36, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var7, var36, 2, var6);
    double var37[2];
    if (input[0] <= 2.430500030517578) {
        if (input[0] <= 2.1290000677108765) {
            if (input[0] <= 2.0554999113082886) {
                double tempArray[2] = { 0.996887301146431, 0.0031126988535690684 };
                memcpy(var37, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9177518848526388, 0.0822481151473612 };
                memcpy(var37, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[1] <= -0.02350000012665987) {
                double tempArray[2] = { 0.6399515188243315, 0.3600484811756685 };
                memcpy(var37, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.7771020624008461, 0.22289793759915388 };
                memcpy(var37, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[6] <= 3.4054999351501465) {
            if (input[0] <= 2.9809999465942383) {
                double tempArray[2] = { 0.4933692608442965, 0.5066307391557034 };
                memcpy(var37, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.06227531021686188, 0.9377246897831382 };
                memcpy(var37, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[6] <= 4.680500030517578) {
                double tempArray[2] = { 0.1735899039699768, 0.8264100960300231 };
                memcpy(var37, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.10261155665227226, 0.8973884433477277 };
                memcpy(var37, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var6, var37, 2, var5);
    double var38[2];
    if (input[8] <= 2.1605000495910645) {
        if (input[0] <= 2.3615000247955322) {
            if (input[1] <= 0.017500000074505806) {
                double tempArray[2] = { 0.9043957494554815, 0.09560425054451852 };
                memcpy(var38, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.970790562371091, 0.029209437628909 };
                memcpy(var38, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[1] <= 0.04050000011920929) {
                double tempArray[2] = { 0.2930177306524704, 0.7069822693475296 };
                memcpy(var38, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.44404578848453785, 0.5559542115154622 };
                memcpy(var38, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[0] <= 2.4609999656677246) {
            if (input[2] <= 0.04649999924004078) {
                double tempArray[2] = { 0.8063445910683665, 0.1936554089316335 };
                memcpy(var38, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8935465002924546, 0.10645349970754533 };
                memcpy(var38, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[2] <= 0.1314999982714653) {
                double tempArray[2] = { 0.11262227989618687, 0.8873777201038131 };
                memcpy(var38, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.29617732091230325, 0.7038226790876968 };
                memcpy(var38, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var5, var38, 2, var4);
    double var39[2];
    if (input[6] <= 3.551500082015991) {
        if (input[0] <= 2.503499984741211) {
            if (input[6] <= 1.6704999804496765) {
                double tempArray[2] = { 0.9888302695573324, 0.011169730442667546 };
                memcpy(var39, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9135485848574639, 0.08645141514253601 };
                memcpy(var39, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.9809999465942383) {
                double tempArray[2] = { 0.45405436469168436, 0.5459456353083156 };
                memcpy(var39, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.0537666467899026, 0.9462333532100974 };
                memcpy(var39, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[0] <= 2.385499954223633) {
            if (input[0] <= 2.177000045776367) {
                double tempArray[2] = { 0.9616316919454256, 0.03836830805457442 };
                memcpy(var39, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.6599287929252325, 0.3400712070747674 };
                memcpy(var39, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[6] <= 4.681500196456909) {
                double tempArray[2] = { 0.16989341742844596, 0.830106582571554 };
                memcpy(var39, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.10954446854663774, 0.8904555314533622 };
                memcpy(var39, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var4, var39, 2, var3);
    double var40[2];
    if (input[0] <= 2.430500030517578) {
        if (input[1] <= 0.025500000454485416) {
            if (input[1] <= -0.04150000028312206) {
                double tempArray[2] = { 0.8176590613986273, 0.18234093860137265 };
                memcpy(var40, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.8828523403891627, 0.1171476596108373 };
                memcpy(var40, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[7] <= 2.9079999923706055) {
                double tempArray[2] = { 0.9572261177485529, 0.04277388225144703 };
                memcpy(var40, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.7950401167031363, 0.2049598832968636 };
                memcpy(var40, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[8] <= 2.072000026702881) {
            if (input[6] <= 3.4574999809265137) {
                double tempArray[2] = { 0.427602523659306, 0.572397476340694 };
                memcpy(var40, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.22694122694122695, 0.773058773058773 };
                memcpy(var40, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[0] <= 2.8289999961853027) {
                double tempArray[2] = { 0.3699687557434295, 0.6300312442565705 };
                memcpy(var40, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.04299264747151597, 0.957007352528484 };
                memcpy(var40, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var3, var40, 2, var2);
    double var41[2];
    if (input[0] <= 2.438499927520752) {
        if (input[7] <= 2.4954999685287476) {
            if (input[0] <= 2.070499897003174) {
                double tempArray[2] = { 0.9949668240779485, 0.005033175922051527 };
                memcpy(var41, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.792615925573218, 0.20738407442678192 };
                memcpy(var41, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[1] <= -0.018499999307096004) {
                double tempArray[2] = { 0.7805842581961985, 0.2194157418038015 };
                memcpy(var41, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.9000285904888974, 0.09997140951110264 };
                memcpy(var41, tempArray, 2 * sizeof(double));
            }
        }
    }
    else {
        if (input[8] <= 2.072000026702881) {
            if (input[0] <= 2.7545000314712524) {
                double tempArray[2] = { 0.5369258589511754, 0.4630741410488246 };
                memcpy(var41, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.05043227665706052, 0.9495677233429395 };
                memcpy(var41, tempArray, 2 * sizeof(double));
            }
        }
        else {
            if (input[6] <= 3.2325000762939453) {
                double tempArray[2] = { 0.3138094867201495, 0.6861905132798505 };
                memcpy(var41, tempArray, 2 * sizeof(double));
            }
            else {
                double tempArray[2] = { 0.12103386404415271, 0.8789661359558473 };
                memcpy(var41, tempArray, 2 * sizeof(double));
            }
        }
    }
    add_vectors(var2, var41, 2, var1);
    mul_vector_number(var1, 0.047619047619047616, 2, var0);
    memcpy(output, var0, 2 * sizeof(double));
}

