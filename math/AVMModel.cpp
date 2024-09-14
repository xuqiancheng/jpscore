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

    // Step 1: update the density of agents based on the number of pedestrians in the surrounding square area （2*2） 
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped1 = allPeds[p];
        Point p1 = ped1->GetPos();
        Room* room1 = building->GetRoom(ped1->GetRoomID());
        SubRoom* subroom1 = room1->GetSubRoom(ped1->GetSubRoomID());
        // Calculate the area of the measurement ares (square area 2*2)
        double radius = 1;
        // the distance to wall on each direction
        vector<Point> directions = { Point(1,0),Point(0,1),Point(-1,0),Point(0,-1) };
        vector<double> spaces;
        for (int j = 0; j < 4; j++)
        {
            Point direction = directions[j];
            double space =radius;
            for (const auto& wall : subroom1->GetAllWalls())
            {
                double distance = GetSpacingWallDirection(ped1, wall, direction) + ped1->GetEllipse().GetBmax();
                space = distance < space ? distance : space;
            }
            for (const auto& goal : subroom1->GetAllTransitions())
            {
                double distance = GetSpacingWallDirection(ped1, *(static_cast<Line*>(goal)), direction) + ped1->GetEllipse().GetBmax();
                space = distance < space ? distance : space;
            }
            for (const auto& goal : subroom1->GetAllCrossings())
            {
                double distance = GetSpacingWallDirection(ped1, *(static_cast<Line*>(goal)), direction) + ped1->GetEllipse().GetBmax();
                space = distance < space ? distance : space;
            }
            spaces.push_back(space);
        }
        double area = (spaces[0] + spaces[2]) * (spaces[1] + spaces[3]);
        // printf("pos=(%f,%f)\n", p1._x, p1._y);
        // printf("spaces=(%f,%f,%f,%f)\n", spaces[0], spaces[1], spaces[2], spaces[3]);
        // printf("area=%f\n", area);
        // calculate the numebr of agents in the square region
        vector<Pedestrian*> neighbours;
        building->GetGrid()->GetNeighbourhood(ped1, neighbours);
        int size = (int)neighbours.size();
        int number = 1;
        for (int i = 0; i < size; i++)
        {
            Pedestrian* ped2 = neighbours[i];
            Point p2 = ped2->GetPos();
            double dis12 = (p2 - p1).Norm();
            Point ep12 = (p2 - p1).Normalized();
            for (int j = 0; j < 4; j++)
            {
                Point direction = directions[j];
                if (dis12 * ep12.ScalarProduct(direction) > spaces[j])
                {
                    //printf("ped2 (%f,%f) not in the region\n", p2._x, p2._y);
                    break;
                }
                if (j == 3)
                {
                   //printf("ped2 (%f,%f) in the region\n", p2._x, p2._y);
                    number++;
                }
            }
        }
        // printf("number=%d\n", number);
        ped1->SetDensity(number / area);
        // printf("density=%f\n\n", ped1->GetDensity());
    }
    
    // Step 2: calculate the direction of movement for each agent
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

    // Step 3: update the direction of movement for each agent
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped = allPeds[p];
        Point direction = resultDir[p];
        ped->SetMoveDirection(direction);
    }

    // Step 4: calculate the speed for each agent
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
        speed = ei * PushSpeed(ped1, spacing, neighbours,room1);
        resultAcc.push_back(speed);
    }

    // Step 5: calculate the pushing efforts from neighbors and walls
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
            // the wall also has contact force, without looks better]
            //contactForce += ForceConRoom(ped1, subroom1);
        } //for i
         resultForce.push_back(contactForce);
    }


    //Step 6: update everything, position, speed
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
        
        // Avoid agents cross the wall
        Room* room = building->GetRoom(ped->GetRoomID());
        double maxSpace = FLT_MAX;
        for (auto it = room->GetAllSubRooms().begin(); it != room->GetAllSubRooms().end(); ++it)
        {
            SubRoom* subroom = room->GetSubRoom(it->first);
            // here to adjust
            double space2wall = GetSpacingRoomDirection(ped, subroom, dirNeu,false);
            maxSpace = space2wall < maxSpace ? space2wall : maxSpace;
        }
        if (vNeu.Norm() * deltaT > maxSpace)
        {
            vNeu = dirNeu * (maxSpace / deltaT);
        }
        // printf("ped pos=(%f,%f)\n", ped->GetPos()._x, ped->GetPos()._y);
        // printf("direction=(%f,%f)\n", dirNeu._x, dirNeu._y);
        // printf("maxSpace=%f\n\n", maxSpace);
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
    double aForce = GetaPush()/10;
    double DForce = GetDPush()*10;
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

double AVMModel::GetSpacingRoomDirection(Pedestrian* ped, SubRoom* subroom, Point direction, bool checkDoor) const
{
    double spacing = FLT_MAX;
    //first the walls
    for (const auto& wall : subroom->GetAllWalls())
    {
        double distance = GetSpacingWallDirection(ped, wall,direction);
        spacing = spacing > distance ? distance : spacing;
    }
    //then the obstacles
    for (const auto& obst : subroom->GetAllObstacles())
    {
        for (const auto& wall : obst->GetAllWalls())
        {
            double distance = GetSpacingWallDirection(ped, wall, direction);
            spacing = spacing > distance ? distance : spacing;
        }
    }
    //and finally the doors
    if (checkDoor == true)
    {
        for (const auto& goal : subroom->GetAllTransitions())
        {
            double distance = GetSpacingWallDirection(ped, *(static_cast<Line*>(goal)), direction);
            spacing = spacing > distance ? distance : spacing;
        }
    }
    return spacing;
}

double AVMModel::GetSpacingWallDirection(Pedestrian* ped, const Line& l, Point direction) const
{
    double spacing = FLT_MAX;
    Point pp = ped->GetPos();
    Point pt = l.ShortestPoint(ped->GetPos());
    Point p1 = l.GetPoint1();
    Point p2 = l.GetPoint2();
    Point dist = pt - pp;
    Point ei_vertical;
    Point ei = direction;
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
    if (A1_result1 * A1_result2 > 0 && A2_result1 * A2_result2 > 0 && A1_result1 * A2_result1 > 0 && A1_result2 * A2_result2 > 0)
    {
        return spacing;
    }
    double effdis = dist.Norm();
    double cosangle = dist.ScalarProduct(ei) / (dist.Norm() * ei.Norm());
    if (cosangle < 0.00001)
    {
        return spacing;
    }
    spacing = effdis / cosangle-b;
    return spacing;
}


double AVMModel::PushSpeed(Pedestrian* ped, double spacing, vector<Pedestrian*> neighbours, Room* room) const
{
    double meanPlevel = ped->GetP0();
    double plevel = ped->GetPlevel();
    // The plevel of pedestrians is predicted by the RF classifier
    // The RF classifier code is generated by m2cgen from python
    // Calculate the features for prediction-----------------------------------------------------------------------
    vector<double> features = vector<double >();
    int N =2; // the value of N will be determined by the result of mechine learning part
    int featureNum = 1 + N * 4;
    features.reserve(2*featureNum);
    NeighborInfoMLFeatures(features, N, ped, neighbours, room);
    /*
    printf("\nsize is %d\n", features.size());
    for (int i = 0; i < features.size(); i++)
    {
        printf("index %d is %f\n", i, features[i]);
    }
    */
    // Predict the pushing level based on the extracted features------------------------------------------------
    double* input = new double[features.size()];
    memcpy(input, &features[0], features.size() * sizeof(double));
    /* check if input is correct
    for (int i = 0; i < features.size(); i++)
    {
        printf("i=%d,feature=%f,input=%f.\n", i, features[i], input[i]);
    }
    */
    double output=-1;
    double *po=&output;
    score(input, po);
    //the prediction result is pushing: score<0.5, nonpushing: score>0.5
    if (output < 0.5)
    {
        plevel = 1;
        ped->SetPlevel(plevel);
    }
    //----------------------------------------------------------------------------------------------------------------------
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
    std::normal_distribution<double> spaceDist = std::normal_distribution<double>(0,0.1);
    std::normal_distribution<double> TDist = std::normal_distribution<double>(0, 0.1);
    extraSpace = extraSpace + spaceDist(gen);
    extraSpace = extraSpace >2* ped->GetEllipse().GetBmax() ? 2 * ped->GetEllipse().GetBmax() : extraSpace;
    T = T + TDist(gen);
    T = T < 0.01 ? 0.01 : T;
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

void AVMModel::NeighborInfoMLFeatures(vector<double>& features, const int N,  Pedestrian* ped, const vector<Pedestrian*> neighbours, Room* room) const
{
    // the first feature is the mean plevel
    double mplevel = ped->GetP0();
    features.push_back(mplevel);
    // ped information
    Point pos1 = ped->GetPos();
    Point desiredDirection = DesireDirection(ped, room);
    SubRoom* subroom = room->GetSubRoom(ped->GetSubRoomID());
    // space, speed, density, plevel 
    vector<double> feaSpaces = vector<double>();
    vector<double> feaSpeeds = vector<double>();
    vector<double> feaDensitys = vector<double>();
    vector<double> feaPlevels = vector<double>();
    // neighbor information
    int size = (int)neighbours.size();
    for (int index = 0; index < N; index++)
    {
        vector<double> spaces = vector<double>();
        vector<double> speeds = vector<double>();
        vector<double> densitys = vector<double>();
        vector<double> plevels = vector<double>();
        double angle = index * M_PI *2 / N;
        //printf("index=%d, angle=%f\n", index, angle);
        Point direction = desiredDirection.Rotate(cos(angle), sin(angle));
        /*
        printf("angle=%f, cos=%f, sin=%f\n", angle, cos(angle), sin(angle));
        printf("desiredDirection=(%f,%f)\n", desiredDirection._x, desiredDirection._y);
        printf("direction=(%f,%f)\n\n", direction._x, direction._y);
        */
        for (int i = 0; i < size; i++)
        {
            Pedestrian* ped2 = neighbours[i];
            Point pos2 = ped2->GetPos();
            double dist12 = (pos2 - pos1).Norm();
            Point ep12 = (pos2 - pos1).Normalized();
            double temp = ep12.ScalarProduct(direction);
            temp = temp > 1 ? 1 : temp;
            temp = temp < -1 ? -1 : temp;
            double angle2dir = acos(temp);
            // ped2 is in the range
            // printf("ped1=(%f,%f), ped2=(%f,%f)\n", pos1._x, pos1._y, pos2._x, pos2._y);
            // printf("direction=(%f,%f)\n ", direction._x, direction._y);
            if (angle2dir < M_PI / N)
            {
                //printf("ped %d in range\n\n", ped2->GetID());
                spaces.push_back(dist12);
                double projSpeed = ped2->GetV().ScalarProduct(desiredDirection);
                //printf("ped2v=(%f,%f)\n", ped2->GetV()._x, ped2->GetV()._y);
                //printf("desiredDirection=(%f,%f)\n", desiredDirection._x, desiredDirection._y);
                //printf("projSpeed = % f\n\n", projSpeed);
                speeds.push_back(projSpeed);
                double density = ped2->GetDensity();
                densitys.push_back(density);
                // Since we don't have the real plevel, this step unifies the dimensions of simulation and experimental data.
                double plevel = ped2->GetPlevel() + 2;
                plevels.push_back(plevel);
            }
        }
        // the space 2 walls
        for (auto it = room->GetAllSubRooms().begin(); it != room->GetAllSubRooms().end(); ++it)
        {
            SubRoom* subroom2 = room->GetSubRoom(it->first);
            double space2wall = GetSpacingRoomDirection(ped, subroom2, direction,true) + ped->GetEllipse().GetBmax();
            spaces.push_back(space2wall);
        }
        double spaceMin = *min_element(spaces.begin(), spaces.end());
        double speedMean = 0;
        double densityMean = 0;
        double plevelMean = 0;
        if (speeds.size() != 0)
        {
            speedMean = MeanVector(speeds);
            densityMean = MeanVector(densitys);
            plevelMean = MeanVector(plevels);
        }
        feaSpaces.push_back(spaceMin);
        feaSpeeds.push_back(speedMean);
        feaDensitys.push_back(densityMean);
        feaPlevels.push_back(plevelMean);
    }
    for (int i = 0; i < feaSpaces.size(); i++)
    {
        features.push_back(feaSpaces[i]);
    }
    for (int i = 0; i < feaSpeeds.size(); i++)
    {
        features.push_back(feaSpeeds[i]);
    }
    for (int i = 0; i < feaDensitys.size(); i++)
    {
        features.push_back(feaDensitys[i]);
    }
    for (int i = 0; i < feaPlevels.size(); i++)
    {
        features.push_back(feaPlevels[i]);
    }
}

void AVMModel::score(double* input, double* output) const{
    double var0[2];
    double var1[2];
    double var2[2];
    double var3[2];
    double var4[2];
    double var5[2];
    double var6[2];
    double var7[2];
    double var8[2];
    if (input[0] <= 2.4404999017715454) {
        if (input[2] <= 0.06149999983608723) {
            if (input[0] <= 2.149500012397766) {
                if (input[0] <= 2.0754998922348022) {
                    if (input[7] <= 1.8394999504089355) {
                        if (input[1] <= -0.009999999776482582) {
                            if (input[8] <= 1.9555000066757202) {
                                if (input[7] <= 1.7085000276565552) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.0199999809265137) {
                                if (input[1] <= -0.003999999957159162) {
                                    double tempArray[2] = { 0.8666666666666667, 0.13333333333333333 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.30899999290704727) {
                                    double tempArray[2] = { 0.96, 0.04 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 3.6584999561309814) {
                            if (input[0] <= 2.0554999113082886) {
                                if (input[7] <= 2.5700000524520874) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9972789115646259, 0.0027210884353741495 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.24199999868869781) {
                                    double tempArray[2] = { 0.9794275491949911, 0.020572450805008944 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.78125, 0.21875 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.114500045776367) {
                                if (input[1] <= -0.0005000000237487257) {
                                    double tempArray[2] = { 0.7619047619047619, 0.23809523809523808 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9856184084372004, 0.014381591562799617 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.059999942779541) {
                                    double tempArray[2] = { 0.9926569356138623, 0.0073430643861377275 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9080459770114943, 0.09195402298850575 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[8] <= 2.513000011444092) {
                        if (input[8] <= 2.1705000400543213) {
                            if (input[1] <= -0.09149999916553497) {
                                if (input[0] <= 2.1200000047683716) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 6.985499858856201) {
                                    double tempArray[2] = { 0.9765060240963855, 0.023493975903614458 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 5.4084999561309814) {
                                if (input[7] <= 2.165000081062317) {
                                    double tempArray[2] = { 0.14130434782608695, 0.8586956521739131 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6305220883534136, 0.36947791164658633 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.2820000648498535) {
                                    double tempArray[2] = { 0.4585987261146497, 0.5414012738853503 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9553072625698324, 0.0446927374301676 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 5.576499938964844) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[7] <= 2.36299991607666) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.237499952316284) {
                    if (input[8] <= 1.8535000085830688) {
                        if (input[6] <= 3.371999979019165) {
                            if (input[7] <= 2.6705000400543213) {
                                if (input[2] <= 0.05949999950826168) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.975, 0.025 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.050999999046325684) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 3.912000060081482) {
                                if (input[7] <= 2.3179999589920044) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[8] <= 2.2109999656677246) {
                            if (input[5] <= 5.6875) {
                                if (input[6] <= 3.097499966621399) {
                                    double tempArray[2] = { 0.7239077669902912, 0.27609223300970875 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.32697947214076245, 0.6730205278592375 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.016000000294297934) {
                                    double tempArray[2] = { 0.6642857142857143, 0.3357142857142857 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9511363636363637, 0.048863636363636366 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 3.1679999828338623) {
                                if (input[0] <= 2.18149995803833) {
                                    double tempArray[2] = { 0.3504531722054381, 0.649546827794562 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9878048780487805, 0.012195121951219513 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 8.861000061035156) {
                                    double tempArray[2] = { 0.9160903642231443, 0.08390963577685569 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2553191489361702, 0.7446808510638298 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= -0.04350000061094761) {
                        if (input[5] <= 7.0945000648498535) {
                            if (input[0] <= 2.4200000762939453) {
                                if (input[7] <= 2.4854999780654907) {
                                    double tempArray[2] = { 0.5855397148676171, 0.41446028513238287 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3180379746835443, 0.6819620253164557 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.322499990463257) {
                                    double tempArray[2] = { 0.3382352941176471, 0.6617647058823529 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9675324675324676, 0.032467532467532464 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 8.314499855041504) {
                                if (input[0] <= 2.3200000524520874) {
                                    double tempArray[2] = { 0.8962264150943396, 0.10377358490566038 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7123050259965338, 0.2876949740034662 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.5450000762939453) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2755102040816326, 0.7244897959183674 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 2.340499997138977) {
                            if (input[5] <= 4.419500112533569) {
                                if (input[2] <= -0.05849999934434891) {
                                    double tempArray[2] = { 0.9368421052631579, 0.06315789473684211 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.5354999899864197) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6934984520123839, 0.3065015479876161 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.253499984741211) {
                                if (input[8] <= 1.9570000171661377) {
                                    double tempArray[2] = { 0.27769110764430577, 0.7223088923556942 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6244239631336406, 0.37557603686635943 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.385499954223633) {
                                    double tempArray[2] = { 0.7196652719665272, 0.2803347280334728 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.41503604531410915, 0.5849639546858908 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[0] <= 2.0649999380111694) {
                if (input[3] <= 0.03849999979138374) {
                    if (input[8] <= 1.8855000138282776) {
                        if (input[4] <= 0.08150000125169754) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[4] <= 0.08249999955296516) {
                                if (input[7] <= 1.6875) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[5] <= 3.1584999561309814) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[8] <= 1.9259999990463257) {
                                if (input[3] <= -0.04649999924004078) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.952755905511811, 0.047244094488188976 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 1.8604999780654907) {
                                    double tempArray[2] = { 0.9468223086900129, 0.05317769130998703 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9859518348623854, 0.01404816513761468 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[4] <= -0.022499999962747097) {
                        if (input[5] <= 4.183500051498413) {
                            if (input[4] <= -0.04350000061094761) {
                                if (input[0] <= 2.0175000429153442) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3333333333333333, 0.6666666666666666 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.8569999933242798) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7894736842105263, 0.21052631578947367 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[0] <= 2.0175000429153442) {
                            if (input[6] <= 4.200000047683716) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[0] <= 1.9950000047683716) {
                                    double tempArray[2] = { 0.8434782608695652, 0.1565217391304348 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 2.864000082015991) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[5] <= 5.518499851226807) {
                                    double tempArray[2] = { 0.9848866498740554, 0.015113350125944584 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9065420560747663, 0.09345794392523364 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[8] <= 1.8139999508857727) {
                    if (input[7] <= 2.009500026702881) {
                        if (input[0] <= 2.1589999198913574) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[0] <= 2.1804999113082886) {
                            if (input[5] <= 2.790500044822693) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[2] <= 0.2264999970793724) {
                                    double tempArray[2] = { 0.9820971867007673, 0.017902813299232736 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8062283737024222, 0.19377162629757785 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.19849999994039536) {
                                if (input[6] <= 2.705999970436096) {
                                    double tempArray[2] = { 0.9940369707811568, 0.005963029218843173 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8656716417910447, 0.13432835820895522 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.21049999445676804) {
                                    double tempArray[2] = { 0.125, 0.875 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9333333333333333, 0.06666666666666667 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[3] <= 0.13650000095367432) {
                        if (input[0] <= 2.3394999504089355) {
                            if (input[8] <= 2.4570000171661377) {
                                if (input[5] <= 4.534499883651733) {
                                    double tempArray[2] = { 0.8512521182451516, 0.14874788175484843 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6940966010733453, 0.30590339892665475 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.07099999859929085) {
                                    double tempArray[2] = { 0.7948717948717948, 0.20512820512820512 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1357142857142857, 0.8642857142857143 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 3.3765000104904175) {
                                if (input[2] <= 0.30550000071525574) {
                                    double tempArray[2] = { 0.6801573917931422, 0.3198426082068578 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4562289562289562, 0.5437710437710438 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.039499999955296516) {
                                    double tempArray[2] = { 0.0842911877394636, 0.9157088122605364 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.39076923076923076, 0.6092307692307692 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[8] <= 2.031499981880188) {
                            if (input[8] <= 1.918500006198883) {
                                if (input[0] <= 2.1545000076293945) {
                                    double tempArray[2] = { 0.2, 0.8 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8930817610062893, 0.1069182389937107 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.11899995803833) {
                                    double tempArray[2] = { 0.8358445678033307, 0.1641554321966693 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9340241796200346, 0.06597582037996545 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.21249999850988388) {
                                if (input[8] <= 2.0674999952316284) {
                                    double tempArray[2] = { 0.15217391304347827, 0.8478260869565217 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.727896066402021, 0.2721039335979791 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.2459999993443489) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9306930693069307, 0.06930693069306931 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        if (input[0] <= 2.7415000200271606) {
            if (input[6] <= 3.4744999408721924) {
                if (input[5] <= 3.249500036239624) {
                    if (input[3] <= 0.4064999967813492) {
                        if (input[4] <= 0.013500000350177288) {
                            if (input[7] <= 2.1790000200271606) {
                                if (input[6] <= 2.896499991416931) {
                                    double tempArray[2] = { 0.4828767123287671, 0.5171232876712328 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9396984924623115, 0.06030150753768844 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.1799999475479126) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 2.3264999389648438) {
                                if (input[0] <= 2.527999997138977) {
                                    double tempArray[2] = { 0.43478260869565216, 0.5652173913043478 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8536585365853658, 0.14634146341463414 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.3824999928474426) {
                                    double tempArray[2] = { 0.30630372492836677, 0.6936962750716332 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9, 0.1 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= 0.02950000064447522) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[4] <= 0.23600000143051147) {
                                if (input[0] <= 2.4915000200271606) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.3795000314712524) {
                                    double tempArray[2] = { 0.8851674641148325, 0.11483253588516747 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.611999988555908) {
                        if (input[1] <= 0.02049999963492155) {
                            if (input[7] <= 2.0190000534057617) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[7] <= 2.634999990463257) {
                                    double tempArray[2] = { 0.6488824801730353, 0.35111751982696465 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.25732899022801303, 0.742671009771987 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 3.371500015258789) {
                                if (input[3] <= 0.057500001043081284) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5561497326203209, 0.44385026737967914 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.0799999237060547) {
                                    double tempArray[2] = { 0.910086004691165, 0.08991399530883502 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7341153470185728, 0.26588465298142716 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 1.5999999642372131) {
                            if (input[3] <= 0.08449999988079071) {
                                if (input[0] <= 2.6750000715255737) {
                                    double tempArray[2] = { 0.5853658536585366, 0.4146341463414634 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9883720930232558, 0.011627906976744186 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.5080000162124634) {
                                    double tempArray[2] = { 0.5333333333333333, 0.4666666666666667 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.030303030303030304, 0.9696969696969697 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.7345000505447388) {
                                if (input[0] <= 2.677000045776367) {
                                    double tempArray[2] = { 0.5102707749766573, 0.48972922502334265 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1023936170212766, 0.8976063829787234 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.1184999980032444) {
                                    double tempArray[2] = { 0.5208333333333334, 0.4791666666666667 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9166666666666666, 0.08333333333333333 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[3] <= 0.43050000071525574) {
                    if (input[0] <= 2.6880000829696655) {
                        if (input[5] <= 5.870499849319458) {
                            if (input[4] <= 0.2864999920129776) {
                                if (input[1] <= 0.07850000262260437) {
                                    double tempArray[2] = { 0.22134467395706764, 0.7786553260429324 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.46194926568758343, 0.5380507343124166 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.3429999351501465) {
                                    double tempArray[2] = { 0.6571428571428571, 0.34285714285714286 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.968421052631579, 0.031578947368421054 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.009499999694526196) {
                                if (input[7] <= 2.409500002861023) {
                                    double tempArray[2] = { 0.5345765345765345, 0.4654234654234654 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3598404255319149, 0.6401595744680851 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 4.555500030517578) {
                                    double tempArray[2] = { 0.7756613756613756, 0.22433862433862434 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4148296593186373, 0.5851703406813628 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.24550000578165054) {
                            if (input[5] <= 6.661999940872192) {
                                if (input[2] <= -0.16499999910593033) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0264797507788162, 0.9735202492211839 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= -0.018499999307096004) {
                                    double tempArray[2] = { 0.49047619047619045, 0.5095238095238095 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.18175853018372704, 0.818241469816273 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= 0.12449999805539846) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[8] <= 2.4375) {
                        if (input[1] <= -0.02600000100210309) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        double tempArray[2] = { 1.0, 0.0 };
                        memcpy(var8, tempArray, 2 * sizeof(double));
                    }
                }
            }
        }
        else {
            if (input[1] <= 0.13749999552965164) {
                if (input[6] <= 2.5434999465942383) {
                    if (input[8] <= 2.0700000524520874) {
                        if (input[5] <= 7.5325000286102295) {
                            if (input[3] <= 0.20000000298023224) {
                                if (input[7] <= 2.48199999332428) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.042328042328042326, 0.9576719576719577 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.971500039100647) {
                                    double tempArray[2] = { 0.44525547445255476, 0.5547445255474452 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 0.9679999947547913) {
                                if (input[7] <= 2.6855000257492065) {
                                    double tempArray[2] = { 0.8780487804878049, 0.12195121951219512 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.14634146341463414, 0.8536585365853658 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 2.343500018119812) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0975609756097561, 0.9024390243902439 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.7400000095367432) {
                            if (input[8] <= 2.174499988555908) {
                                if (input[0] <= 2.91949999332428) {
                                    double tempArray[2] = { 0.9, 0.1 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3300970873786408, 0.6699029126213593 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.0025000000605359674) {
                                    double tempArray[2] = { 0.19763513513513514, 0.8023648648648649 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.36670416197975253, 0.6332958380202475 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.3765000104904175) {
                                if (input[2] <= 0.21849999576807022) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[8] <= 2.989500045776367) {
                        if (input[0] <= 2.853999972343445) {
                            if (input[4] <= 0.26900000870227814) {
                                if (input[7] <= 2.6894999742507935) {
                                    double tempArray[2] = { 0.18532374100719423, 0.8146762589928057 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.06619504800404244, 0.9338049519959576 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.06550000235438347) {
                                    double tempArray[2] = { 0.5281173594132029, 0.4718826405867971 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.03614457831325301, 0.963855421686747 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= 0.08749999850988388) {
                                if (input[4] <= 0.28199999034404755) {
                                    double tempArray[2] = { 0.0010878104792409501, 0.998912189520759 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.04868913857677903, 0.951310861423221 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.33949999511241913) {
                                    double tempArray[2] = { 0.013367080597959345, 0.9866329194020407 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.07160555004955402, 0.928394449950446 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.805500030517578) {
                            if (input[2] <= 0.18450000137090683) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[1] <= -0.08650000020861626) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3225806451612903, 0.6774193548387096 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= -0.15299999713897705) {
                                if (input[5] <= 6.419000148773193) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 6.42300009727478) {
                                    double tempArray[2] = { 0.027868852459016394, 0.9721311475409836 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[7] <= 2.1454999446868896) {
                    if (input[1] <= 0.14750000089406967) {
                        if (input[8] <= 2.166499972343445) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[4] <= 0.1055000014603138) {
                                if (input[8] <= 2.2135000228881836) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[3] <= 0.018499999307096004) {
                            if (input[0] <= 2.89900004863739) {
                                if (input[6] <= 2.354499936103821) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6551724137931034, 0.3448275862068966 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var8, tempArray, 2 * sizeof(double));
                        }
                    }
                }
                else {
                    if (input[0] <= 2.91949999332428) {
                        if (input[6] <= 2.9364999532699585) {
                            if (input[7] <= 2.6994999647140503) {
                                if (input[0] <= 2.7610000371932983) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.962877030162413, 0.037122969837587005 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[2] <= -0.06749999895691872) {
                                if (input[0] <= 2.887500047683716) {
                                    double tempArray[2] = { 0.9444444444444444, 0.05555555555555555 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.1600000038743019) {
                                    double tempArray[2] = { 0.2, 0.8 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.047244094488188976, 0.952755905511811 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.627500057220459) {
                            if (input[8] <= 2.1309999227523804) {
                                if (input[0] <= 2.9674999713897705) {
                                    double tempArray[2] = { 0.6808510638297872, 0.3191489361702128 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.027093596059113302, 0.9729064039408867 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 2.2489999532699585) {
                                    double tempArray[2] = { 0.8031496062992126, 0.1968503937007874 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.25180897250361794, 0.748191027496382 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= -0.20899999886751175) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var8, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[3] <= 0.28550000488758087) {
                                    double tempArray[2] = { 0.04487179487179487, 0.9551282051282052 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var8, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    double var9[2];
    if (input[0] <= 2.438499927520752) {
        if (input[1] <= 0.021499999798834324) {
            if (input[2] <= 0.11649999767541885) {
                if (input[0] <= 2.174499988555908) {
                    if (input[7] <= 2.009500026702881) {
                        if (input[6] <= 2.6679999828338623) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[0] <= 2.149500012397766) {
                                if (input[1] <= -0.018499999307096004) {
                                    double tempArray[2] = { 0.7297921478060047, 0.2702078521939954 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9301397205588823, 0.06986027944111776 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.09950000047683716) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.17117117117117117, 0.8288288288288288 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= -0.05849999934434891) {
                            if (input[0] <= 2.1675000190734863) {
                                if (input[0] <= 2.072999954223633) {
                                    double tempArray[2] = { 0.997332147621165, 0.0026678523788350376 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9799833194328608, 0.020016680567139282 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.802500009536743) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.0754998922348022) {
                                if (input[3] <= 0.05049999989569187) {
                                    double tempArray[2] = { 0.9841269841269841, 0.015873015873015872 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9950031230480949, 0.004996876951905059 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.8535000085830688) {
                                    double tempArray[2] = { 0.9986111111111111, 0.001388888888888889 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7545038167938931, 0.24549618320610686 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[5] <= 5.327500104904175) {
                        if (input[7] <= 2.8214999437332153) {
                            if (input[6] <= 3.0449999570846558) {
                                if (input[7] <= 2.5049999952316284) {
                                    double tempArray[2] = { 0.713855421686747, 0.286144578313253 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5143853530950305, 0.48561464690496947 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.0390000343322754) {
                                    double tempArray[2] = { 0.09420289855072464, 0.9057971014492754 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3938931297709924, 0.6061068702290077 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 1.878499984741211) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 4.165999889373779) {
                                    double tempArray[2] = { 0.038461538461538464, 0.9615384615384616 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3384615384615385, 0.6615384615384615 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= -0.05650000087916851) {
                            if (input[5] <= 7.5909998416900635) {
                                if (input[6] <= 2.430500030517578) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.16853932584269662, 0.8314606741573034 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.01000000024214387) {
                                    double tempArray[2] = { 0.8867924528301887, 0.11320754716981132 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.26785714285714285, 0.7321428571428571 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 5.625499963760376) {
                                if (input[6] <= 2.715499997138977) {
                                    double tempArray[2] = { 0.825211176088369, 0.17478882391163092 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6245358681122828, 0.3754641318877172 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.6829999685287476) {
                                    double tempArray[2] = { 0.7832310838445807, 0.2167689161554192 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.991304347826087, 0.008695652173913044 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[8] <= 1.8149999976158142) {
                    if (input[0] <= 2.0679999589920044) {
                        double tempArray[2] = { 1.0, 0.0 };
                        memcpy(var9, tempArray, 2 * sizeof(double));
                    }
                    else {
                        if (input[7] <= 2.009500026702881) {
                            if (input[0] <= 2.1589999198913574) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[3] <= 0.3409999907016754) {
                                if (input[8] <= 1.7319999933242798) {
                                    double tempArray[2] = { 0.9979508196721312, 0.0020491803278688526 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8703703703703703, 0.12962962962962962 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[7] <= 2.367500066757202) {
                        if (input[1] <= -0.07349999994039536) {
                            if (input[2] <= 0.12150000035762787) {
                                if (input[0] <= 2.0785000324249268) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 1.9164999723434448) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9273743016759777, 0.07262569832402235 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.309499979019165) {
                                if (input[7] <= 1.7319999933242798) {
                                    double tempArray[2] = { 0.1, 0.9 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9553025511279781, 0.04469744887202193 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.0480000972747803) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6041666666666666, 0.3958333333333333 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[3] <= 0.13250000029802322) {
                            if (input[8] <= 2.341499924659729) {
                                if (input[5] <= 3.6735000610351562) {
                                    double tempArray[2] = { 0.5447154471544715, 0.45528455284552843 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7747105966162066, 0.2252894033837934 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.02649998664856) {
                                    double tempArray[2] = { 0.6164383561643836, 0.3835616438356164 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9809725158562368, 0.019027484143763214 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 5.061500072479248) {
                                if (input[8] <= 2.159500002861023) {
                                    double tempArray[2] = { 0.9161290322580645, 0.08387096774193549 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6739130434782609, 0.32608695652173914 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.8415000438690186) {
                                    double tempArray[2] = { 0.9926793557833089, 0.007320644216691069 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8783783783783784, 0.12162162162162163 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[0] <= 2.070499897003174) {
                if (input[0] <= 2.0175000429153442) {
                    if (input[2] <= -0.09050000086426735) {
                        if (input[8] <= 2.395500063896179) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[7] <= 2.419000029563904) {
                                if (input[1] <= 0.08949999883770943) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.774999976158142) {
                                    double tempArray[2] = { 0.7972972972972973, 0.20270270270270271 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[1] <= 0.02650000061839819) {
                            if (input[3] <= 0.019499999471008778) {
                                if (input[3] <= 0.018499999307096004) {
                                    double tempArray[2] = { 0.9885714285714285, 0.011428571428571429 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8571428571428571, 0.14285714285714285 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.022499999962747097) {
                                    double tempArray[2] = { 0.9962825278810409, 0.0037174721189591076 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 3.869499921798706) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[0] <= 1.9950000047683716) {
                                    double tempArray[2] = { 0.9949044585987261, 0.005095541401273885 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[7] <= 2.0334999561309814) {
                        if (input[2] <= 0.13450000435113907) {
                            if (input[1] <= 0.13899999856948853) {
                                if (input[4] <= 0.051500000059604645) {
                                    double tempArray[2] = { 0.8166666666666667, 0.18333333333333332 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[1] <= 0.17399999499320984) {
                                if (input[0] <= 2.0644999742507935) {
                                    double tempArray[2] = { 0.11627906976744186, 0.8837209302325582 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.049999997951090336) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.0625) {
                            if (input[6] <= 2.04449999332428) {
                                if (input[5] <= 5.382499933242798) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8545454545454545, 0.14545454545454545 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= -0.0625) {
                                    double tempArray[2] = { 0.8333333333333334, 0.16666666666666666 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.8855000734329224) {
                                if (input[6] <= 2.212000012397766) {
                                    double tempArray[2] = { 0.6785714285714286, 0.32142857142857145 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9764837625979843, 0.023516237402015677 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.04350000061094761) {
                                    double tempArray[2] = { 0.2727272727272727, 0.7272727272727273 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.3289999961853027) {
                    if (input[8] <= 2.2020000219345093) {
                        if (input[4] <= 0.6745000183582306) {
                            if (input[6] <= 5.2225000858306885) {
                                if (input[5] <= 7.364000082015991) {
                                    double tempArray[2] = { 0.8710901152387087, 0.12890988476129125 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[7] <= 2.4485000371932983) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[8] <= 2.1714999675750732) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 5.444000005722046) {
                            if (input[2] <= -0.08749999850988388) {
                                if (input[1] <= 0.24049999564886093) {
                                    double tempArray[2] = { 0.9585798816568047, 0.04142011834319527 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.8609999418258667) {
                                    double tempArray[2] = { 0.6897935779816514, 0.31020642201834864 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2028985507246377, 0.7971014492753623 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= -0.06750000268220901) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 3.4884999990463257) {
                                    double tempArray[2] = { 0.7474402730375427, 0.2525597269624573 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8970212765957447, 0.10297872340425532 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[4] <= 0.34300000965595245) {
                        if (input[3] <= 0.1614999994635582) {
                            if (input[6] <= 4.326500177383423) {
                                if (input[2] <= 0.3244999945163727) {
                                    double tempArray[2] = { 0.6892464013547841, 0.3107535986452159 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.38687782805429866, 0.6131221719457014 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.4200000762939453) {
                                    double tempArray[2] = { 0.26666666666666666, 0.7333333333333333 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7317073170731707, 0.2682926829268293 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= -0.02949999924749136) {
                                if (input[8] <= 2.1304999589920044) {
                                    double tempArray[2] = { 0.6981132075471698, 0.3018867924528302 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1919191919191919, 0.8080808080808081 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 6.682500123977661) {
                                    double tempArray[2] = { 0.9108433734939759, 0.0891566265060241 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.0740000009536743) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[7] <= 2.5570000410079956) {
                                if (input[3] <= -0.0010000001639127731) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9801587301587301, 0.01984126984126984 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.5695000886917114) {
                                    double tempArray[2] = { 0.023809523809523808, 0.9761904761904762 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8733153638814016, 0.12668463611859837 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        if (input[0] <= 2.7615000009536743) {
            if (input[1] <= 0.022499999962747097) {
                if (input[0] <= 2.618499994277954) {
                    if (input[2] <= 0.02650000061839819) {
                        if (input[5] <= 5.672499895095825) {
                            if (input[4] <= 0.18650000542402267) {
                                if (input[1] <= -0.060499999672174454) {
                                    double tempArray[2] = { 0.3205574912891986, 0.6794425087108014 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.09874917709019092, 0.9012508229098091 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.18050000071525574) {
                                    double tempArray[2] = { 0.8297872340425532, 0.1702127659574468 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3333333333333333, 0.6666666666666666 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.472999930381775) {
                                if (input[8] <= 2.072499990463257) {
                                    double tempArray[2] = { 0.7744107744107744, 0.2255892255892256 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.16581632653061223, 0.8341836734693877 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 5.371500015258789) {
                                    double tempArray[2] = { 0.3333333333333333, 0.6666666666666666 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7152145643693107, 0.2847854356306892 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 2.93149995803833) {
                            if (input[4] <= 0.07850000262260437) {
                                if (input[2] <= 0.14299999922513962) {
                                    double tempArray[2] = { 0.10112359550561797, 0.898876404494382 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6082474226804123, 0.3917525773195876 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.34299999475479126) {
                                    double tempArray[2] = { 0.05396825396825397, 0.946031746031746 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6923076923076923, 0.3076923076923077 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.3550000190734863) {
                                if (input[8] <= 2.0549999475479126) {
                                    double tempArray[2] = { 0.6852818371607515, 0.3147181628392484 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5490033222591362, 0.45099667774086377 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.04050000011920929) {
                                    double tempArray[2] = { 0.17959183673469387, 0.8204081632653061 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.45531197301854975, 0.5446880269814502 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 1.6305000185966492) {
                        if (input[7] <= 2.4730000495910645) {
                            if (input[5] <= 7.136499881744385) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[1] <= -0.027500000782310963) {
                            if (input[5] <= 7.542500019073486) {
                                if (input[3] <= 0.018499999307096004) {
                                    double tempArray[2] = { 0.4530456852791878, 0.5469543147208121 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.19563636363636364, 0.8043636363636364 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.628499984741211) {
                                    double tempArray[2] = { 0.13142375737152484, 0.8685762426284751 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.009324009324009324, 0.9906759906759907 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 1.9369999766349792) {
                                if (input[8] <= 1.8855000138282776) {
                                    double tempArray[2] = { 0.509090909090909, 0.4909090909090909 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.06949999928474426) {
                                    double tempArray[2] = { 0.5378548895899053, 0.46214511041009465 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2923076923076923, 0.7076923076923077 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.5460000038146973) {
                    if (input[4] <= 0.32250000536441803) {
                        if (input[6] <= 4.694499969482422) {
                            if (input[2] <= 0.22749999910593033) {
                                if (input[5] <= 5.871500015258789) {
                                    double tempArray[2] = { 0.5966057441253264, 0.40339425587467365 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8865740740740741, 0.11342592592592593 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.424999952316284) {
                                    double tempArray[2] = { 0.36904761904761907, 0.6309523809523809 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9153846153846154, 0.08461538461538462 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.4644999504089355) {
                                if (input[8] <= 2.3639999628067017) {
                                    double tempArray[2] = { 0.5568862275449101, 0.4431137724550898 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.19008264462809918, 0.8099173553719008 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.396499991416931) {
                                    double tempArray[2] = { 0.21052631578947367, 0.7894736842105263 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= 0.0625) {
                            if (input[0] <= 2.4529999494552612) {
                                if (input[6] <= 4.451499938964844) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.3375000059604645) {
                                    double tempArray[2] = { 0.9536679536679536, 0.04633204633204633 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.13043478260869565, 0.8695652173913043 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 6.648999929428101) {
                                if (input[2] <= 0.18249999731779099) {
                                    double tempArray[2] = { 0.9682926829268292, 0.03170731707317073 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7909836065573771, 0.20901639344262296 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 5.354000091552734) {
                        if (input[8] <= 1.9070000052452087) {
                            if (input[7] <= 2.1890000104904175) {
                                if (input[1] <= 0.13300000131130219) {
                                    double tempArray[2] = { 0.6808510638297872, 0.3191489361702128 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.815999984741211) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7234513274336283, 0.27654867256637167 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 2.7975000143051147) {
                                if (input[2] <= 0.45650000870227814) {
                                    double tempArray[2] = { 0.3331608898085877, 0.6668391101914123 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5864022662889519, 0.41359773371104813 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.28849999606609344) {
                                    double tempArray[2] = { 0.5178735105407882, 0.4821264894592117 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9071274298056156, 0.09287257019438445 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.2810000032186508) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[2] <= -0.08999999985098839) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[1] <= 0.13350000232458115) {
                if (input[4] <= 0.2814999967813492) {
                    if (input[2] <= 4.335999965667725) {
                        if (input[3] <= -0.13749999552965164) {
                            if (input[2] <= 0.01250000111758709) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[2] <= 0.17549999803304672) {
                                if (input[6] <= 2.9665000438690186) {
                                    double tempArray[2] = { 0.10301302931596092, 0.8969869706840391 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.01750438436848549, 0.9824956156315146 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.919000029563904) {
                                    double tempArray[2] = { 0.22208281053952322, 0.7779171894604768 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0224642614023145, 0.9775357385976855 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.5334999561309814) {
                            if (input[5] <= 7.335500001907349) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                    }
                }
                else {
                    if (input[0] <= 2.9809999465942383) {
                        if (input[8] <= 2.6179999113082886) {
                            if (input[4] <= 0.6169999837875366) {
                                if (input[0] <= 2.784500002861023) {
                                    double tempArray[2] = { 0.7227722772277227, 0.27722772277227725 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3329694323144105, 0.6670305676855895 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.07349999994039536) {
                                    double tempArray[2] = { 0.9381443298969072, 0.061855670103092786 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5507246376811594, 0.4492753623188406 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.472000002861023) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[0] <= 3.027500033378601) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[6] <= 6.525500059127808) {
                                if (input[0] <= 3.134500026702881) {
                                    double tempArray[2] = { 0.4890282131661442, 0.5109717868338558 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.014768624876928127, 0.9852313751230719 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.02000000048428774) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[7] <= 2.1454999446868896) {
                    if (input[1] <= 0.14549999684095383) {
                        if (input[0] <= 2.927500009536743) {
                            if (input[4] <= 0.058999999426305294) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[6] <= 3.1829999685287476) {
                            if (input[0] <= 2.927500009536743) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                    }
                }
                else {
                    if (input[7] <= 2.606500029563904) {
                        if (input[3] <= 0.19849999994039536) {
                            if (input[1] <= 0.28550000488758087) {
                                if (input[6] <= 1.712499976158142) {
                                    double tempArray[2] = { 0.9285714285714286, 0.07142857142857142 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.11578947368421053, 0.8842105263157894 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.962499976158142) {
                                    double tempArray[2] = { 0.8490566037735849, 0.1509433962264151 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.30666666666666664, 0.6933333333333334 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= 0.39800000190734863) {
                                if (input[7] <= 2.416000008583069) {
                                    double tempArray[2] = { 0.5728155339805825, 0.42718446601941745 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.32679738562091504, 0.673202614379085 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.3525000065565109) {
                                    double tempArray[2] = { 0.795774647887324, 0.20422535211267606 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.25, 0.75 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.45400001108646393) {
                            if (input[3] <= 0.11349999904632568) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var9, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[5] <= 5.608000040054321) {
                                    double tempArray[2] = { 0.2792452830188679, 0.720754716981132 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.02100840336134454, 0.9789915966386554 };
                                    memcpy(var9, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var9, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
        }
    }
    add_vectors(var8, var9, 2, var7);
    double var10[2];
    if (input[0] <= 2.4404999017715454) {
        if (input[7] <= 2.4954999685287476) {
            if (input[1] <= 0.060499999672174454) {
                if (input[8] <= 1.7944999933242798) {
                    if (input[5] <= 7.513499975204468) {
                        if (input[0] <= 2.0460000038146973) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var10, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[7] <= 2.009500026702881) {
                                if (input[3] <= 0.2005000039935112) {
                                    double tempArray[2] = { 0.7933333333333333, 0.20666666666666667 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.06666666666666667, 0.9333333333333333 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 2.6924999952316284) {
                                    double tempArray[2] = { 0.980922650017343, 0.019077349982656956 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8742724097788126, 0.12572759022118743 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        double tempArray[2] = { 0.0, 1.0 };
                        memcpy(var10, tempArray, 2 * sizeof(double));
                    }
                }
                else {
                    if (input[0] <= 2.070499897003174) {
                        if (input[7] <= 1.8449999690055847) {
                            if (input[6] <= 4.1214998960494995) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 4.648999929428101) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 4.188500165939331) {
                                if (input[0] <= 2.0554999113082886) {
                                    double tempArray[2] = { 0.9979237317461416, 0.0020762682538583983 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9394930498773508, 0.06050695012264922 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 4.632499933242798) {
                                    double tempArray[2] = { 0.907608695652174, 0.09239130434782608 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 5.492000102996826) {
                            if (input[5] <= 4.11050009727478) {
                                if (input[2] <= 0.06750000268220901) {
                                    double tempArray[2] = { 0.5804899387576553, 0.4195100612423447 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8112582781456954, 0.18874172185430463 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.2615000009536743) {
                                    double tempArray[2] = { 0.2478017585931255, 0.7521982414068745 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6511175898931001, 0.3488824101068999 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= -0.05550000071525574) {
                                if (input[3] <= -0.057499999180436134) {
                                    double tempArray[2] = { 0.23809523809523808, 0.7619047619047619 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9322493224932249, 0.06775067750677506 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 6.064499855041504) {
                                    double tempArray[2] = { 0.7675263349073738, 0.23247366509262624 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4391891891891892, 0.5608108108108109 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[1] <= 0.10949999839067459) {
                    if (input[6] <= 2.7484999895095825) {
                        if (input[3] <= 0.09449999779462814) {
                            if (input[4] <= 0.4060000032186508) {
                                if (input[0] <= 2.35699999332428) {
                                    double tempArray[2] = { 0.978923311325546, 0.021076688674454037 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.08823529411764706, 0.9117647058823529 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[4] <= 0.38099999725818634) {
                                if (input[3] <= 0.24849999696016312) {
                                    double tempArray[2] = { 0.8900891169322172, 0.10991088306778288 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9629629629629629, 0.037037037037037035 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.13099999725818634) {
                                    double tempArray[2] = { 0.6875, 0.3125 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.2605000734329224) {
                            if (input[0] <= 1.6009999513626099) {
                                if (input[5] <= 5.112499952316284) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5471698113207547, 0.4528301886792453 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.0199999809265137) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9455864570737605, 0.05441354292623942 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 1.881500005722046) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 4.797500133514404) {
                                    double tempArray[2] = { 0.8401122019635343, 0.15988779803646563 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1891891891891892, 0.8108108108108109 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[8] <= 2.2454999685287476) {
                        if (input[2] <= 0.22350000590085983) {
                            if (input[3] <= 0.001500000071246177) {
                                if (input[0] <= 2.0199999809265137) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6766169154228856, 0.32338308457711445 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 1.25) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.990057803468208, 0.009942196531791908 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.121999979019165) {
                                if (input[0] <= 2.0175000429153442) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8837837837837837, 0.11621621621621622 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.009500026702881) {
                                    double tempArray[2] = { 0.38271604938271603, 0.6172839506172839 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7896296296296297, 0.21037037037037037 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 3.969499945640564) {
                            if (input[0] <= 2.0719999074935913) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[0] <= 2.2605000734329224) {
                                    double tempArray[2] = { 0.6666666666666666, 0.3333333333333333 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9531568228105907, 0.04684317718940937 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.2790000438690186) {
                                if (input[6] <= 4.090500116348267) {
                                    double tempArray[2] = { 0.7692307692307693, 0.23076923076923078 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.27149999886751175) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[0] <= 2.1654999256134033) {
                if (input[5] <= 5.807500123977661) {
                    if (input[7] <= 2.9884999990463257) {
                        if (input[8] <= 2.1230000257492065) {
                            if (input[5] <= 5.693499803543091) {
                                if (input[0] <= 2.149999976158142) {
                                    double tempArray[2] = { 0.998910527032548, 0.001089472967451995 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9212598425196851, 0.07874015748031496 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.013500000350177288) {
                                    double tempArray[2] = { 0.9431524547803618, 0.056847545219638244 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.1320000886917114) {
                                if (input[0] <= 2.0649999380111694) {
                                    double tempArray[2] = { 0.9870561483295435, 0.012943851670456533 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8414872798434442, 0.15851272015655576 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.2480000257492065) {
                                    double tempArray[2] = { 0.1488095238095238, 0.8511904761904762 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6615969581749049, 0.33840304182509506 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[8] <= 1.9039999842643738) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var10, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[8] <= 2.1394999027252197) {
                                if (input[8] <= 2.1234999895095825) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.1610000133514404) {
                        if (input[6] <= 7.534499883651733) {
                            if (input[0] <= 2.137500047683716) {
                                if (input[7] <= 2.6160000562667847) {
                                    double tempArray[2] = { 0.9951792708647182, 0.004820729135281711 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.950000047683716) {
                                    double tempArray[2] = { 0.9754273504273504, 0.024572649572649572 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5833333333333334, 0.4166666666666667 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.6980000734329224) {
                                if (input[7] <= 2.694000005722046) {
                                    double tempArray[2] = { 0.8484848484848485, 0.15151515151515152 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.15789473684210525, 0.8421052631578947 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.7445000410079956) {
                                    double tempArray[2] = { 0.9178082191780822, 0.0821917808219178 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= -0.006500000134110451) {
                            if (input[5] <= 6.376500129699707) {
                                if (input[2] <= -0.022499999962747097) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.17647058823529413, 0.8235294117647058 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[1] <= -0.0690000019967556) {
                                if (input[5] <= 6.801500082015991) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8732394366197183, 0.1267605633802817 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 2.9739999771118164) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7916666666666666, 0.20833333333333334 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[1] <= -0.04450000077486038) {
                    if (input[8] <= 2.5740000009536743) {
                        if (input[5] <= 7.398499965667725) {
                            if (input[4] <= 0.2344999983906746) {
                                if (input[8] <= 2.4220000505447388) {
                                    double tempArray[2] = { 0.41877411058722674, 0.5812258894127732 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.15289982425307558, 0.8471001757469244 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.3785001039505005) {
                                    double tempArray[2] = { 0.51875, 0.48125 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9841269841269841, 0.015873015873015872 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= 0.06750000268220901) {
                                if (input[0] <= 2.3200000524520874) {
                                    double tempArray[2] = { 0.8661417322834646, 0.13385826771653545 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5836298932384342, 0.41637010676156583 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.10449999943375587) {
                                    double tempArray[2] = { 0.3533834586466165, 0.6466165413533834 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.846500039100647) {
                            if (input[3] <= -0.05949999950826168) {
                                if (input[6] <= 4.236000061035156) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.7300000190734863) {
                                    double tempArray[2] = { 0.7943262411347518, 0.20567375886524822 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9783783783783784, 0.021621621621621623 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.3279999494552612) {
                                if (input[2] <= -0.04950000159442425) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.1444999948143959) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[3] <= -0.010499999858438969) {
                        if (input[4] <= 0.001500000071246177) {
                            if (input[8] <= 2.640999913215637) {
                                if (input[7] <= 2.575500011444092) {
                                    double tempArray[2] = { 0.009615384615384616, 0.9903846153846154 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.26666666666666666, 0.7333333333333333 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[7] <= 2.9179999828338623) {
                                if (input[2] <= 0.02350000012665987) {
                                    double tempArray[2] = { 0.4735099337748344, 0.5264900662251656 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8969072164948454, 0.10309278350515463 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.13449999690055847) {
                                    double tempArray[2] = { 0.08823529411764706, 0.9117647058823529 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.75, 0.25 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[8] <= 1.8945000171661377) {
                            if (input[8] <= 1.5855000019073486) {
                                if (input[0] <= 2.2475000619888306) {
                                    double tempArray[2] = { 0.989247311827957, 0.010752688172043012 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.19148936170212766, 0.8085106382978723 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.1745000034570694) {
                                    double tempArray[2] = { 0.9374358974358974, 0.06256410256410257 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4, 0.6 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.856000065803528) {
                                if (input[8] <= 2.2300000190734863) {
                                    double tempArray[2] = { 0.5684210526315789, 0.43157894736842106 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.784513933205701, 0.21548606679429907 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.485499858856201) {
                                    double tempArray[2] = { 0.315450643776824, 0.6845493562231759 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.62402496099844, 0.37597503900156004 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        if (input[0] <= 2.8289999961853027) {
            if (input[6] <= 3.4714999198913574) {
                if (input[7] <= 2.634999990463257) {
                    if (input[0] <= 2.503499984741211) {
                        if (input[2] <= 0.2639999985694885) {
                            if (input[8] <= 2.5609999895095825) {
                                if (input[0] <= 2.447999954223633) {
                                    double tempArray[2] = { 0.5354330708661418, 0.4645669291338583 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.895866454689984, 0.1041335453100159 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.600000023841858) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.10869565217391304, 0.8913043478260869 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.4674999713897705) {
                                if (input[7] <= 2.150499939918518) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7065527065527065, 0.2934472934472934 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[5] <= 3.130500078201294) {
                            if (input[1] <= 0.13250000029802322) {
                                if (input[4] <= 0.06549999862909317) {
                                    double tempArray[2] = { 0.4254658385093168, 0.5745341614906833 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1300251256281407, 0.8699748743718593 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 2.5325000286102295) {
                                    double tempArray[2] = { 0.7833089311859444, 0.21669106881405564 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.36671575846833576, 0.6332842415316642 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.011500000022351742) {
                                if (input[5] <= 3.6390000581741333) {
                                    double tempArray[2] = { 0.6697713801862828, 0.3302286198137172 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.397342012476268, 0.602657987523732 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= -0.022499999962747097) {
                                    double tempArray[2] = { 0.010869565217391304, 0.9891304347826086 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.671803127874885, 0.328196872125115 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.7445000410079956) {
                        if (input[2] <= 0.4545000046491623) {
                            if (input[2] <= -0.036000000312924385) {
                                if (input[4] <= 0.27949999272823334) {
                                    double tempArray[2] = { 0.9791666666666666, 0.020833333333333332 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.125, 0.875 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.011500000022351742) {
                                    double tempArray[2] = { 0.16956521739130434, 0.8304347826086956 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5025906735751295, 0.49740932642487046 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.575500011444092) {
                                if (input[0] <= 2.618499994277954) {
                                    double tempArray[2] = { 0.22916666666666666, 0.7708333333333334 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9784366576819407, 0.0215633423180593 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.7044999897480011) {
                                    double tempArray[2] = { 0.030612244897959183, 0.9693877551020408 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[1] <= -0.15399999916553497) {
                            if (input[0] <= 2.8065000772476196) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[2] <= 4.3424999713897705) {
                                if (input[0] <= 2.7730000019073486) {
                                    double tempArray[2] = { 0.09278350515463918, 0.9072164948453608 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.6880000829696655) {
                    if (input[1] <= 0.08349999785423279) {
                        if (input[5] <= 5.8974997997283936) {
                            if (input[0] <= 2.4609999656677246) {
                                if (input[2] <= 0.03749999962747097) {
                                    double tempArray[2] = { 0.2822857142857143, 0.7177142857142857 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7952755905511811, 0.2047244094488189 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.2094999998807907) {
                                    double tempArray[2] = { 0.16786689843555996, 0.83213310156444 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9188311688311688, 0.08116883116883117 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.6160000562667847) {
                                if (input[6] <= 4.594499826431274) {
                                    double tempArray[2] = { 0.627430373095113, 0.37256962690488704 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4057971014492754, 0.5942028985507246 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.9354999661445618) {
                                    double tempArray[2] = { 0.6040955631399317, 0.39590443686006827 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2504109589041096, 0.7495890410958904 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.0845000743865967) {
                            if (input[5] <= 4.598999977111816) {
                                if (input[0] <= 2.462499976158142) {
                                    double tempArray[2] = { 0.8809523809523809, 0.11904761904761904 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.16304347826086957, 0.8369565217391305 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[7] <= 2.7300000190734863) {
                                if (input[6] <= 4.433000087738037) {
                                    double tempArray[2] = { 0.8654618473895582, 0.13453815261044177 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.48066298342541436, 0.5193370165745856 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.4749999940395355) {
                                    double tempArray[2] = { 0.27979274611398963, 0.7202072538860104 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[5] <= 5.478500127792358) {
                        if (input[2] <= -0.07050000131130219) {
                            if (input[5] <= 3.630500078201294) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[5] <= 5.409999847412109) {
                                if (input[3] <= 0.2800000011920929) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.015723270440251572, 0.9842767295597484 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.2044999971985817) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.2775000035762787) {
                            if (input[0] <= 2.759999990463257) {
                                if (input[5] <= 8.153500080108643) {
                                    double tempArray[2] = { 0.27619760479041916, 0.7238023952095808 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.06181015452538632, 0.9381898454746137 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.6894999742507935) {
                                    double tempArray[2] = { 0.20591715976331362, 0.7940828402366864 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.010542168674698794, 0.9894578313253012 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.78249990940094) {
                                if (input[7] <= 2.6435000896453857) {
                                    double tempArray[2] = { 0.45263157894736844, 0.5473684210526316 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9831460674157303, 0.016853932584269662 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[5] <= 5.697499990463257) {
                if (input[6] <= 2.58050000667572) {
                    if (input[7] <= 2.6785000562667847) {
                        if (input[3] <= 0.18949999660253525) {
                            if (input[5] <= 2.3174999952316284) {
                                if (input[7] <= 2.2200000286102295) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6666666666666666, 0.3333333333333333 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.6029999256134033) {
                                    double tempArray[2] = { 0.10946196660482375, 0.8905380333951762 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.64, 0.36 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 1.6349999904632568) {
                                if (input[1] <= 0.22450000047683716) {
                                    double tempArray[2] = { 0.07926829268292683, 0.9207317073170732 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6470588235294118, 0.35294117647058826 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.9809999465942383) {
                                    double tempArray[2] = { 0.8988505747126436, 0.10114942528735632 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.26360338573155984, 0.7363966142684402 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        double tempArray[2] = { 0.0, 1.0 };
                        memcpy(var10, tempArray, 2 * sizeof(double));
                    }
                }
                else {
                    if (input[0] <= 2.9524999856948853) {
                        if (input[7] <= 2.6729999780654907) {
                            if (input[4] <= 0.29999999701976776) {
                                if (input[3] <= -0.06949999928474426) {
                                    double tempArray[2] = { 0.8666666666666667, 0.13333333333333333 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.14381270903010032, 0.8561872909698997 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.5759999752044678) {
                                    double tempArray[2] = { 0.14285714285714285, 0.8571428571428571 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8478260869565217, 0.15217391304347827 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.3899999856948853) {
                                if (input[8] <= 2.3279999494552612) {
                                    double tempArray[2] = { 0.0196078431372549, 0.9803921568627451 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3877551020408163, 0.6122448979591837 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.6985000371932983) {
                                    double tempArray[2] = { 0.16666666666666666, 0.8333333333333334 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 3.124500036239624) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var10, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[0] <= 3.347499966621399) {
                                if (input[8] <= 2.361999988555908) {
                                    double tempArray[2] = { 0.015151515151515152, 0.9848484848484849 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.27522935779816515, 0.7247706422018348 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.853999972343445) {
                    if (input[3] <= -0.031499999575316906) {
                        if (input[8] <= 2.4524999856948853) {
                            if (input[1] <= -0.028999999165534973) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var10, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[2] <= -0.07750000059604645) {
                            if (input[8] <= 2.3589999675750732) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var10, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 5.765499830245972) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.23529411764705882, 0.7647058823529411 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= 0.10649999976158142) {
                                if (input[2] <= -0.04649999924004078) {
                                    double tempArray[2] = { 0.09696969696969697, 0.9030303030303031 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 5.603000164031982) {
                                    double tempArray[2] = { 0.11228070175438597, 0.887719298245614 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[8] <= 3.0575000047683716) {
                        if (input[0] <= 2.905500054359436) {
                            if (input[4] <= 0.3489999920129776) {
                                if (input[0] <= 2.896499991416931) {
                                    double tempArray[2] = { 0.004553734061930784, 0.9954462659380692 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.10037174721189591, 0.8996282527881041 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.893999934196472) {
                                    double tempArray[2] = { 0.15384615384615385, 0.8461538461538461 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9242424242424242, 0.07575757575757576 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= 0.06549999862909317) {
                                if (input[8] <= 2.9884999990463257) {
                                    double tempArray[2] = { 0.0005647803004631198, 0.9994352196995369 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.007751937984496124, 0.9922480620155039 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.07850000262260437) {
                                    double tempArray[2] = { 0.050189393939393936, 0.9498106060606061 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0009655616350177019, 0.9990344383649823 };
                                    memcpy(var10, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.7549999952316284) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var10, tempArray, 2 * sizeof(double));
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var10, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
        }
    }
    add_vectors(var7, var10, 2, var6);
    double var11[2];
    if (input[5] <= 4.888499975204468) {
        if (input[0] <= 2.3554999828338623) {
            if (input[7] <= 2.5325000286102295) {
                if (input[1] <= 0.05250000022351742) {
                    if (input[0] <= 2.070499897003174) {
                        if (input[0] <= 2.0394999980926514) {
                            if (input[6] <= 4.187000036239624) {
                                if (input[6] <= 3.8545000553131104) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9735744089012517, 0.02642559109874826 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 1.9950000047683716) {
                                    double tempArray[2] = { 0.40397350993377484, 0.5960264900662252 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.0399999618530273) {
                                if (input[8] <= 1.690500020980835) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7966101694915254, 0.2033898305084746 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.403499960899353) {
                                    double tempArray[2] = { 0.9860834990059643, 0.013916500994035786 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7439024390243902, 0.25609756097560976 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[8] <= 1.7639999985694885) {
                            if (input[0] <= 2.2695000171661377) {
                                if (input[0] <= 2.0929999351501465) {
                                    double tempArray[2] = { 0.7017543859649122, 0.2982456140350877 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9887920298879203, 0.0112079701120797 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 2.294000029563904) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5894039735099338, 0.4105960264900662 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 3.5770000219345093) {
                                if (input[6] <= 2.0994999408721924) {
                                    double tempArray[2] = { 0.8686708860759493, 0.13132911392405064 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6690837178642056, 0.3309162821357943 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.017500000074505806) {
                                    double tempArray[2] = { 0.547244094488189, 0.452755905511811 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.176, 0.824 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 2.6994999647140503) {
                        if (input[0] <= 2.070499897003174) {
                            if (input[7] <= 2.0089999437332153) {
                                if (input[0] <= 2.0175000429153442) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.25925925925925924, 0.7407407407407407 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.340000033378601) {
                                if (input[8] <= 2.226499915122986) {
                                    double tempArray[2] = { 0.8622198879551821, 0.13778011204481794 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6305732484076433, 0.36942675159235666 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.1394999027252197) {
                                    double tempArray[2] = { 0.8737864077669902, 0.1262135922330097 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.48392857142857143, 0.5160714285714286 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.0394999980926514) {
                            if (input[1] <= 0.057499999180436134) {
                                if (input[7] <= 2.3660000562667847) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9649122807017544, 0.03508771929824561 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[6] <= 4.0970001220703125) {
                                if (input[8] <= 1.7735000252723694) {
                                    double tempArray[2] = { 0.5384615384615384, 0.46153846153846156 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9038611039129308, 0.09613889608706919 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.07899999991059303) {
                                    double tempArray[2] = { 0.8888888888888888, 0.1111111111111111 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0625, 0.9375 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[2] <= -0.007500000298023224) {
                    if (input[6] <= 3.7200000286102295) {
                        if (input[1] <= -0.08349999785423279) {
                            if (input[0] <= 2.242000102996826) {
                                if (input[5] <= 4.725500106811523) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.149999976158142) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[7] <= 2.806999921798706) {
                                    double tempArray[2] = { 0.7568306010928961, 0.24316939890710382 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.24468085106382978, 0.7553191489361702 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.871000051498413) {
                            if (input[7] <= 2.591499924659729) {
                                if (input[2] <= -0.09950000047683716) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.09090909090909091, 0.9090909090909091 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 4.0945000648498535) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5977011494252874, 0.40229885057471265 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 4.665999889373779) {
                                if (input[3] <= -0.14400000125169754) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.046875, 0.953125 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.1679999828338623) {
                                    double tempArray[2] = { 0.13333333333333333, 0.8666666666666667 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7777777777777778, 0.2222222222222222 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= 0.019499999471008778) {
                        if (input[2] <= 0.27650000154972076) {
                            if (input[0] <= 2.1054999828338623) {
                                if (input[5] <= 4.8374998569488525) {
                                    double tempArray[2] = { 0.983739837398374, 0.016260162601626018 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8292682926829268, 0.17073170731707318 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.7890000343322754) {
                                    double tempArray[2] = { 0.4296577946768061, 0.5703422053231939 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9102564102564102, 0.08974358974358974 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.075500011444092) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[0] <= 2.270500063896179) {
                                    double tempArray[2] = { 0.8775510204081632, 0.12244897959183673 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[1] <= 0.18949999660253525) {
                            if (input[8] <= 2.2165000438690186) {
                                if (input[2] <= 0.11050000041723251) {
                                    double tempArray[2] = { 0.9477693144722524, 0.05223068552774755 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.822429906542056, 0.17757009345794392 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.125) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8819444444444444, 0.11805555555555555 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= 0.5684999823570251) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 1.278499960899353) {
                                    double tempArray[2] = { 0.56, 0.44 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[6] <= 3.512500047683716) {
                if (input[0] <= 2.9809999465942383) {
                    if (input[6] <= 2.871500015258789) {
                        if (input[6] <= 1.881500005722046) {
                            if (input[3] <= 0.2435000017285347) {
                                if (input[6] <= 1.6355000138282776) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4248927038626609, 0.575107296137339 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.9289999902248383) {
                                    double tempArray[2] = { 0.9152046783625731, 0.0847953216374269 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.19298245614035087, 0.8070175438596491 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.3445000648498535) {
                                if (input[5] <= 2.3244999647140503) {
                                    double tempArray[2] = { 0.689161554192229, 0.310838445807771 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.30549510337323177, 0.6945048966267682 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.135500006377697) {
                                    double tempArray[2] = { 0.6818580192813322, 0.31814198071866784 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.48548812664907653, 0.5145118733509235 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.1344999074935913) {
                            if (input[7] <= 2.009999990463257) {
                                if (input[5] <= 3.083999991416931) {
                                    double tempArray[2] = { 0.7881136950904393, 0.21188630490956073 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2693409742120344, 0.7306590257879656 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.122499942779541) {
                                    double tempArray[2] = { 0.8960055096418733, 0.10399449035812672 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6890243902439024, 0.31097560975609756 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= -0.007500000298023224) {
                                if (input[8] <= 2.10699999332428) {
                                    double tempArray[2] = { 0.4074074074074074, 0.5925925925925926 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.10857142857142857, 0.8914285714285715 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.033000001683831215) {
                                    double tempArray[2] = { 0.07766990291262135, 0.9223300970873787 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6697574893009985, 0.3302425106990014 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 2.4714999198913574) {
                        if (input[0] <= 3.027500033378601) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var11, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[3] <= 0.44050000607967377) {
                                if (input[0] <= 3.2195000648498535) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.10997442455242967, 0.8900255754475703 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.4429999589920044) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var11, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[3] <= 0.135500006377697) {
                                if (input[7] <= 2.5779999494552612) {
                                    double tempArray[2] = { 0.2465753424657534, 0.7534246575342466 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.2604999989271164) {
                                    double tempArray[2] = { 0.010135135135135136, 0.9898648648648649 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.12195121951219512, 0.8780487804878049 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[6] <= 5.3420000076293945) {
                    if (input[0] <= 2.5329999923706055) {
                        if (input[4] <= 0.07499999925494194) {
                            if (input[3] <= 0.20549999922513962) {
                                if (input[6] <= 4.581000089645386) {
                                    double tempArray[2] = { 0.5722070844686649, 0.42779291553133514 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.034722222222222224, 0.9652777777777778 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.3760000467300415) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[8] <= 1.972000002861023) {
                                    double tempArray[2] = { 0.7, 0.3 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9922779922779923, 0.007722007722007722 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 2.8644999265670776) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var11, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[4] <= 0.5925000011920929) {
                                if (input[5] <= 2.9499999284744263) {
                                    double tempArray[2] = { 0.8048780487804879, 0.1951219512195122 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.06940170940170941, 0.9305982905982906 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.971500039100647) {
                                    double tempArray[2] = { 0.8854166666666666, 0.11458333333333333 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.430999994277954) {
                        if (input[3] <= 0.3995000123977661) {
                            if (input[3] <= 0.153499998152256) {
                                if (input[0] <= 2.4200000762939453) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var11, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[6] <= 7.7204999923706055) {
                            if (input[3] <= 0.39249999821186066) {
                                if (input[1] <= 0.09349999949336052) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.002, 0.998 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.09449999779462814) {
                                    double tempArray[2] = { 0.6451612903225806, 0.3548387096774194 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.037037037037037035, 0.9629629629629629 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var11, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
        }
    }
    else {
        if (input[8] <= 2.0225000381469727) {
            if (input[7] <= 2.7209999561309814) {
                if (input[0] <= 2.277999997138977) {
                    if (input[0] <= 2.1654999256134033) {
                        if (input[5] <= 8.76200008392334) {
                            if (input[1] <= 0.02850000001490116) {
                                if (input[3] <= 0.2735000103712082) {
                                    double tempArray[2] = { 0.9811550151975684, 0.01884498480243161 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7380952380952381, 0.2619047619047619 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.0484999418258667) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.992989165073295, 0.007010834926704908 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var11, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[5] <= 5.310499906539917) {
                            if (input[8] <= 1.6460000276565552) {
                                if (input[2] <= 0.14549999684095383) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8125, 0.1875 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.291499972343445) {
                                    double tempArray[2] = { 0.10552763819095477, 0.8944723618090452 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7262247838616714, 0.2737752161383285 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 1.7595000267028809) {
                                if (input[6] <= 2.7144999504089355) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9701492537313433, 0.029850746268656716 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.003999999957159162) {
                                    double tempArray[2] = { 0.4533333333333333, 0.5466666666666666 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9462025316455697, 0.05379746835443038 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= -0.022499999962747097) {
                        if (input[2] <= 4.607000112533569) {
                            if (input[8] <= 0.5) {
                                if (input[0] <= 2.9010000228881836) {
                                    double tempArray[2] = { 0.9659090909090909, 0.03409090909090909 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.7619999647140503) {
                                    double tempArray[2] = { 0.5255673222390318, 0.4744326777609682 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.013048016701461378, 0.9869519832985386 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var11, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[5] <= 5.421499967575073) {
                            if (input[3] <= 0.08650000020861626) {
                                if (input[8] <= 0.7854999899864197) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.11383928571428571, 0.8861607142857143 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.3899999856948853) {
                                    double tempArray[2] = { 0.6428571428571429, 0.35714285714285715 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.16129032258064516, 0.8387096774193549 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= 0.024500000290572643) {
                                if (input[2] <= -0.05949999950826168) {
                                    double tempArray[2] = { 0.5696022727272727, 0.4303977272727273 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.36857280153772226, 0.6314271984622778 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.753499984741211) {
                                    double tempArray[2] = { 0.7431363424848767, 0.2568636575151233 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.016, 0.984 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[7] <= 2.8209999799728394) {
                    if (input[0] <= 2.4264999628067017) {
                        if (input[6] <= 3.2144999504089355) {
                            if (input[0] <= 2.365000009536743) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[7] <= 2.7384999990463257) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9054054054054054, 0.0945945945945946 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.7350000143051147) {
                                if (input[0] <= 2.2549999952316284) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.26666666666666666, 0.7333333333333333 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= -0.05949999950826168) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9190535491905355, 0.08094645080946451 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= 0.523499995470047) {
                            if (input[0] <= 2.7619999647140503) {
                                if (input[5] <= 7.552500009536743) {
                                    double tempArray[2] = { 0.09278350515463918, 0.9072164948453608 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.7690000534057617) {
                                if (input[0] <= 2.6750000715255737) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 3.9359999895095825) {
                        if (input[5] <= 5.563499927520752) {
                            if (input[8] <= 1.9455000162124634) {
                                if (input[0] <= 2.4220000505447388) {
                                    double tempArray[2] = { 0.87, 0.13 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.08000000193715096) {
                                    double tempArray[2] = { 0.041666666666666664, 0.9583333333333334 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.39893617021276595, 0.601063829787234 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.6945000886917114) {
                                if (input[8] <= 1.9135000109672546) {
                                    double tempArray[2] = { 0.3076923076923077, 0.6923076923076923 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9393939393939394, 0.06060606060606061 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.833999991416931) {
                            if (input[6] <= 4.029999852180481) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[4] <= -0.027999999932944775) {
                                if (input[1] <= -0.0010000020265579224) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[1] <= 0.0035000001080334187) {
                if (input[0] <= 2.4609999656677246) {
                    if (input[0] <= 2.177000045776367) {
                        if (input[7] <= 2.215499997138977) {
                            if (input[6] <= 4.436000108718872) {
                                if (input[2] <= 0.21400000154972076) {
                                    double tempArray[2] = { 0.1, 0.9 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.247499942779541) {
                                    double tempArray[2] = { 0.5, 0.5 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9733333333333334, 0.02666666666666667 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 3.1644999980926514) {
                                if (input[5] <= 6.627000093460083) {
                                    double tempArray[2] = { 0.7247838616714697, 0.27521613832853026 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9659090909090909, 0.03409090909090909 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.044000029563904) {
                                    double tempArray[2] = { 0.5476190476190477, 0.4523809523809524 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9792772302249179, 0.020722769775082132 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.26749999821186066) {
                            if (input[8] <= 2.671500086784363) {
                                if (input[0] <= 2.3785001039505005) {
                                    double tempArray[2] = { 0.6224148118967895, 0.37758518810321057 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4371584699453552, 0.5628415300546448 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.8450000286102295) {
                                    double tempArray[2] = { 0.8807947019867549, 0.11920529801324503 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9937888198757764, 0.006211180124223602 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.3774999380111694) {
                                if (input[7] <= 2.6399999856948853) {
                                    double tempArray[2] = { 0.6620689655172414, 0.33793103448275863 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2222222222222222, 0.7777777777777778 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.02049999963492155) {
                                    double tempArray[2] = { 0.9831081081081081, 0.016891891891891893 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.868663594470046, 0.1313364055299539 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 3.1614999771118164) {
                        if (input[7] <= 2.652000069618225) {
                            if (input[4] <= 0.34049999713897705) {
                                if (input[0] <= 2.8100000619888306) {
                                    double tempArray[2] = { 0.4365924491771539, 0.5634075508228461 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.03555045871559633, 0.9644495412844036 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.6815000176429749) {
                                    double tempArray[2] = { 0.5661375661375662, 0.43386243386243384 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.014925373134328358, 0.9850746268656716 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= 0.011500000022351742) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 3.069999933242798) {
                                    double tempArray[2] = { 0.10337552742616034, 0.8966244725738397 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4, 0.6 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[8] <= 2.041000008583069) {
                            if (input[7] <= 2.5075000524520874) {
                                if (input[4] <= 0.02850000001490116) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9166666666666666, 0.08333333333333333 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.6204999685287476) {
                                    double tempArray[2] = { 0.4444444444444444, 0.5555555555555556 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= -0.13849999755620956) {
                                if (input[6] <= 5.130499839782715) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 8.496500015258789) {
                                    double tempArray[2] = { 0.08897893714912836, 0.9110210628508716 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.013252555850056797, 0.9867474441499432 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[8] <= 2.1864999532699585) {
                    if (input[6] <= 4.581499814987183) {
                        if (input[2] <= -0.12049999833106995) {
                            if (input[3] <= 0.14749999716877937) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[3] <= -0.03150000050663948) {
                                if (input[6] <= 4.2225000858306885) {
                                    double tempArray[2] = { 0.9295774647887324, 0.07042253521126761 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.533500075340271) {
                                    double tempArray[2] = { 0.8115727002967359, 0.1884272997032641 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1062015503875969, 0.8937984496124031 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.09050000086426735) {
                            if (input[8] <= 2.1234999895095825) {
                                if (input[0] <= 2.447000026702881) {
                                    double tempArray[2] = { 0.8653846153846154, 0.1346153846153846 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.07407407407407407, 0.9259259259259259 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[3] <= 0.17099999636411667) {
                                if (input[0] <= 2.7640000581741333) {
                                    double tempArray[2] = { 0.8055555555555556, 0.19444444444444445 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var11, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 4.728500127792358) {
                        if (input[6] <= 3.3954999446868896) {
                            if (input[2] <= -0.05949999950826168) {
                                if (input[0] <= 2.833000063896179) {
                                    double tempArray[2] = { 0.8752783964365256, 0.12472160356347439 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.08888888888888889, 0.9111111111111111 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.28550000488758087) {
                                    double tempArray[2] = { 0.6180348040077343, 0.38196519599226575 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3811074918566775, 0.6188925081433225 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.7415000200271606) {
                                if (input[4] <= 0.1655000001192093) {
                                    double tempArray[2] = { 0.7630188679245283, 0.2369811320754717 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.955625, 0.044375 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.6424999237060547) {
                                    double tempArray[2] = { 0.0911983032873807, 0.9088016967126193 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.02054794520547945, 0.9794520547945206 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.385499954223633) {
                            if (input[0] <= 2.0959999561309814) {
                                if (input[3] <= -0.019499999471008778) {
                                    double tempArray[2] = { 0.9351351351351351, 0.06486486486486487 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.24899999052286148) {
                                    double tempArray[2] = { 0.6921723834652594, 0.30782761653474056 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9464285714285714, 0.05357142857142857 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.896499991416931) {
                                if (input[7] <= 2.9200000762939453) {
                                    double tempArray[2] = { 0.2651215375918598, 0.7348784624081401 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.5115000009536743) {
                                    double tempArray[2] = { 0.01069937369519833, 0.9893006263048016 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3709677419354839, 0.6290322580645161 };
                                    memcpy(var11, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    add_vectors(var6, var11, 2, var5);
    double var12[2];
    if (input[6] <= 3.549499988555908) {
        if (input[5] <= 6.019500017166138) {
            if (input[7] <= 2.4774999618530273) {
                if (input[6] <= 1.784500002861023) {
                    if (input[0] <= 2.91949999332428) {
                        if (input[8] <= 2.1950000524520874) {
                            if (input[0] <= 2.1649999618530273) {
                                if (input[0] <= 2.0175000429153442) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9739336492890995, 0.026066350710900472 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.1804999113082886) {
                                    double tempArray[2] = { 0.8733654507914659, 0.12663454920853406 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9704280155642023, 0.029571984435797664 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.13599999994039536) {
                                if (input[0] <= 2.1169999837875366) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.06451612903225806, 0.9354838709677419 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.029500000178813934) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= 0.2905000001192093) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var12, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[6] <= 1.6380000114440918) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[8] <= 2.093000054359436) {
                                    double tempArray[2] = { 0.17307692307692307, 0.8269230769230769 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9803921568627451, 0.0196078431372549 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[8] <= 1.7960000038146973) {
                        if (input[1] <= -0.14649999886751175) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var12, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[5] <= 0.8510000109672546) {
                                if (input[8] <= 1.7735000252723694) {
                                    double tempArray[2] = { 0.12195121951219512, 0.8780487804878049 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.6385000944137573) {
                                    double tempArray[2] = { 0.968741038141669, 0.03125896185833094 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.012987012987012988, 0.987012987012987 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.215499997138977) {
                            if (input[0] <= 2.059499979019165) {
                                if (input[6] <= 2.5824999809265137) {
                                    double tempArray[2] = { 0.9971479500891266, 0.00285204991087344 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9999627643729521, 0.00003723562704795949 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 2.752500057220459) {
                                    double tempArray[2] = { 0.6525704809286899, 0.3474295190713101 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8909287257019438, 0.10907127429805616 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.1234999895095825) {
                                if (input[7] <= 2.0295000076293945) {
                                    double tempArray[2] = { 0.46831530139103555, 0.5316846986089645 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6884646628757108, 0.3115353371242892 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.503499984741211) {
                                    double tempArray[2] = { 0.6640949554896143, 0.3359050445103858 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3315771583956492, 0.6684228416043508 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[1] <= 0.07349999994039536) {
                    if (input[0] <= 2.656999945640564) {
                        if (input[1] <= 0.00849999999627471) {
                            if (input[0] <= 2.1320000886917114) {
                                if (input[0] <= 2.0554999113082886) {
                                    double tempArray[2] = { 0.9968096419709322, 0.0031903580290677065 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8260869565217391, 0.17391304347826086 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 2.2829999923706055) {
                                    double tempArray[2] = { 0.4604651162790698, 0.5395348837209303 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.651564185544768, 0.3484358144552319 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.0959999561309814) {
                                if (input[0] <= 2.0649999380111694) {
                                    double tempArray[2] = { 0.9967558799675588, 0.0032441200324412004 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8670886075949367, 0.13291139240506328 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.12150000035762787) {
                                    double tempArray[2] = { 0.7858439201451906, 0.21415607985480944 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5846833578792342, 0.41531664212076586 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.676500082015991) {
                            if (input[8] <= 2.356500029563904) {
                                if (input[2] <= 0.41599999368190765) {
                                    double tempArray[2] = { 0.07643312101910828, 0.9235668789808917 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.47928994082840237, 0.5207100591715976 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.5609999895095825) {
                                    double tempArray[2] = { 0.3886138613861386, 0.6113861386138614 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.07799442896935933, 0.9220055710306406 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.7445000410079956) {
                                if (input[4] <= 1.1829999685287476) {
                                    double tempArray[2] = { 0.1643835616438356, 0.8356164383561644 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.8240000009536743) {
                                    double tempArray[2] = { 0.0028835063437139563, 0.9971164936562861 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0975609756097561, 0.9024390243902439 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[7] <= 2.8019999265670776) {
                        if (input[6] <= 1.8264999985694885) {
                            if (input[5] <= 4.102499961853027) {
                                if (input[0] <= 2.6395000219345093) {
                                    double tempArray[2] = { 0.8802816901408451, 0.11971830985915492 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.26749999821186066) {
                                    double tempArray[2] = { 0.9853911404335532, 0.014608859566446749 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6637931034482759, 0.33620689655172414 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= 0.3774999976158142) {
                                if (input[2] <= 0.4675000011920929) {
                                    double tempArray[2] = { 0.8265011547344111, 0.17349884526558892 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5265700483091788, 0.47342995169082125 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.5135000050067902) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3375, 0.6625 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 2.249500036239624) {
                            if (input[5] <= 2.8865000009536743) {
                                if (input[8] <= 2.065999984741211) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 1.121500015258789) {
                                    double tempArray[2] = { 0.2657894736842105, 0.7342105263157894 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5076923076923077, 0.49230769230769234 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= -0.0020000000367872417) {
                                if (input[1] <= 0.17750000208616257) {
                                    double tempArray[2] = { 0.4117647058823529, 0.5882352941176471 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.38750000298023224) {
                                    double tempArray[2] = { 0.8509433962264151, 0.1490566037735849 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.125, 0.875 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[0] <= 2.614500045776367) {
                if (input[6] <= 2.6535000801086426) {
                    if (input[2] <= 0.2874999940395355) {
                        if (input[8] <= 2.083500027656555) {
                            if (input[2] <= 0.27950000762939453) {
                                if (input[5] <= 6.417500019073486) {
                                    double tempArray[2] = { 0.9258893280632411, 0.0741106719367589 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9718456725755996, 0.028154327424400417 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.371500015258789) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1, 0.9 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.1389999389648438) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[5] <= 7.488499879837036) {
                                    double tempArray[2] = { 0.6675938803894298, 0.33240611961057026 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9393939393939394, 0.06060606060606061 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.7045000791549683) {
                            if (input[8] <= 2.100000023841858) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 1.8419999480247498) {
                                    double tempArray[2] = { 0.8636363636363636, 0.13636363636363635 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.7125000953674316) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.1644999980926514) {
                        if (input[3] <= 0.27250000834465027) {
                            if (input[0] <= 2.1109999418258667) {
                                if (input[5] <= 6.079499959945679) {
                                    double tempArray[2] = { 0.9827586206896551, 0.017241379310344827 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.044999999925494194) {
                                    double tempArray[2] = { 0.48484848484848486, 0.5151515151515151 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9766081871345029, 0.023391812865497075 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 2.79449999332428) {
                                if (input[8] <= 2.197000026702881) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[6] <= 3.2899999618530273) {
                            if (input[7] <= 2.7519999742507935) {
                                if (input[7] <= 2.6344999074935913) {
                                    double tempArray[2] = { 0.4889397406559878, 0.5110602593440122 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.20964360587002095, 0.790356394129979 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.5590000003576279) {
                                    double tempArray[2] = { 0.9428571428571428, 0.05714285714285714 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.56, 0.44 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= -0.042000001296401024) {
                                if (input[3] <= -0.005500000203028321) {
                                    double tempArray[2] = { 0.875, 0.125 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4669811320754717, 0.5330188679245284 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.384500026702881) {
                                    double tempArray[2] = { 0.907563025210084, 0.09243697478991597 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6644736842105263, 0.3355263157894737 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[6] <= 1.56850004196167) {
                    if (input[2] <= 4.656000137329102) {
                        if (input[5] <= 7.15149998664856) {
                            if (input[6] <= 1.5390000343322754) {
                                if (input[2] <= 2.101499915122986) {
                                    double tempArray[2] = { 0.4397590361445783, 0.5602409638554217 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.437999963760376) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.13199999928474426) {
                                if (input[0] <= 2.7690000534057617) {
                                    double tempArray[2] = { 0.8802177858439202, 0.11978221415607986 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5592105263157895, 0.4407894736842105 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 1.7665000557899475) {
                                    double tempArray[2] = { 0.03225806451612903, 0.967741935483871 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3, 0.7 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        double tempArray[2] = { 0.0, 1.0 };
                        memcpy(var12, tempArray, 2 * sizeof(double));
                    }
                }
                else {
                    if (input[0] <= 2.8289999961853027) {
                        if (input[8] <= 2.9744999408721924) {
                            if (input[4] <= 0.06650000065565109) {
                                if (input[1] <= 0.0025000000605359674) {
                                    double tempArray[2] = { 0.18840579710144928, 0.8115942028985508 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6490066225165563, 0.3509933774834437 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.7345000505447388) {
                                    double tempArray[2] = { 0.03429602888086643, 0.9657039711191335 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.30991735537190085, 0.6900826446280992 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= 0.09750000014901161) {
                                if (input[3] <= 0.011500000022351742) {
                                    double tempArray[2] = { 0.12727272727272726, 0.8727272727272727 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8057851239669421, 0.19421487603305784 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 7.4100000858306885) {
                                    double tempArray[2] = { 0.13286713286713286, 0.8671328671328671 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7209302325581395, 0.27906976744186046 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 6.452500104904175) {
                            if (input[1] <= -0.05550000071525574) {
                                if (input[6] <= 2.909999966621399) {
                                    double tempArray[2] = { 0.53125, 0.46875 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.07142857142857142, 0.9285714285714286 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.11549999937415123) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.01818181818181818, 0.9818181818181818 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.112499952316284) {
                                if (input[5] <= 7.570500135421753) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.08249999955296516) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0013351134846461949, 0.9986648865153538 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        if (input[8] <= 2.013000011444092) {
            if (input[6] <= 4.601500034332275) {
                if (input[0] <= 2.3669999837875366) {
                    if (input[0] <= 2.149999976158142) {
                        if (input[1] <= -0.09449999779462814) {
                            if (input[0] <= 2.0934998989105225) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[5] <= 4.425500154495239) {
                                if (input[3] <= 0.10050000250339508) {
                                    double tempArray[2] = { 0.9314720812182741, 0.06852791878172589 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.252500057220459) {
                                    double tempArray[2] = { 0.9929453262786596, 0.007054673721340388 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 5.588500022888184) {
                            if (input[1] <= 0.02350000012665987) {
                                if (input[0] <= 2.305999994277954) {
                                    double tempArray[2] = { 0.1331168831168831, 0.8668831168831169 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8142857142857143, 0.18571428571428572 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.7639999985694885) {
                                    double tempArray[2] = { 0.3492063492063492, 0.6507936507936508 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8381742738589212, 0.16182572614107885 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.274999976158142) {
                                if (input[1] <= 0.03649999923072755) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.343000054359436) {
                                    double tempArray[2] = { 0.9959595959595959, 0.00404040404040404 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6976744186046512, 0.3023255813953488 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[5] <= 6.023499965667725) {
                        if (input[1] <= 0.10050000250339508) {
                            if (input[0] <= 2.3734999895095825) {
                                if (input[7] <= 2.4884999990463257) {
                                    double tempArray[2] = { 0.34285714285714286, 0.6571428571428571 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.7489999532699585) {
                                    double tempArray[2] = { 0.2328830926874709, 0.7671169073125291 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= 0.24400000274181366) {
                                if (input[7] <= 2.1160000562667847) {
                                    double tempArray[2] = { 0.26666666666666666, 0.7333333333333333 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8384615384615385, 0.16153846153846155 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.8634999990463257) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= 0.00849999999627471) {
                            if (input[8] <= 1.9424999952316284) {
                                if (input[4] <= 0.05849999934434891) {
                                    double tempArray[2] = { 0.7, 0.3 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.935064935064935, 0.06493506493506493 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.9704999923706055) {
                                    double tempArray[2] = { 0.5029940119760479, 0.49700598802395207 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.24007561436672967, 0.7599243856332704 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.007500000298023224) {
                                if (input[2] <= 0.21299999952316284) {
                                    double tempArray[2] = { 0.1036745406824147, 0.8963254593175853 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5705128205128205, 0.42948717948717946 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.468999981880188) {
                                    double tempArray[2] = { 0.8688524590163934, 0.13114754098360656 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.47107438016528924, 0.5289256198347108 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.385499954223633) {
                    if (input[0] <= 2.268499970436096) {
                        double tempArray[2] = { 1.0, 0.0 };
                        memcpy(var12, tempArray, 2 * sizeof(double));
                    }
                    else {
                        if (input[1] <= -0.07349999994039536) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var12, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[2] <= 0.09699999913573265) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= -0.07349999994039536) {
                        if (input[7] <= 2.6214998960494995) {
                            if (input[3] <= 0.0015000002458691597) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.591999888420105) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[6] <= 6.480499982833862) {
                            if (input[7] <= 2.490499973297119) {
                                if (input[6] <= 4.605499982833862) {
                                    double tempArray[2] = { 0.36363636363636365, 0.6363636363636364 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.01901743264659271, 0.9809825673534073 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.5195000171661377) {
                                    double tempArray[2] = { 0.5892857142857143, 0.4107142857142857 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.13654618473895583, 0.8634538152610441 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var12, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
        }
        else {
            if (input[5] <= 3.4954999685287476) {
                if (input[3] <= 0.013999999966472387) {
                    if (input[6] <= 4.7769999504089355) {
                        if (input[2] <= 0.045499999076128006) {
                            if (input[8] <= 2.452500104904175) {
                                if (input[0] <= 2.4980000257492065) {
                                    double tempArray[2] = { 0.6794871794871795, 0.32051282051282054 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[3] <= 0.010999999940395355) {
                                if (input[0] <= 2.8315000534057617) {
                                    double tempArray[2] = { 0.9802371541501976, 0.019762845849802372 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[6] <= 6.1545000076293945) {
                            if (input[6] <= 5.228500127792358) {
                                if (input[0] <= 2.319000005722046) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= -0.05949999950826168) {
                                    double tempArray[2] = { 0.1, 0.9 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.506999969482422) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= 0.07149999961256981) {
                        if (input[2] <= 0.0625) {
                            if (input[6] <= 4.796500205993652) {
                                if (input[0] <= 2.0199999809265137) {
                                    double tempArray[2] = { 0.9712230215827338, 0.02877697841726619 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.13402061855670103, 0.865979381443299 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.24249999970197678) {
                                    double tempArray[2] = { 0.09803921568627451, 0.9019607843137255 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.5509999990463257) {
                                if (input[8] <= 2.084999918937683) {
                                    double tempArray[2] = { 0.5945945945945946, 0.40540540540540543 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9829867674858223, 0.017013232514177693 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.28450000286102295) {
                            if (input[8] <= 2.3700000047683716) {
                                if (input[5] <= 2.7195000648498535) {
                                    double tempArray[2] = { 0.8607318405243036, 0.13926815947569635 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9732664995822891, 0.026733500417710943 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.031000018119812) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7631578947368421, 0.23684210526315788 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var12, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
            else {
                if (input[6] <= 7.250499963760376) {
                    if (input[4] <= 0.006500000134110451) {
                        if (input[0] <= 2.3669999837875366) {
                            if (input[8] <= 2.361999988555908) {
                                if (input[5] <= 4.35699987411499) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9102990033222591, 0.08970099667774087 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.2855000495910645) {
                                    double tempArray[2] = { 0.9927652733118971, 0.007234726688102894 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5283540802213001, 0.47164591977869985 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= -0.14750000089406967) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[5] <= 8.148499965667725) {
                                    double tempArray[2] = { 0.13734335839598996, 0.86265664160401 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.03861003861003861, 0.9613899613899614 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[1] <= -0.02049999963492155) {
                            if (input[0] <= 2.4609999656677246) {
                                if (input[0] <= 2.177000045776367) {
                                    double tempArray[2] = { 0.9427038626609442, 0.057296137339055794 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5268518518518519, 0.47314814814814815 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.0425000190734863) {
                                    double tempArray[2] = { 0.49310344827586206, 0.506896551724138 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.07578111513896081, 0.9242188848610392 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.1054999977350235) {
                                if (input[5] <= 8.55049991607666) {
                                    double tempArray[2] = { 0.3909062898895198, 0.6090937101104802 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.16005665722379603, 0.839943342776204 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.127500057220459) {
                                    double tempArray[2] = { 0.7857142857142857, 0.21428571428571427 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.49515669515669514, 0.5048433048433049 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= -0.0755000002682209) {
                        if (input[1] <= -0.16849999874830246) {
                            if (input[0] <= 2.6480000019073486) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.8075000047683716) {
                                if (input[7] <= 2.502500057220459) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.99822695035461, 0.0017730496453900709 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.45550000667572) {
                            if (input[5] <= 7.051500082015991) {
                                if (input[2] <= -0.07750000059604645) {
                                    double tempArray[2] = { 0.7714285714285715, 0.22857142857142856 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.17777777777777778, 0.8222222222222222 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.677500009536743) {
                                    double tempArray[2] = { 0.8780487804878049, 0.12195121951219512 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.896499991416931) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var12, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[3] <= 0.27400000393390656) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var12, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    add_vectors(var5, var12, 2, var4);
    double var13[2];
    if (input[6] <= 3.6085000038146973) {
        if (input[0] <= 2.503499984741211) {
            if (input[1] <= 0.025500000454485416) {
                if (input[6] <= 1.6704999804496765) {
                    if (input[3] <= 0.19049999862909317) {
                        if (input[4] <= 0.4829999953508377) {
                            if (input[8] <= 1.899999976158142) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[0] <= 2.274999976158142) {
                                    double tempArray[2] = { 0.9936708860759493, 0.006329113924050633 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7971014492753623, 0.2028985507246377 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= 0.6245000064373016) {
                                if (input[2] <= 0.1849999949336052) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[5] <= 4.36050009727478) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var13, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[6] <= 1.237500011920929) {
                                if (input[0] <= 2.416000008583069) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6, 0.4 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.5544999837875366) {
                                    double tempArray[2] = { 0.037037037037037035, 0.9629629629629629 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[7] <= 2.458500027656555) {
                        if (input[0] <= 2.215499997138977) {
                            if (input[0] <= 2.070499897003174) {
                                if (input[0] <= 2.0554999113082886) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9486486486486486, 0.051351351351351354 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.09349999949336052) {
                                    double tempArray[2] = { 0.24347826086956523, 0.7565217391304347 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9031872509960159, 0.09681274900398407 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 2.553499937057495) {
                                if (input[5] <= 4.976000070571899) {
                                    double tempArray[2] = { 0.7338638373121131, 0.2661361626878868 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9021739130434783, 0.09782608695652174 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.450000047683716) {
                                    double tempArray[2] = { 0.5009885330170027, 0.49901146698299725 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8839009287925697, 0.11609907120743033 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.1320000886917114) {
                            if (input[7] <= 2.5735000371932983) {
                                if (input[3] <= 0.29350000619888306) {
                                    double tempArray[2] = { 0.9486301369863014, 0.05136986301369863 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5416666666666666, 0.4583333333333333 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 1.634499967098236) {
                                    double tempArray[2] = { 0.9239130434782609, 0.07608695652173914 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9974076474400518, 0.002592352559948153 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.3214999437332153) {
                                if (input[8] <= 1.5855000019073486) {
                                    double tempArray[2] = { 0.08, 0.92 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5671169269734065, 0.4328830730265935 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.746000051498413) {
                                    double tempArray[2] = { 0.7083128381701919, 0.29168716182980814 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9376947040498442, 0.06230529595015576 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[8] <= 1.7245000004768372) {
                    if (input[0] <= 2.0679999589920044) {
                        double tempArray[2] = { 1.0, 0.0 };
                        memcpy(var13, tempArray, 2 * sizeof(double));
                    }
                    else {
                        if (input[0] <= 2.1804999113082886) {
                            if (input[2] <= 1.4194999933242798) {
                                if (input[0] <= 2.166499972343445) {
                                    double tempArray[2] = { 0.9058823529411765, 0.09411764705882353 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2857142857142857, 0.7142857142857143 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[2] <= 0.24650000035762787) {
                                if (input[2] <= 0.14150000363588333) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9779411764705882, 0.022058823529411766 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.04150000028312206) {
                                    double tempArray[2] = { 0.09090909090909091, 0.9090909090909091 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9898477157360406, 0.01015228426395939 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 2.850499987602234) {
                        if (input[8] <= 2.114500045776367) {
                            if (input[8] <= 1.7300000190734863) {
                                if (input[3] <= 0.15600000321865082) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.070499897003174) {
                                    double tempArray[2] = { 0.9977969762419007, 0.002203023758099352 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8038513763732873, 0.19614862362671276 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.0649999380111694) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[2] <= -0.14350000023841858) {
                                    double tempArray[2] = { 0.9897959183673469, 0.01020408163265306 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7162004662004662, 0.2837995337995338 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.1414999961853027) {
                            if (input[4] <= -0.022499999962747097) {
                                if (input[7] <= 2.0190000534057617) {
                                    double tempArray[2] = { 0.425531914893617, 0.574468085106383 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9821428571428571, 0.017857142857142856 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.7300000190734863) {
                                    double tempArray[2] = { 0.8333333333333334, 0.16666666666666666 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9978028605893503, 0.002197139410649664 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.2760000228881836) {
                                if (input[7] <= 2.009999990463257) {
                                    double tempArray[2] = { 0.6870503597122302, 0.3129496402877698 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9087688219663419, 0.0912311780336581 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.9144999980926514) {
                                    double tempArray[2] = { 0.5529801324503312, 0.4470198675496689 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7734311328443357, 0.2265688671556642 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[5] <= 3.742500066757202) {
                if (input[5] <= 3.249500036239624) {
                    if (input[0] <= 2.9809999465942383) {
                        if (input[0] <= 2.765500068664551) {
                            if (input[1] <= 0.1405000016093254) {
                                if (input[3] <= 0.1744999960064888) {
                                    double tempArray[2] = { 0.35009548058561424, 0.6499045194143858 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.131496062992126, 0.868503937007874 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.194000005722046) {
                                    double tempArray[2] = { 0.5895741556534508, 0.4104258443465492 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.19101123595505617, 0.8089887640449438 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 3.0175000429153442) {
                                if (input[2] <= 0.5320000052452087) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.42857142857142855, 0.5714285714285714 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.009500000393018126) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.863013698630137, 0.136986301369863 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 3.027500033378601) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var13, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[6] <= 1.781499981880188) {
                                if (input[8] <= 2.0770000219345093) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8235294117647058, 0.17647058823529413 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 2.5149999856948853) {
                                    double tempArray[2] = { 0.16129032258064516, 0.8387096774193549 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[2] <= -0.00849999999627471) {
                        if (input[8] <= 2.1160000562667847) {
                            if (input[8] <= 2.0140000581741333) {
                                if (input[4] <= 0.015500000212341547) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8797814207650273, 0.12021857923497267 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.8984999656677246) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.7844998836517334) {
                                if (input[2] <= -0.03749999962747097) {
                                    double tempArray[2] = { 0.5172413793103449, 0.4827586206896552 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9811320754716981, 0.018867924528301886 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.1055000014603138) {
                                    double tempArray[2] = { 0.16666666666666666, 0.8333333333333334 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[3] <= 0.12250000238418579) {
                            if (input[7] <= 2.1959999799728394) {
                                if (input[0] <= 2.7204999923706055) {
                                    double tempArray[2] = { 0.7839195979899497, 0.21608040201005024 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.5549999475479126) {
                                    double tempArray[2] = { 0.022099447513812154, 0.9779005524861878 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2962962962962963, 0.7037037037037037 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.13700000196695328) {
                                if (input[4] <= 0.05949999950826168) {
                                    double tempArray[2] = { 0.6387665198237885, 0.36123348017621143 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1793400286944046, 0.8206599713055954 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.04649999924004078) {
                                    double tempArray[2] = { 0.21739130434782608, 0.782608695652174 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6005434782608695, 0.39945652173913043 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[8] <= 1.7735000252723694) {
                    if (input[5] <= 7.551500082015991) {
                        if (input[2] <= 0.0625) {
                            if (input[0] <= 2.7799999713897705) {
                                if (input[3] <= 0.12200000137090683) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6086956521739131, 0.391304347826087 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.7864999771118164) {
                                if (input[1] <= 0.04350000061094761) {
                                    double tempArray[2] = { 0.8292682926829268, 0.17073170731707318 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.45251396648044695, 0.547486033519553 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.8240000009536743) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.10476190476190476, 0.8952380952380953 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.541499972343445) {
                            if (input[5] <= 7.575000047683716) {
                                if (input[1] <= -0.12800000049173832) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5, 0.5 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.004000000073574483) {
                                    double tempArray[2] = { 0.09090909090909091, 0.9090909090909091 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= 0.08000000193715096) {
                                if (input[7] <= 2.64300000667572) {
                                    double tempArray[2] = { 0.9890909090909091, 0.01090909090909091 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7291666666666666, 0.2708333333333333 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.9524999856948853) {
                        if (input[7] <= 2.6239999532699585) {
                            if (input[8] <= 1.8359999656677246) {
                                if (input[6] <= 2.7660000324249268) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.01509433962264151, 0.9849056603773585 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.2044999599456787) {
                                    double tempArray[2] = { 0.4972793943695292, 0.5027206056304708 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.265379113018598, 0.734620886981402 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.614500045776367) {
                                if (input[7] <= 2.677999973297119) {
                                    double tempArray[2] = { 0.43478260869565216, 0.5652173913043478 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.3605000078678131) {
                                    double tempArray[2] = { 0.09186496956281129, 0.9081350304371887 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5135869565217391, 0.48641304347826086 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 2.3615000247955322) {
                            if (input[5] <= 4.874500036239624) {
                                if (input[0] <= 3.291999936103821) {
                                    double tempArray[2] = { 0.4588235294117647, 0.5411764705882353 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.135500006377697) {
                                    double tempArray[2] = { 0.0033783783783783786, 0.9966216216216216 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.09774436090225563, 0.9022556390977443 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 3.027500033378601) {
                                if (input[6] <= 2.3669999837875366) {
                                    double tempArray[2] = { 0.5, 0.5 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 3.0575000047683716) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.04394812680115274, 0.9560518731988472 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        if (input[1] <= 0.05949999950826168) {
            if (input[0] <= 2.38100004196167) {
                if (input[7] <= 2.831499934196472) {
                    if (input[4] <= -0.05949999950826168) {
                        if (input[5] <= 6.643500089645386) {
                            if (input[8] <= 2.362499952316284) {
                                if (input[1] <= 0.01000000024214387) {
                                    double tempArray[2] = { 0.29411764705882354, 0.7058823529411765 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.1480000019073486) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= -0.17599999904632568) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[5] <= 7.922500133514404) {
                                    double tempArray[2] = { 0.9111111111111111, 0.08888888888888889 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.582089552238806, 0.417910447761194 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.1545000076293945) {
                            if (input[8] <= 2.0169999599456787) {
                                if (input[1] <= -0.09750000014901161) {
                                    double tempArray[2] = { 0.9066339066339066, 0.09336609336609336 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9961984793917567, 0.0038015206082432974 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.2790000438690186) {
                                    double tempArray[2] = { 0.8253373313343328, 0.17466266866566715 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.981686541737649, 0.018313458262350937 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 5.259000062942505) {
                                if (input[6] <= 3.958500027656555) {
                                    double tempArray[2] = { 0.46601941747572817, 0.5339805825242718 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.11274509803921569, 0.8872549019607843 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.03849999979138374) {
                                    double tempArray[2] = { 0.7290322580645161, 0.2709677419354839 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8331524688008681, 0.16684753119913184 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 5.9019999504089355) {
                        if (input[4] <= 0.24799999594688416) {
                            if (input[2] <= 0.09349999949336052) {
                                if (input[1] <= -0.00849999999627471) {
                                    double tempArray[2] = { 0.46115288220551376, 0.5388471177944862 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7616666666666667, 0.23833333333333334 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.16499999910593033) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5789473684210527, 0.42105263157894735 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= -0.14150000363588333) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[8] <= 2.5665000677108765) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7959183673469388, 0.20408163265306123 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[7] <= 3.131500005722046) {
                            if (input[7] <= 3.097000002861023) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 7.884999990463257) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 7.871500253677368) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
            }
            else {
                if (input[4] <= 0.32250000536441803) {
                    if (input[0] <= 2.768999934196472) {
                        if (input[8] <= 2.531999945640564) {
                            if (input[7] <= 2.2695000171661377) {
                                if (input[8] <= 1.930999994277954) {
                                    double tempArray[2] = { 0.011594202898550725, 0.9884057971014493 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.22140672782874618, 0.7785932721712538 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.3945000171661377) {
                                    double tempArray[2] = { 0.5491898148148148, 0.4508101851851852 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.34244514106583074, 0.6575548589341693 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= -0.025500000454485416) {
                                if (input[8] <= 2.978999972343445) {
                                    double tempArray[2] = { 0.06824591088550479, 0.9317540891144952 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7777777777777778, 0.2222222222222222 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.916499972343445) {
                                    double tempArray[2] = { 0.2550693703308431, 0.7449306296691569 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9259259259259259, 0.07407407407407407 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.8499999046325684) {
                            if (input[6] <= 5.2804999351501465) {
                                if (input[4] <= 0.2939999997615814) {
                                    double tempArray[2] = { 0.05482775351770985, 0.9451722464822901 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7368421052631579, 0.2631578947368421 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.10999999940395355) {
                                    double tempArray[2] = { 0.8875, 0.1125 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.303886925795053, 0.696113074204947 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 3.0429999828338623) {
                                if (input[7] <= 2.472999930381775) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0037025008119519324, 0.9962974991880481 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 3.106500029563904) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.971500039100647) {
                        if (input[8] <= 2.753999948501587) {
                            if (input[8] <= 2.559000015258789) {
                                if (input[2] <= 0.04649999924004078) {
                                    double tempArray[2] = { 0.5690721649484536, 0.4309278350515464 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8189845474613686, 0.18101545253863136 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.8125) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8571428571428571, 0.14285714285714285 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= 0.2760000079870224) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[3] <= 0.22950000315904617) {
                            if (input[1] <= -0.10050000250339508) {
                                if (input[5] <= 6.890500068664551) {
                                    double tempArray[2] = { 0.22826086956521738, 0.7717391304347826 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.1875) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.03597122302158273, 0.9640287769784173 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var13, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
        }
        else {
            if (input[5] <= 3.2204999923706055) {
                if (input[8] <= 2.350499987602234) {
                    if (input[0] <= 2.462499976158142) {
                        if (input[8] <= 1.7319999933242798) {
                            if (input[0] <= 2.1589999198913574) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[6] <= 5.743499994277954) {
                                if (input[8] <= 2.309499979019165) {
                                    double tempArray[2] = { 0.9942112879884226, 0.005788712011577424 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7222222222222222, 0.2777777777777778 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.2435000017285347) {
                            if (input[0] <= 2.8634999990463257) {
                                if (input[5] <= 2.837000012397766) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8333333333333334, 0.16666666666666666 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[3] <= 0.3154999911785126) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 4.442000150680542) {
                        if (input[5] <= 2.753999948501587) {
                            if (input[0] <= 2.5080000162124634) {
                                if (input[8] <= 2.6304999589920044) {
                                    double tempArray[2] = { 0.7777777777777778, 0.2222222222222222 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[7] <= 2.166499972343445) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[2] <= 0.025500000454485416) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[3] <= -0.05599999986588955) {
                            if (input[3] <= -0.0645000021904707) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var13, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.6880000829696655) {
                    if (input[2] <= 0.018499999307096004) {
                        if (input[4] <= 0.2135000005364418) {
                            if (input[7] <= 2.978999972343445) {
                                if (input[2] <= -0.0494999997317791) {
                                    double tempArray[2] = { 0.7291325695581015, 0.27086743044189854 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8413936557462298, 0.15860634425377015 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.10500000044703484) {
                                    double tempArray[2] = { 0.975, 0.025 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2904761904761905, 0.7095238095238096 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.14300000667572) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[4] <= 0.3684999942779541) {
                                    double tempArray[2] = { 0.9754601226993865, 0.024539877300613498 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8463302752293578, 0.1536697247706422 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[1] <= 0.09149999916553497) {
                            if (input[5] <= 5.697499990463257) {
                                if (input[0] <= 2.152500033378601) {
                                    double tempArray[2] = { 0.9970760233918129, 0.0029239766081871343 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4954128440366973, 0.5045871559633027 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.56850004196167) {
                                    double tempArray[2] = { 0.9595015576323987, 0.040498442367601244 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.13333333333333333, 0.8666666666666667 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 6.523499965667725) {
                                if (input[1] <= 0.3490000069141388) {
                                    double tempArray[2] = { 0.9618705035971223, 0.038129496402877695 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5333333333333333, 0.4666666666666667 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 4.578500151634216) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[8] <= 2.940000057220459) {
                        if (input[3] <= 0.19749999791383743) {
                            if (input[3] <= -0.08749999850988388) {
                                if (input[0] <= 2.9234999418258667) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.884999990463257) {
                                    double tempArray[2] = { 0.054385964912280704, 0.9456140350877194 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0008116883116883117, 0.9991883116883117 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 4.317499876022339) {
                                if (input[6] <= 4.24150013923645) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.018867924528301886, 0.9811320754716981 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.05650000087916851) {
                                    double tempArray[2] = { 0.3838383838383838, 0.6161616161616161 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0711340206185567, 0.9288659793814433 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= 0.05949999950826168) {
                            if (input[3] <= 0.226500004529953) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var13, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[2] <= 0.04600000008940697) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.75, 0.25 };
                                    memcpy(var13, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var13, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
        }
    }
    add_vectors(var4, var13, 2, var3);
    double var14[2];
    if (input[6] <= 3.54449999332428) {
        if (input[0] <= 2.503499984741211) {
            if (input[7] <= 2.4854999780654907) {
                if (input[1] <= 0.060499999672174454) {
                    if (input[8] <= 1.7670000195503235) {
                        if (input[8] <= 1.5634999871253967) {
                            if (input[0] <= 2.166499972343445) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[0] <= 2.174499988555908) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= 0.3400000035762787) {
                                if (input[6] <= 2.652500033378601) {
                                    double tempArray[2] = { 0.9856601731601732, 0.01433982683982684 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.905341446923597, 0.09465855307640297 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.215499997138977) {
                            if (input[8] <= 2.197499990463257) {
                                if (input[4] <= -0.019499999471008778) {
                                    double tempArray[2] = { 0.927710843373494, 0.07228915662650602 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9801774327696147, 0.019822567230385363 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.42649999260902405) {
                                    double tempArray[2] = { 0.8552927927927928, 0.1447072072072072 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4791666666666667, 0.5208333333333334 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= 0.19250000268220901) {
                                if (input[8] <= 2.090000033378601) {
                                    double tempArray[2] = { 0.4965986394557823, 0.5034013605442177 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6702833031946956, 0.3297166968053044 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.009999990463257) {
                                    double tempArray[2] = { 0.35714285714285715, 0.6428571428571429 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7912488605287147, 0.20875113947128532 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[6] <= 1.6324999928474426) {
                        if (input[8] <= 2.305500030517578) {
                            if (input[0] <= 2.1649999618530273) {
                                if (input[2] <= 0.2654999941587448) {
                                    double tempArray[2] = { 0.9952170713760118, 0.004782928623988227 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9996625421822273, 0.0003374578177727784 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.3135000020265579) {
                                    double tempArray[2] = { 0.9371508379888268, 0.06284916201117319 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7142857142857143, 0.2857142857142857 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= 0.3149999976158142) {
                                if (input[4] <= 0.39249999076128006) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3333333333333333, 0.6666666666666666 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.10349999740719795) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= 0.3034999966621399) {
                            if (input[3] <= -0.03150000050663948) {
                                if (input[0] <= 2.059499979019165) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.34285714285714286, 0.6571428571428571 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 2.9644999504089355) {
                                    double tempArray[2] = { 0.9602745550628683, 0.03972544493713175 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9779129459471119, 0.02208705405288814 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 1.9045000076293945) {
                                if (input[7] <= 2.2820000648498535) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9052044609665427, 0.09479553903345725 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.340999960899353) {
                                    double tempArray[2] = { 0.9408202587393338, 0.059179741260666115 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4696296296296296, 0.5303703703703704 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[1] <= 0.021499999798834324) {
                    if (input[8] <= 1.4164999723434448) {
                        if (input[0] <= 2.445499897003174) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[0] <= 2.4549999237060547) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[8] <= 1.5855000019073486) {
                            if (input[0] <= 2.2440000772476196) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 2.1320000886917114) {
                                if (input[0] <= 2.0554999113082886) {
                                    double tempArray[2] = { 0.9964997307485192, 0.0035002692514808833 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9042338709677419, 0.09576612903225806 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.4200000762939453) {
                                    double tempArray[2] = { 0.6233696623190995, 0.3766303376809005 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4661016949152542, 0.5338983050847458 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= 0.11749999970197678) {
                        if (input[4] <= 0.5355000197887421) {
                            if (input[0] <= 2.1654999256134033) {
                                if (input[6] <= 2.1024999618530273) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9765213224724485, 0.02347867752755151 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 6.551000118255615) {
                                    double tempArray[2] = { 0.7775164628410159, 0.222483537158984 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3038461538461538, 0.6961538461538461 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 2.340499997138977) {
                                if (input[3] <= 0.22949999570846558) {
                                    double tempArray[2] = { 0.497737556561086, 0.502262443438914 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8522727272727273, 0.14772727272727273 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.11550000309944153) {
                                    double tempArray[2] = { 0.25, 0.75 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9455587392550143, 0.054441260744985676 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.5099999904632568) {
                            if (input[3] <= 0.10850000008940697) {
                                if (input[6] <= 1.7695000171661377) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9260846277450455, 0.07391537225495447 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.184000015258789) {
                                    double tempArray[2] = { 0.9974152785755313, 0.002584721424468696 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8870967741935484, 0.11290322580645161 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.1050000190734863) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[4] <= 0.8215000033378601) {
                                    double tempArray[2] = { 0.5459183673469388, 0.45408163265306123 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8787878787878788, 0.12121212121212122 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[0] <= 2.9809999465942383) {
                if (input[8] <= 2.1125000715255737) {
                    if (input[1] <= -0.02650000061839819) {
                        if (input[3] <= 0.21550000458955765) {
                            if (input[6] <= 3.2510000467300415) {
                                if (input[7] <= 2.0640000104904175) {
                                    double tempArray[2] = { 0.8936170212765957, 0.10638297872340426 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.39560931899641577, 0.6043906810035843 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 3.733500123023987) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.08904109589041095, 0.910958904109589 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= 0.6564999967813492) {
                                if (input[5] <= 4.7815001010894775) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.375, 0.625 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.6344999969005585) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 3.137500047683716) {
                            if (input[5] <= 2.353500008583069) {
                                if (input[2] <= 0.029500000178813934) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.79957805907173, 0.20042194092827004 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 1.527999997138977) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.29666254635352285, 0.7033374536464772 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= -0.011500000022351742) {
                                if (input[5] <= 7.832000017166138) {
                                    double tempArray[2] = { 0.08333333333333333, 0.9166666666666666 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.7515000104904175) {
                                    double tempArray[2] = { 0.7402634593356243, 0.25973654066437574 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3475298126064736, 0.6524701873935264 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= 0.0494999997317791) {
                        if (input[8] <= 2.986999988555908) {
                            if (input[0] <= 2.677000045776367) {
                                if (input[4] <= 0.22950000315904617) {
                                    double tempArray[2] = { 0.24985056784219964, 0.7501494321578004 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7724719101123596, 0.22752808988764045 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.18450000137090683) {
                                    double tempArray[2] = { 0.1190176322418136, 0.8809823677581864 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.30142857142857143, 0.6985714285714286 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= 0.3525000065565109) {
                                if (input[1] <= -0.039499999955296516) {
                                    double tempArray[2] = { 0.6608391608391608, 0.33916083916083917 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2538860103626943, 0.7461139896373057 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[2] <= -0.08049999922513962) {
                            if (input[7] <= 2.4285000562667847) {
                                if (input[3] <= 0.10200000554323196) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 3.415000081062317) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.875, 0.125 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.618499994277954) {
                                if (input[7] <= 2.843500018119812) {
                                    double tempArray[2] = { 0.3439922480620155, 0.6560077519379846 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 4.927999973297119) {
                                    double tempArray[2] = { 0.7586520947176685, 0.24134790528233152 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.33182503770739064, 0.6681749622926093 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[0] <= 3.027500033378601) {
                    double tempArray[2] = { 0.0, 1.0 };
                    memcpy(var14, tempArray, 2 * sizeof(double));
                }
                else {
                    if (input[6] <= 2.4614999294281006) {
                        if (input[8] <= 2.4854999780654907) {
                            if (input[7] <= 2.4635000228881836) {
                                if (input[4] <= 0.3110000044107437) {
                                    double tempArray[2] = { 0.2980132450331126, 0.7019867549668874 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7285714285714285, 0.2714285714285714 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.0859999991953373) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.20918367346938777, 0.7908163265306123 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[3] <= 0.24249999970197678) {
                                if (input[8] <= 3.0145000219345093) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[5] <= 5.550500154495239) {
                            if (input[8] <= 2.4079999923706055) {
                                if (input[8] <= 2.197000026702881) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.08545727136431784, 0.9145427286356822 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.077499866485596) {
                                    double tempArray[2] = { 0.05090909090909091, 0.9490909090909091 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.33633633633633636, 0.6636636636636637 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 3.0575000047683716) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[8] <= 3.0145000219345093) {
                                    double tempArray[2] = { 0.004816955684007707, 0.9951830443159922 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        if (input[5] <= 3.4954999685287476) {
            if (input[0] <= 2.534500002861023) {
                if (input[0] <= 2.3174999952316284) {
                    if (input[0] <= 2.0199999809265137) {
                        if (input[1] <= -0.05650000087916851) {
                            if (input[5] <= 3.305500030517578) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[6] <= 4.198499917984009) {
                            if (input[2] <= 0.08049999922513962) {
                                if (input[7] <= 2.39300000667572) {
                                    double tempArray[2] = { 0.9558823529411765, 0.04411764705882353 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.21052631578947367, 0.7894736842105263 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[8] <= 2.0169999599456787) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[8] <= 2.3914999961853027) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[2] <= 0.10050000250339508) {
                        if (input[4] <= 0.28949999809265137) {
                            if (input[6] <= 4.648999929428101) {
                                if (input[5] <= 3.2235000133514404) {
                                    double tempArray[2] = { 0.8557692307692307, 0.14423076923076922 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.2719999998807907) {
                                    double tempArray[2] = { 0.018867924528301886, 0.9811320754716981 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2, 0.8 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[5] <= 3.2494999170303345) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
            else {
                if (input[7] <= 2.1024999618530273) {
                    if (input[5] <= 3.34499990940094) {
                        double tempArray[2] = { 0.0, 1.0 };
                        memcpy(var14, tempArray, 2 * sizeof(double));
                    }
                    else {
                        if (input[7] <= 1.9039999842643738) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[6] <= 4.697999954223633) {
                                if (input[6] <= 4.696500062942505) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
                else {
                    if (input[7] <= 2.118000030517578) {
                        if (input[2] <= 0.016499999910593033) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[4] <= 0.04050000011920929) {
                            if (input[2] <= 0.042500000447034836) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[4] <= 0.023500001057982445) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var14, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
        }
        else {
            if (input[1] <= 0.005499999970197678) {
                if (input[0] <= 2.3734999895095825) {
                    if (input[5] <= 5.243499994277954) {
                        if (input[2] <= 0.06149999983608723) {
                            if (input[7] <= 2.118000030517578) {
                                if (input[8] <= 2.01800000667572) {
                                    double tempArray[2] = { 0.7898089171974523, 0.21019108280254778 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.26180257510729615, 0.7381974248927039 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.06549999862909317) {
                                    double tempArray[2] = { 0.4797047970479705, 0.5202952029520295 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7644557823129252, 0.23554421768707484 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.0195000171661377) {
                                if (input[4] <= 0.08249999955296516) {
                                    double tempArray[2] = { 0.9891304347826086, 0.010869565217391304 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8522727272727273, 0.14772727272727273 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.07600000128149986) {
                                    double tempArray[2] = { 0.2972972972972973, 0.7027027027027027 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7872340425531915, 0.2127659574468085 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.177000045776367) {
                            if (input[0] <= 2.059999942779541) {
                                if (input[6] <= 3.847499966621399) {
                                    double tempArray[2] = { 0.9926315789473684, 0.007368421052631579 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9996222851746931, 0.00037771482530689327 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.256500005722046) {
                                    double tempArray[2] = { 0.6148648648648649, 0.38513513513513514 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9575569358178054, 0.042443064182194616 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.4264999628067017) {
                                if (input[8] <= 2.6074999570846558) {
                                    double tempArray[2] = { 0.9101796407185628, 0.08982035928143713 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2, 0.8 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= -0.04450000077486038) {
                                    double tempArray[2] = { 0.22972972972972974, 0.7702702702702703 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6669484361792054, 0.3330515638207946 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.7730000019073486) {
                        if (input[8] <= 2.531999945640564) {
                            if (input[4] <= 0.3125) {
                                if (input[6] <= 6.5299999713897705) {
                                    double tempArray[2] = { 0.3252405459834415, 0.6747594540165586 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.705685618729097, 0.29431438127090304 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.445499897003174) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6403269754768393, 0.35967302452316074 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 2.972000002861023) {
                                if (input[4] <= 0.3009999990463257) {
                                    double tempArray[2] = { 0.10076891946580332, 0.8992310805341966 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7559523809523809, 0.24404761904761904 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 8.54699993133545) {
                                    double tempArray[2] = { 0.8958333333333334, 0.10416666666666667 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.3464999943971634) {
                            if (input[8] <= 2.9494999647140503) {
                                if (input[3] <= -0.10450000315904617) {
                                    double tempArray[2] = { 0.9565217391304348, 0.043478260869565216 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.008499810513778355, 0.9915001894862217 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.16899999976158142) {
                                    double tempArray[2] = { 0.9629629629629629, 0.037037037037037035 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.07264957264957266, 0.9273504273504274 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.971500039100647) {
                                if (input[4] <= 0.460999995470047) {
                                    double tempArray[2] = { 0.2743362831858407, 0.7256637168141593 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9574468085106383, 0.0425531914893617 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.10050000250339508) {
                                    double tempArray[2] = { 0.05263157894736842, 0.9473684210526315 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.005964214711729622, 0.9940357852882704 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[6] <= 4.696500062942505) {
                    if (input[4] <= 0.20750000327825546) {
                        if (input[0] <= 2.384500026702881) {
                            if (input[7] <= 2.9394999742507935) {
                                if (input[8] <= 2.037999987602234) {
                                    double tempArray[2] = { 0.9762168141592921, 0.023783185840707963 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8667152221412965, 0.13328477785870357 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.1514999866485596) {
                                    double tempArray[2] = { 0.2, 0.8 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.6080000400543213) {
                                if (input[1] <= 0.08950000256299973) {
                                    double tempArray[2] = { 0.29138576779026215, 0.7086142322097378 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5503649635036496, 0.44963503649635034 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 6.628999948501587) {
                                    double tempArray[2] = { 0.0660377358490566, 0.9339622641509434 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2614213197969543, 0.7385786802030457 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[2] <= 0.34549999237060547) {
                            if (input[1] <= 0.03150000050663948) {
                                if (input[0] <= 2.659500002861023) {
                                    double tempArray[2] = { 0.9174311926605505, 0.08256880733944955 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.734499931335449) {
                                    double tempArray[2] = { 0.9536213468869124, 0.04637865311308768 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.022847100175746926, 0.9771528998242531 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.7014999389648438) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var14, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 4.317999839782715) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.625, 0.375 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[4] <= 0.016500000841915607) {
                        if (input[8] <= 2.3625000715255737) {
                            if (input[5] <= 4.291499853134155) {
                                if (input[0] <= 2.6780000925064087) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= -0.1314999982714653) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.46321070234113715, 0.5367892976588629 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 8.492499828338623) {
                                if (input[3] <= 0.04649999924004078) {
                                    double tempArray[2] = { 0.3938879456706282, 0.6061120543293718 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7175572519083969, 0.2824427480916031 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.40749990940094) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[4] <= 0.24949999898672104) {
                            if (input[0] <= 2.3524999618530273) {
                                if (input[0] <= 2.0959999561309814) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6735395189003437, 0.32646048109965636 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.271000012755394) {
                                    double tempArray[2] = { 0.08704780361757106, 0.912952196382429 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7741935483870968, 0.22580645161290322 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 6.042999982833862) {
                                if (input[8] <= 2.521000027656555) {
                                    double tempArray[2] = { 0.4035493827160494, 0.5964506172839507 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6735218508997429, 0.3264781491002571 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.07500000298023224) {
                                    double tempArray[2] = { 0.8362573099415205, 0.16374269005847952 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.15151515151515152, 0.8484848484848485 };
                                    memcpy(var14, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    add_vectors(var3, var14, 2, var2);
    double var15[2];
    if (input[8] <= 2.188499927520752) {
        if (input[0] <= 2.3615000247955322) {
            if (input[0] <= 2.1200000047683716) {
                if (input[7] <= 2.0089999437332153) {
                    if (input[0] <= 2.0175000429153442) {
                        if (input[8] <= 2.1100000143051147) {
                            if (input[1] <= 0.004500000039115548) {
                                if (input[6] <= 4.319999933242798) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7333333333333333, 0.26666666666666666 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[1] <= -0.017500000074505806) {
                                if (input[7] <= 1.8779999613761902) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8152173913043478, 0.18478260869565216 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 3.5485000610351562) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.970873786407767, 0.02912621359223301 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 2.99399995803833) {
                            if (input[0] <= 2.0619999170303345) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[8] <= 1.8195000290870667) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9479674796747968, 0.05203252032520325 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.023000000044703484) {
                                if (input[8] <= 1.690500020980835) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8949579831932774, 0.10504201680672269 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 1.7335000038146973) {
                                    double tempArray[2] = { 0.9545454545454546, 0.045454545454545456 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.35018050541516244, 0.6498194945848376 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[1] <= -0.09750000014901161) {
                        if (input[1] <= -0.1054999977350235) {
                            if (input[8] <= 2.040000081062317) {
                                if (input[1] <= -0.1145000010728836) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9960629921259843, 0.003937007874015748 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.0479999780654907) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 3.9390000104904175) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[6] <= 4.114500045776367) {
                                    double tempArray[2] = { 0.13043478260869565, 0.8695652173913043 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9565217391304348, 0.043478260869565216 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[1] <= 0.06949999928474426) {
                            if (input[4] <= 0.49300000071525574) {
                                if (input[7] <= 2.9665000438690186) {
                                    double tempArray[2] = { 0.9933383651100974, 0.006661634889902537 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7701149425287356, 0.22988505747126436 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.0509999990463257) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.38095238095238093, 0.6190476190476191 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.958500027656555) {
                                if (input[4] <= 0.5360000133514404) {
                                    double tempArray[2] = { 0.9997794550088218, 0.00022054499117820034 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9906832298136646, 0.009316770186335404 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.023999999277293682) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[1] <= 0.017500000074505806) {
                    if (input[5] <= 5.481499910354614) {
                        if (input[6] <= 3.593999981880188) {
                            if (input[8] <= 1.7889999747276306) {
                                if (input[3] <= 0.12049999833106995) {
                                    double tempArray[2] = { 0.7964470762398224, 0.20355292376017764 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9841628959276018, 0.01583710407239819 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.1804999113082886) {
                                    double tempArray[2] = { 0.8549382716049383, 0.14506172839506173 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5332339791356184, 0.4667660208643815 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.2640000581741333) {
                                if (input[2] <= -0.08150000125169754) {
                                    double tempArray[2] = { 0.75, 0.25 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.05921052631578947, 0.9407894736842105 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= -0.02850000001490116) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.34210526315789475, 0.6578947368421053 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[1] <= -0.05550000071525574) {
                            if (input[0] <= 2.340000033378601) {
                                if (input[7] <= 2.6480000019073486) {
                                    double tempArray[2] = { 0.978675645342312, 0.02132435465768799 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7810650887573964, 0.21893491124260356 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[8] <= 1.5354999899864197) {
                                if (input[0] <= 2.305999994277954) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.16129032258064516, 0.8387096774193549 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= -0.035499999299645424) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9060518731988473, 0.09394812680115273 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[2] <= 0.2344999983906746) {
                        if (input[1] <= 0.11950000002980232) {
                            if (input[4] <= 0.5744999945163727) {
                                if (input[7] <= 1.816499948501587) {
                                    double tempArray[2] = { 0.2, 0.8 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8405331620362936, 0.15946683796370645 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.2850000858306885) {
                                    double tempArray[2] = { 0.045454545454545456, 0.9545454545454546 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 2.6535000801086426) {
                                if (input[6] <= 2.2924998998641968) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.516320474777448, 0.4836795252225519 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.6809999942779541) {
                                    double tempArray[2] = { 0.9413317688471693, 0.05866823115283074 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.42105263157894735, 0.5789473684210527 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.340999960899353) {
                            if (input[8] <= 2.0364999771118164) {
                                if (input[0] <= 2.240499973297119) {
                                    double tempArray[2] = { 0.8235574630424416, 0.1764425369575584 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9805309734513274, 0.019469026548672566 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 2.9140000343322754) {
                                    double tempArray[2] = { 0.3668639053254438, 0.6331360946745562 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.697508896797153, 0.302491103202847 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.3480000495910645) {
                                if (input[8] <= 2.083500027656555) {
                                    double tempArray[2] = { 0.48655256723716384, 0.5134474327628362 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[1] <= 0.04350000061094761) {
                if (input[0] <= 2.759999990463257) {
                    if (input[2] <= 0.5940000116825104) {
                        if (input[6] <= 3.475499987602234) {
                            if (input[8] <= 2.111999988555908) {
                                if (input[5] <= 2.9299999475479126) {
                                    double tempArray[2] = { 0.1590909090909091, 0.8409090909090909 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6203389830508474, 0.37966101694915255 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= -0.0005000000237487257) {
                                    double tempArray[2] = { 0.6705882352941176, 0.32941176470588235 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.19508196721311474, 0.8049180327868852 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= -0.12049999833106995) {
                                if (input[0] <= 2.416000008583069) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.3419354838709677, 0.6580645161290323 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.1269999742507935) {
                                    double tempArray[2] = { 0.40034529061960483, 0.5996547093803951 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2079107505070994, 0.7920892494929006 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 1.8845000267028809) {
                            if (input[5] <= 7.105000019073486) {
                                if (input[1] <= -0.005499999970197678) {
                                    double tempArray[2] = { 0.8904109589041096, 0.1095890410958904 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4489795918367347, 0.5510204081632653 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.7794999778270721) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8846153846153846, 0.11538461538461539 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.135499954223633) {
                                if (input[7] <= 2.0329999923706055) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9534883720930233, 0.046511627906976744 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[2] <= 0.7069999873638153) {
                                    double tempArray[2] = { 0.696969696969697, 0.30303030303030304 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[4] <= 0.5819999873638153) {
                        if (input[6] <= 0.6330000162124634) {
                            if (input[2] <= 4.60450005531311) {
                                if (input[1] <= -0.07499999925494194) {
                                    double tempArray[2] = { 0.6378378378378379, 0.3621621621621622 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.05, 0.95 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[2] <= 0.5144999921321869) {
                                if (input[0] <= 2.7730000019073486) {
                                    double tempArray[2] = { 0.1715686274509804, 0.8284313725490197 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.006479933110367893, 0.9935200668896321 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.5175000429153442) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 3.140500068664551) {
                            if (input[6] <= 1.55349999666214) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var15, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.9524999856948853) {
                    if (input[4] <= 0.29649999737739563) {
                        if (input[0] <= 2.6880000829696655) {
                            if (input[2] <= -0.09349999949336052) {
                                if (input[3] <= 0.02949999924749136) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8841698841698842, 0.11583011583011583 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 2.3709999322891235) {
                                    double tempArray[2] = { 0.7738095238095238, 0.2261904761904762 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5562735595045772, 0.4437264404954227 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[1] <= 0.2305000051856041) {
                                if (input[6] <= 2.063499927520752) {
                                    double tempArray[2] = { 0.6967213114754098, 0.30327868852459017 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.1340057636887608, 0.8659942363112392 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[4] <= 0.07049999758601189) {
                                    double tempArray[2] = { 0.8333333333333334, 0.16666666666666666 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[6] <= 1.6765000224113464) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var15, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[7] <= 2.24399995803833) {
                                if (input[1] <= 0.28450000286102295) {
                                    double tempArray[2] = { 0.5621761658031088, 0.4378238341968912 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9649122807017544, 0.03508771929824561 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.2995000034570694) {
                                    double tempArray[2] = { 0.8323353293413174, 0.16766467065868262 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.40540540540540543, 0.5945945945945946 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[4] <= 0.5180000066757202) {
                        if (input[2] <= 0.6770000159740448) {
                            if (input[1] <= 0.3955000042915344) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[1] <= 0.4294999986886978) {
                                    double tempArray[2] = { 0.5454545454545454, 0.45454545454545453 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.3289999961853027) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[7] <= 2.6605000495910645) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2857142857142857, 0.7142857142857143 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[5] <= 3.2975000143051147) {
                            if (input[3] <= 0.43150000274181366) {
                                if (input[2] <= 0.27650000154972076) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9166666666666666, 0.08333333333333333 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[0] <= 3.059000015258789) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[1] <= 0.22949999570846558) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        if (input[0] <= 2.4609999656677246) {
            if (input[2] <= 0.05250000022351742) {
                if (input[4] <= 0.26350000500679016) {
                    if (input[5] <= 6.963500022888184) {
                        if (input[0] <= 2.1320000886917114) {
                            if (input[2] <= -0.08449999988079071) {
                                if (input[7] <= 2.4045000076293945) {
                                    double tempArray[2] = { 0.7666666666666667, 0.23333333333333334 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9455882352941176, 0.054411764705882354 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.0464999675750732) {
                                    double tempArray[2] = { 0.9936708860759493, 0.006329113924050633 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9318181818181818, 0.06818181818181818 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 4.093999862670898) {
                                if (input[0] <= 2.18149995803833) {
                                    double tempArray[2] = { 0.31494252873563217, 0.6850574712643678 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6827195467422096, 0.31728045325779036 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.371500015258789) {
                                    double tempArray[2] = { 0.6603375527426161, 0.339662447257384 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.34384581690757776, 0.6561541830924222 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.177000045776367) {
                            if (input[0] <= 2.0679999589920044) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[8] <= 2.6735000610351562) {
                                    double tempArray[2] = { 0.9511316872427984, 0.048868312757201646 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5797101449275363, 0.42028985507246375 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[7] <= 2.472999930381775) {
                                if (input[2] <= -0.014500000048428774) {
                                    double tempArray[2] = { 0.877221324717286, 0.12277867528271405 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6428571428571429, 0.35714285714285715 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.7015000581741333) {
                                    double tempArray[2] = { 0.5181040157998683, 0.4818959842001317 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7813712807244502, 0.2186287192755498 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[8] <= 2.5665000677108765) {
                        if (input[6] <= 6.226999998092651) {
                            if (input[8] <= 2.3799999952316284) {
                                if (input[3] <= 0.07849999889731407) {
                                    double tempArray[2] = { 0.5578947368421052, 0.4421052631578947 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8717434869739479, 0.1282565130260521 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.1185000017285347) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9886363636363636, 0.011363636363636364 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var15, tempArray, 2 * sizeof(double));
                        }
                    }
                    else {
                        if (input[7] <= 2.1109999418258667) {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var15, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[0] <= 2.1675000190734863) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[4] <= 1.0359999537467957) {
                                    double tempArray[2] = { 0.6644067796610169, 0.33559322033898303 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (input[7] <= 2.3070000410079956) {
                    if (input[2] <= 0.3084999918937683) {
                        if (input[0] <= 2.166499972343445) {
                            if (input[4] <= 0.3645000010728836) {
                                if (input[7] <= 2.240499973297119) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9897959183673469, 0.01020408163265306 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 3.6955000162124634) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= 0.42149999737739563) {
                                if (input[7] <= 1.9445000290870667) {
                                    double tempArray[2] = { 0.1774193548387097, 0.8225806451612904 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9212598425196851, 0.07874015748031496 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[5] <= 3.3865000009536743) {
                            if (input[1] <= 0.22950000315904617) {
                                if (input[6] <= 1.7290000319480896) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9962825278810409, 0.0037174721189591076 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.2109999656677246) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.083500027656555) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[5] <= 3.4234999418258667) {
                                    double tempArray[2] = { 0.14285714285714285, 0.8571428571428571 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[4] <= 0.24450000375509262) {
                        if (input[1] <= -0.03150000050663948) {
                            if (input[8] <= 2.590499997138977) {
                                if (input[7] <= 2.6729999780654907) {
                                    double tempArray[2] = { 0.45229007633587787, 0.5477099236641222 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7575757575757576, 0.24242424242424243 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= -0.032500000670552254) {
                                    double tempArray[2] = { 0.9620689655172414, 0.03793103448275862 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.375, 0.625 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 2.0649999380111694) {
                                if (input[7] <= 2.315500020980835) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9986553115194979, 0.001344688480502017 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.853500008583069) {
                                    double tempArray[2] = { 0.7513846153846154, 0.24861538461538463 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.35978835978835977, 0.6402116402116402 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[1] <= 0.10649999976158142) {
                            if (input[5] <= 5.169500112533569) {
                                if (input[3] <= 0.37300001084804535) {
                                    double tempArray[2] = { 0.774468085106383, 0.225531914893617 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.3554999828338623) {
                                    double tempArray[2] = { 0.9970887918486172, 0.002911208151382824 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.9164345403899722, 0.08356545961002786 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[2] <= 0.48399999737739563) {
                                if (input[6] <= 2.409999966621399) {
                                    double tempArray[2] = { 0.9642184557438794, 0.035781544256120526 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[7] <= 2.603999972343445) {
                                    double tempArray[2] = { 0.2857142857142857, 0.7142857142857143 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8666666666666667, 0.13333333333333333 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (input[6] <= 3.1295000314712524) {
                if (input[3] <= 0.33650000393390656) {
                    if (input[0] <= 2.9809999465942383) {
                        if (input[7] <= 2.634999990463257) {
                            if (input[4] <= 0.1834999993443489) {
                                if (input[5] <= 3.1880000829696655) {
                                    double tempArray[2] = { 0.07166123778501629, 0.9283387622149837 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4522058823529412, 0.5477941176470589 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.878499984741211) {
                                    double tempArray[2] = { 0.6710280373831776, 0.32897196261682243 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4186046511627907, 0.5813953488372093 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[5] <= 7.573500156402588) {
                                if (input[2] <= 0.02049999963492155) {
                                    double tempArray[2] = { 0.004784688995215311, 0.9952153110047847 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2639593908629442, 0.7360406091370558 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 8.079500198364258) {
                                    double tempArray[2] = { 0.630057803468208, 0.3699421965317919 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0273972602739726, 0.9726027397260274 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 3.1095000505447388) {
                            if (input[0] <= 3.027500033378601) {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                if (input[5] <= 6.499499797821045) {
                                    double tempArray[2] = { 0.7916666666666666, 0.20833333333333334 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[0] <= 3.2975000143051147) {
                                if (input[3] <= 0.19449999928474426) {
                                    double tempArray[2] = { 0.20707070707070707, 0.7929292929292929 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.7463414634146341, 0.25365853658536586 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[1] <= 0.29350000619888306) {
                                    double tempArray[2] = { 0.010514018691588784, 0.9894859813084113 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.6923076923076923, 0.3076923076923077 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[0] <= 2.7445000410079956) {
                        if (input[5] <= 2.5015000104904175) {
                            double tempArray[2] = { 1.0, 0.0 };
                            memcpy(var15, tempArray, 2 * sizeof(double));
                        }
                        else {
                            if (input[0] <= 2.7345000505447388) {
                                if (input[1] <= 0.04350000061094761) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.46017699115044247, 0.5398230088495575 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                    }
                    else {
                        if (input[7] <= 2.4320000410079956) {
                            if (input[0] <= 2.927500009536743) {
                                double tempArray[2] = { 1.0, 0.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            double tempArray[2] = { 0.0, 1.0 };
                            memcpy(var15, tempArray, 2 * sizeof(double));
                        }
                    }
                }
            }
            else {
                if (input[0] <= 2.853999972343445) {
                    if (input[4] <= 0.31450000405311584) {
                        if (input[2] <= -0.05949999950826168) {
                            if (input[1] <= -0.021499999798834324) {
                                if (input[1] <= -0.07349999994039536) {
                                    double tempArray[2] = { 0.3573243014394581, 0.642675698560542 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.15739644970414202, 0.842603550295858 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[6] <= 6.111999988555908) {
                                    double tempArray[2] = { 0.5640535372848948, 0.4359464627151052 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[4] <= -0.09149999916553497) {
                                if (input[2] <= 0.008999999612569809) {
                                    double tempArray[2] = { 0.12195121951219512, 0.8780487804878049 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 3.971500039100647) {
                                    double tempArray[2] = { 0.5423728813559322, 0.4576271186440678 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.15645412130637637, 0.8435458786936236 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[8] <= 2.322499990463257) {
                            if (input[0] <= 2.678499937057495) {
                                if (input[2] <= 0.04599999915808439) {
                                    double tempArray[2] = { 0.06896551724137931, 0.9310344827586207 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                double tempArray[2] = { 0.0, 1.0 };
                                memcpy(var15, tempArray, 2 * sizeof(double));
                            }
                        }
                        else {
                            if (input[3] <= 0.3034999966621399) {
                                if (input[7] <= 2.7795000076293945) {
                                    double tempArray[2] = { 0.8511326860841424, 0.1488673139158576 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5632911392405063, 0.43670886075949367 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[0] <= 2.744499921798706) {
                                    double tempArray[2] = { 0.9728506787330317, 0.027149321266968326 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.4166666666666667, 0.5833333333333334 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
                else {
                    if (input[4] <= 0.18949999660253525) {
                        if (input[6] <= 3.3604999780654907) {
                            if (input[1] <= 0.15650000423192978) {
                                if (input[4] <= 0.14150000363588333) {
                                    double tempArray[2] = { 0.001336898395721925, 0.9986631016042781 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.06, 0.94 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[8] <= 2.466499924659729) {
                                    double tempArray[2] = { 0.9666666666666667, 0.03333333333333333 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[8] <= 3.0429999828338623) {
                                if (input[1] <= 0.17150000482797623) {
                                    double tempArray[2] = { 0.002543558438255119, 0.9974564415617448 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.060669456066945605, 0.9393305439330544 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.15399999916553497) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.8333333333333334, 0.16666666666666666 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                    else {
                        if (input[0] <= 2.971500039100647) {
                            if (input[4] <= 0.2995000034570694) {
                                if (input[0] <= 2.905500054359436) {
                                    double tempArray[2] = { 0.13175675675675674, 0.8682432432432432 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[5] <= 5.141000032424927) {
                                    double tempArray[2] = { 1.0, 0.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.2585858585858586, 0.7414141414141414 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                        else {
                            if (input[6] <= 6.658999919891357) {
                                if (input[0] <= 3.027500033378601) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.03184713375796178, 0.9681528662420382 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                            else {
                                if (input[3] <= 0.15400000289082527) {
                                    double tempArray[2] = { 0.0, 1.0 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                                else {
                                    double tempArray[2] = { 0.5333333333333333, 0.4666666666666667 };
                                    memcpy(var15, tempArray, 2 * sizeof(double));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    add_vectors(var2, var15, 2, var1);
    mul_vector_number(var1, 0.125, 2, var0);
    memcpy(output, var0, 2 * sizeof(double));
}
