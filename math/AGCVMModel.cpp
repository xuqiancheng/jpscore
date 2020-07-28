/**
* \file       AGCVM.cpp
* \date        Jun 26, 2019
* \version     v0.8
* \copyright   <2009-2015> Forschungszentrum Jülich GmbH. All rights reserved.
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
* Anticipation generalized collision-free velocity model: Qiancheng (8)
*
*
**/

# define NOMINMAX
#include "../pedestrian/Pedestrian.h"
#include "../mpi/LCGrid.h"
#include "../geometry/Wall.h"
#include "../geometry/SubRoom.h"

#include "AGCVMModel.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads()  1
#endif

using std::vector;
using std::string;

int ShowInfo = 0;
int TestID1 = 72;
int TestID2 = -1;

int IfDrill = 0;
int IfCorrection = 0;
int VerticalNoInf = 0;

AGCVMModel::AGCVMModel(std::shared_ptr<DirectionStrategy> dir,
    double aped, double Dped, double awall, double Dwall,
    double Ts, double Td, int GCVM,
    double lb, double rb, double ub, double db, double co,
    int Anticipation, int Cooperation, int AttracForce, int Push,
    double AntiT, double CoopT, double CoreSize)
{
    _direction = dir;

    // Force_rep_PED Parameter
    _aPed = aped;
    _DPed = Dped;

    // Force_rep_WALL Parameter
    _aWall = awall;
    _DWall = Dwall;

    // GCVM Parameter
    _Ts = Ts; // Speed module
    _Td = Td; // Direction module
    _GCVM = GCVM; // If using GCVM

    // Boundary Case
    _LeftBoundary = lb;
    _RightBoundary = rb;
    _UpBoundary = ub;
    _DownBoundary = db;
    _CutOff = co;

    _Anticipation = Anticipation;
    _Cooperation = Cooperation;
    _AttractiveForce = AttracForce;
    _PushingForce = Push;
    _AntiTime = AntiT;
    _CoopTime = CoopT;
    _CoreSize = CoreSize;
}


AGCVMModel::~AGCVMModel()
{

}

string AGCVMModel::GetDescription()
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
bool AGCVMModel::Init(Building* building)
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
                "ERROR:\tGCVMModel::Init() can not initialise route. ped %d is deleted in Room %d %d.\n", ped->GetID(), ped->GetRoomID(), ped->GetSubRoomID());
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

        if (GetPushing())
        {
            ped->SetCore(GetCoreSize());
        }
        else
        {
            ped->SetCore(ped->GetEllipse().GetBmax());
        }
    }
    return true;
}

// Compute the status of pedestrians at next timestep
void AGCVMModel::ComputeNextTimeStep(double current, double deltaT, Building* building, int periodic)
{
    const vector< Pedestrian* >& allPeds = building->GetAllPedestrians();
    vector<Pedestrian*> pedsToRemove;
    pedsToRemove.reserve(500);
    unsigned long nSize;
    nSize = allPeds.size();
    //---------------------------------------------------------------
    vector< Point > normal_acc = vector<Point >();
    normal_acc.reserve(nSize);

    vector< Point > result_dir = vector<Point >();
    result_dir.reserve(nSize);

    int start = 0;
    int end = nSize - 1;

    // Calculate direction for each pedestrians
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped1 = allPeds[p];
        Point p1 = ped1->GetPos();
        Room* room = building->GetRoom(ped1->GetRoomID());
        SubRoom* subroom = room->GetSubRoom(ped1->GetSubRoomID());
        // Effect from neighbours
        Point IniDirection = DesireDirection(ped1, room);//desired moving direction, direction3, and using core_size here.
        vector<Pedestrian*> neighbours;
        building->GetGrid()->GetNeighbourhood(ped1, neighbours);
        int size = (int)neighbours.size();
        Point repPedTurn = Point(0, 0);
        Point repPedPush = Point(0, 0);
        // Calculate the influence from neighbours
        for (int i = 0; i < size; i++)
        {
            Pedestrian* ped2 = neighbours[i];
            Point p2 = ped2->GetPos();
            //Check if two pedestrians can see each other
            SubRoom* subroom2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());
            vector<SubRoom*> emptyVector;
            emptyVector.push_back(subroom);
            emptyVector.push_back(subroom2);
            bool isVisible = building->IsVisible(p1, p2, emptyVector, false);
            if (!isVisible)
            {
                //continue;// we can delete it to make different simulations
            }
            // whole building, not only the pedestrians in same subroom or next subroom
            if (GetPushing() == 1)
            {
                repPedPush += ForceRepPedPush(ped1, ped2, building, periodic);//new method
            }
            if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom->IsDirectlyConnectedWith(subroom2))
            {
                Point Force = ForceRepPed(ped1, ped2, building, periodic);//new method
                repPedTurn += Force;
            }
        } //for i

        Point repWall = ForceRepRoom(ped1, subroom);
        /*
        // the influence from next subroom should be considered
        if (ped1->GetExitLine()->DistTo(p1) < ped1->GetEllipse().GetBmax())
        {
            for (const auto & subr : subroom->GetNeighbors())
            {
                repWall = repWall + ForceRepRoom(ped1, subr);
            }
        }
        */

        Point direction;
        Point a_direction = ped1->GetMoveDirection();
        Point d_direction = IniDirection + repPedTurn + repWall + repPedPush;
        // calculate new direction
        if (ped1->GetTryCoop() == 1)
        {
            //d_direction = a_direction;
        }
        if (IfDrill == 1)
        {
            bool drill = Drill(ped1, neighbours, building, subroom, IniDirection, periodic);
            if (drill == true)
            {
                d_direction = IniDirection.Normalized();
            }
        }
        if (IfCorrection == 1)
        {
            d_direction = CorrectD(ped1, d_direction, subroom);
        }

        Point AccTu = Point(0, 0);
        double angle_tau = GetTd();
        d_direction = d_direction.Normalized();
        AccTu = (d_direction.Normalized()*ped1->GetV0Norm() - ped1->GetV()) / angle_tau;
        direction = ped1->GetV() + AccTu * deltaT;
        direction = direction.Normalized();
        //direction = IniDirection.Normalized();
        //direction = CorrectD(ped1, direction, subroom);

        if (GetGCVMU() == 0)
        {
            direction = d_direction.Normalized();//original method
        }

        // Test code
        if (ShowInfo && (ped1->GetID() == TestID1 || ped1->GetID() == TestID2))
        {
            if (current > 0 && current < 5)
            {
                printf("TEST: \ttime(%f), ID(%d), Pos(%f,%f), e0(%f,%f), ei(%f,%f)\n",
                    current, ped1->GetID(), ped1->GetPos()._x, ped1->GetPos()._y, IniDirection._x, IniDirection._y, direction._x, direction._y);
                printf("\trepPedTurn(%f,%f), repPedPush(%f,%f), repWall(%f,%f)\n\n",
                    repPedTurn._x, repPedTurn._y, repPedPush._x, repPedPush._y, repWall._x, repWall._y);
            }
        }
        result_dir.push_back(direction);

    }

    // Update direction of each pedestrian
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped = allPeds[p];
        Point direction = result_dir[p];
        ped->SetMoveDirection(direction);
    }

    // Calculate speed
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped1 = allPeds[p];
        Room* room = building->GetRoom(ped1->GetRoomID());
        SubRoom* subroom = room->GetSubRoom(ped1->GetSubRoomID());

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
            if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom->IsDirectlyConnectedWith(subroom2) || 1)
            {
                spacings.push_back(GetSpacing(ped1, ped2, periodic));
            }
        }//for i

        // Calculate min spacing
        std::sort(spacings.begin(), spacings.end(), sort_pred_agcvm());
        double spacing = spacings.size() == 0 ? 100 : spacings[0].first;
        double spacing_wall = GetSpacingRoom(ped1, subroom);
        // neighbour subroom needs to be considered
        for (const auto & subr : subroom->GetNeighbors())
        {
            double swn = GetSpacingRoom(ped1, subr);
            spacing_wall = spacing_wall > swn ? swn : spacing_wall;
        }
        spacing = spacing < spacing_wall ? spacing : spacing_wall;

        // Optimal speed function
        Point speed;
        Point ei = ped1->GetMoveDirection();
        speed = ei * OptimalSpeed(ped1, spacing);
        normal_acc.push_back(speed);
    } // for p

    // Update speed of each pedestrian
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped = allPeds[p];
        Point speed = normal_acc[p];
        ped->SetV(speed);;
    }

    // Cooperation
    if (GetCooperation() == 1)
    {
        for (int p = start; p <= end; ++p)
        {
            Pedestrian* ped1 = allPeds[p];
            ped1->SetTryCoop(0);
            vector<Pedestrian*> neighbours;
            building->GetGrid()->GetNeighbourhood(ped1, neighbours);
            int size = (int)neighbours.size();
            vector<my_pair> ttcs = vector<my_pair>();
            ttcs.reserve(size);
            for (int i = 0; i < size; i++)
            {
                Pedestrian* ped2 = neighbours[i];
                my_pair ttc = JudgeCollision(ped1, ped2, building, periodic);
                ttc.second = i; //I don't want to use ID here
                ttcs.push_back(ttc);
            }
            std::sort(ttcs.begin(), ttcs.end(), sort_pred_agcvm());
            double mttc = ttcs.size() == 0 ? FLT_MAX : ttcs[0].first;
            if (mttc < GetCoopT())
            {
                Pedestrian* ped_mttc = neighbours[ttcs[0].second];
                ped1->SetMTTCP(ttcs[0].second);
                double coop1 = ped1->GetCooperation();
                double coop2 = ped_mttc->GetCooperation();
                double s1 = GetSpacing(ped1, ped_mttc, periodic).first;
                double s2 = GetSpacing(ped_mttc, ped1, periodic).first;
                /*
                if (ped1->GetV().Norm() < J_EPS&&ped_mttc->GetV().Norm() > J_EPS)
                {
                    ped1->SetTryCoop(1);
                }
                */
                if ((s1 < 100 && s2>100) || (coop1 > coop2&&s1 > 100 && s2 > 100) || (coop1 > coop2&&s1 < 100 && s2 < 100))
                {
                    ped1->SetTryCoop(1);
                }
            }
        }
    }

    //Update everything
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped = allPeds[p];
        Point v_neu = normal_acc[p];
        if (ped->GetTryCoop() == 1)
        {
            v_neu = Point(0, 0);
        }
        Point dir_neu = result_dir[p];
        UpdatePed(ped, v_neu, dir_neu, deltaT, periodic);
    }

    // remove the pedestrians that have left the building
    for (unsigned int p = 0; p < pedsToRemove.size(); p++)
    {
        building->DeletePedestrian(pedsToRemove[p]);
    }
    pedsToRemove.clear();
}

/*----------Functions important----------*/
Point AGCVMModel::DesireDirection(Pedestrian* ped, Room* room) const
{
    const Point target = this->GetDirection()->GetTarget(room, ped);
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

Point AGCVMModel::ForceRepPed(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const
{
    Point F_rep(0.0, 0.0);
    Point p1 = ped1->GetPos();
    Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();
    Point distp12 = p2 - p1;
    double Distance = distp12.Norm();
    Point ep12;
    if (Distance >= J_EPS)
    {
        ep12 = distp12.Normalized();
    }
    else
    {
        printf("ERROR: \tIn AGCVMModel::forcePedPed() ep12 can not be calculated!!!\n");
        printf("ped1 %d  ped2 %d\n", ped1->GetID(), ped2->GetID());
        printf("ped1 at (%f, %f), ped2 at (%f, %f)\n", p1._x, p1._y, p2._x, p2._y);
        exit(EXIT_FAILURE);
    }
    // Distance between boundaries
    double dist = Distance - (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());

    if (GetGCVMU() == 0)// Using CVM
    {
        double R_ij = _aPed * exp((-dist) / _DPed);
        F_rep = ep12 * (-R_ij);
        return F_rep;
    }

    //Vision area
    Point e1 = ped1->GetMoveDirection();
    Room* room1 = building->GetRoom(ped1->GetRoomID());
    Point d1 = DesireDirection(ped1, room1);

    Point e2 = ped2->GetMoveDirection();
    Room* room2 = building->GetRoom(ped2->GetRoomID());
    Point d2 = DesireDirection(ped2, room2);

    double condition1 = d1.ScalarProduct(ep12);
    double condition2 = e1.ScalarProduct(ep12);

    //Anticipation
    double S_Gap = (e1.ScalarProduct(ep12)*ped1->GetV0Norm() - e2.ScalarProduct(ep12)*ped2->GetV0Norm());
    //S_Gap = (ei.ScalarProduct(ep12) - ei2.ScalarProduct(ep12));
    //S_Gap = (ei.ScalarProduct(ep12)*ped1->GetV().Norm() - ei2.ScalarProduct(ep12)*ped2->GetV().Norm());
    double multi_d = d1.ScalarProduct(d2);
    double beta = (1 - multi_d);
    double Dis_Gap = GetAnticipation() ? S_Gap * GetAntiT()* beta : 0;// Anticipation

    if (condition1 >= 0 || condition2 >= 0)
    {
        Point infd = GetInfDirection(d1, ep12);
        //Using anticipation here
        //optiin 1:
        //infd = GetInfDirection(d1, d2*ped2->GetV0Norm(), ep12, Distance);
        //option 2: Better than option 1, especially in crossing scenarios.
        infd = GetInfDirection(d1, e2*ped2->GetV0Norm(), ep12, Distance);
        double condition3 = d1.ScalarProduct(e2);// ped2 move in the same direction of ped1's e0;		
        if (GetAttracForce() && condition1 > 0 && condition2 > 0 && condition3 > 0 && S_Gap < 0 && dist>0)
        {
            double R_dist = dist + Dis_Gap;
            R_dist = R_dist < 0 ? 0 : R_dist;
            double R_ij = _aPed * exp((-R_dist) / _DPed);
            F_rep = infd * (-R_ij);
        }
        else
        {
            double R_dist = dist - Dis_Gap;
            R_dist = R_dist < 0 ? 0 : R_dist;
            double R_ij = _aPed * exp((-R_dist) / _DPed);
            F_rep = infd * R_ij;
        }
    }

    // I don't know how to judge it. It may be useful. I don't like it, can be removed.
    if (VerticalNoInf == 1)
    {
        F_rep = F_rep * abs(d1.ScalarProduct(d2));
    }
    return F_rep;
}

Point AGCVMModel::ForceRepPedPush(Pedestrian * ped1, Pedestrian * ped2, Building * building, int periodic) const
{
    Point F_rep(0.0, 0.0);
    Point p1 = ped1->GetPos();
    Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();
    Point distp12 = p2 - p1;
    double Distance = distp12.Norm();
    Point ep12;
    if (Distance >= J_EPS)
    {
        ep12 = distp12.Normalized();
    }
    else
    {
        printf("ERROR: \tIn AGCVMModel::forcePedPed() ep12 can not be calculated!!!\n");
        printf("ped1 %d  ped2 %d\n", ped1->GetID(), ped2->GetID());
        printf("ped1 at (%f, %f), ped2 at (%f, %f)\n", p1._x, p1._y, p2._x, p2._y);
        exit(EXIT_FAILURE);
    }
    Point ei = ped1->GetMoveDirection();
    Point ei2 = ped2->GetMoveDirection();

    // Distance between boundaries
    double dist = Distance - (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());
    if (dist <= J_EPS)
    {
        // May be a better function needed here.
        double R_ij = _aPed * exp((-dist) / _DPed);
        F_rep = ep12 * (-R_ij);
    }
    return F_rep;
}

Point AGCVMModel::ForceRepRoom(Pedestrian* ped, SubRoom* subroom) const
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
    // and finally the closed doors
    /*
    for (const auto & goal : subroom->GetAllTransitions())
    {
        if (!goal->IsOpen())
        {
            f += ForceRepWall(ped, *(static_cast<Line*>(goal)), centroid, inside);
        }

        //door is open, but it's not my door (has influence)
        int uid1 = goal->GetUniqueID();
        int uid2 = ped->GetExitIndex();
        if((uid1 != uid2) && (goal->IsOpen()==true ))
        {
            f +=  ForceRepWall(ped,*(static_cast<Line*>(goal)), centroid, inside, e0);
        }
    }
    for (const auto & crossline : subroom->GetAllCrossings())
    {
        //door is open, but it's not my door (has influence)
        int uid1 = crossline->GetUniqueID();
        int uid2 = ped->GetExitLine()->GetUniqueID();
        //printf("\npedid=%d, uid1=%d, uid2=%d", ped->GetID(), uid1, uid2);
        if((uid1 != uid2))
        {
            f +=  ForceRepWall(ped,*(static_cast<Line*>(crossline)), centroid, inside);
        }
    }
    */
}

Point AGCVMModel::ForceRepWall(Pedestrian* ped, const Line& w, const Point& centroid, bool inside) const
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
    if (GetGCVMU() == 0)//all
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

    Point inf_direction = GetInfDirection(e0, e_iw);
    double effdis = Distance - ped->GetEllipse().GetBmax();//Using circle now.
    //effdis = effdis < 0 ? 0 : effdis;
    double R_iw = _aWall * exp((-effdis) / _DWall);
    if (GetGCVMU() == 1)
    {
        F_wrep = inf_direction * R_iw;//new method
    }
    else
    {
        F_wrep = e_iw * R_iw*-1;//original method
    }
    return F_wrep;
}

my_pair AGCVMModel::GetSpacing(Pedestrian* ped1, Pedestrian* ped2, int periodic) const
{
    Point ei = ped1->GetMoveDirection();
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
    double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
    double l = GetPushing() == true ? 2 * GetCoreSize() : (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());
    double condition2 = ei.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
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

double AGCVMModel::GetSpacingRoom(Pedestrian* ped, SubRoom* subroom) const
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

double AGCVMModel::GetSpacingWall(Pedestrian* ped, const Line& l) const
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
    double b = GetPushing() ? ped->GetCore() : ped->GetEllipse().GetBmax();
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

double AGCVMModel::OptimalSpeed(Pedestrian* ped, double spacing) const
{
    double v0 = ped->GetV0Norm();
    double T = GetTs();

    double speed = (spacing) / T;
    speed = (speed > 0) ? speed : 0;
    speed = (speed < v0) ? speed : v0;
    return speed;
}


/*----------Functions helpful----------*/
Point AGCVMModel::GetPosPeriodic(Pedestrian* ped1, Pedestrian* ped2) const
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

Point AGCVMModel::GetInfDirection(Point e0, Point ep12) const
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

Point AGCVMModel::GetInfDirection(Point d1, Point d2, Point ep12, double s12) const
{
    Point newep12 = ep12 * s12 + d2 * GetAntiT();
    return GetInfDirection(d1, newep12);
}

void AGCVMModel::UpdatePed(Pedestrian* ped, Point speed, Point direction, double deltaT, int periodic)
{
    Point pos_neu = ped->GetPos() + speed * deltaT;
    Point e0 = ped->GetLastE0();
    Point MD = direction;
    Point MS = speed;
    //backward_movement
    if (direction.ScalarProduct(e0) < 0)
    {
        MD = MD * -1;
        MS = Point(0, 0);
    }
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

/*----------Functions may helpful----------*/
my_pair AGCVMModel::JudgeCollision(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const
{
    double ttc = FLT_MAX;
    double At = GetCoopT();
    //At = 1;
    Point p1 = ped1->GetPos();
    Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();
    Point  e1 = ped1->GetMoveDirection();
    //e1 = e1 * ped1->GetV0Norm();
    Point  e2 = ped2->GetMoveDirection();
    //e2 = e2 * ped2->GetV0Norm();

    Room* room1 = building->GetRoom(ped1->GetRoomID());
    SubRoom* subroom1 = room1->GetSubRoom(ped1->GetSubRoomID());
    Point eo1 = DesireDirection(ped1, room1);
    Room* room2 = building->GetRoom(ped2->GetRoomID());
    SubRoom* subroom2 = room2->GetSubRoom(ped2->GetSubRoomID());
    Point eo2 = DesireDirection(ped2, room2);
    //e1 = eo1;
    //e2 = eo2;
    //e1 = eo1* ped1->GetV0Norm();
    //e2 = eo2 * ped2->GetV0Norm();


    //e1 = ped1->GetV();
    //e2 = ped2->GetV();
    Point distp12 = p2 - p1; //ped1 ---> ped2
    double distance = distp12.Norm();
    Point e1m = ped1->GetMoveDirection();
    Point e2m = ped2->GetMoveDirection();
    double condition1 = e1m.ScalarProduct(distp12);
    //condition1 = e1.ScalarProduct(distp12);
    double condition2 = e2m.ScalarProduct(distp12);
    //condition2 = e2.ScalarProduct(distp12);
    if (0 && condition1 < 0 || condition2 > 0)
    {
        return my_pair(ttc, ped2->GetID());
    }
    JEllipse Eped1 = ped1->GetEllipse();
    JEllipse Eped2 = ped2->GetEllipse();
    double r = ped1->GetEllipse().GetBmax();
    r = 0.2;

    double a = (e1._x - e2._x)*(e1._x - e2._x) + (e1._y - e2._y)*(e1._y - e2._y);
    double b = 2 * ((p1._x - p2._x)*(e1._x - e2._x) + (p1._y - p2._y)*(e1._y - e2._y));
    double c = distance * distance - 4 * r * r;
    double delta = b * b - 4 * a * c;

    // Distance is zero, collision
    /*
    if (dist-0.4 < 0)
    {
        collision = 1;
        return collision;
    }
    */
    if (delta == 0)
    {
        double t = -b / (2 * a);
        if (t > J_EPS && t < At)
        {
            ttc = t;
        }
    }
    else if (delta > 0)
    {
        double sd = sqrt(delta);
        double t1 = (-b - sd) / (2 * a);
        double t2 = (-b + sd) / (2 * a);
        if ((t1 > J_EPS && t1 < At))
        {
            ttc = t1;
        }
        /*
        else if (t1 <= 0 && t2 > J_EPS && t2 < At)
        {
            ttc = t2;
        }
        */
    }
    return my_pair(ttc, ped2->GetID());
}

bool AGCVMModel::Drill(Pedestrian*ped1, vector<Pedestrian*> neighbours, Building* building, SubRoom* subroom, Point e0, int periodic) const
{
    int size = (int)neighbours.size();
    bool drill = true;
    for (int i = 0; i < size; i++)
    {
        Pedestrian* ped2 = neighbours[i];
        SubRoom* subroom2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());
        Point p2 = periodic ? GetPosPeriodic(ped1, ped2) : ped2->GetPos();
        Point ei = ped1->GetMoveDirection();
        Point distp12 = p2 - ped1->GetPos(); //ped1 ---> ped2
        double Distance = distp12.Norm();
        Point ep12;
        if (Distance >= J_EPS)
        {
            ep12 = distp12.Normalized();
        }
        else
        {
            printf("ERROR: \tin AGCVMModel::Drill() ep12 can not be calculated!!!\n");
            exit(EXIT_FAILURE);
        }

        //Judge conllision
        double condition1 = e0.ScalarProduct(ep12); // < e_i , e_ij > should be positive
        double l = ped1->GetCore() + ped2->GetCore();
        double condition2 = e0.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
        condition2 = (condition2 > 0) ? condition2 : -condition2; // abs
        if ((condition1 >= 0) && (condition2 <= l / Distance))
        {
            drill = false;
            return drill;
        }
    }
    if (DrillRoom(ped1, subroom, e0) == false)
    {
        drill = false;
        return drill;
    }
    for (const auto & subr : subroom->GetNeighbors())
    {
        if (DrillRoom(ped1, subr, e0) == false)
        {
            drill = false;
            return drill;
        }
    }
    return drill;
}

bool AGCVMModel::DrillRoom(Pedestrian* ped, SubRoom* subroom, Point e0) const
{
    bool drill = true;
    for (const auto & wall : subroom->GetAllWalls())
    {
        if (DrillWall(ped, e0, wall) == false)
        {
            drill = false;
            return drill;
        }
    }
    //then the obstacles
    for (const auto & obst : subroom->GetAllObstacles())
    {
        for (const auto & wall : obst->GetAllWalls())
        {
            if (DrillWall(ped, e0, wall) == false)
            {
                drill = false;
                return drill;
            }
        }
    }
    //and finally the closed doors
    for (const auto & goal : subroom->GetAllTransitions())
    {
        if (!goal->IsOpen())
        {
            if (DrillWall(ped, e0, *(static_cast<Line*>(goal))) == false)
            {
                drill = false;
                return drill;
            }
        }
    }
}

bool AGCVMModel::DrillWall(Pedestrian* ped, Point e0, const Line& l) const
{
    Point pp = ped->GetPos();
    Point pt = l.ShortestPoint(ped->GetPos());
    Point p1 = l.GetPoint1();
    Point p2 = l.GetPoint2();
    Point dist = pt - pp;
    Point e0_vertical = e0.Rotate(0, 1);
    //--------------------------------------------------------
    double b = ped->GetCore();
    //--------------------------------------------------------
    Point A1 = pp + e0_vertical * b;
    Point A2 = pp - e0_vertical * b;
    Point p1_A1 = p1 - A1;
    Point p2_A1 = p2 - A1;
    Point p12 = p1 - p2;
    double A1_result1 = e0.CrossProduct(p1_A1);
    double A1_result2 = e0.CrossProduct(p2_A1);
    double result3 = e0.ScalarProduct(dist);
    Point p1_A2 = p1 - A2;
    Point p2_A2 = p2 - A2;
    double A2_result1 = e0.CrossProduct(p1_A2);
    double A2_result2 = e0.CrossProduct(p2_A2);
    if (result3 <= 0)
    {
        return true;
    }
    if (A1_result1*A1_result2 > 0 && A2_result1*A2_result2 > 0 && A1_result1*A2_result1 > 0 && A1_result2*A2_result2 > 0)
    {
        return true;
    }
    double effdis = dist.Norm() - b;
    double cosangle = dist.ScalarProduct(e0) / (dist.Norm()*e0.Norm());
    if (cosangle < J_EPS)
    {
        return true;
    }
    return false;
}

Point AGCVMModel::CorrectD(Pedestrian *ped, Point d_direction, SubRoom* subroom) const
{
    d_direction = d_direction.Normalized();
    Point direction = d_direction;
    for (const auto & wall : subroom->GetAllWalls())
    {
        direction = CorrectDWall(ped, d_direction, wall);
        if ((d_direction - direction).Norm() > J_EPS)
        {
            return direction;
        }
    }
    //then the obstacles
    for (const auto & obst : subroom->GetAllObstacles())
    {
        for (const auto & wall : obst->GetAllWalls())
        {
            direction = CorrectDWall(ped, d_direction, wall);
            if ((d_direction - direction).Norm() > J_EPS)
            {
                return direction;
            }
        }
    }
    //and finally the closed doors
    for (const auto & goal : subroom->GetAllTransitions())
    {
        if (!goal->IsOpen())
        {
            direction = CorrectDWall(ped, d_direction, *(static_cast<Line*>(goal)));
            if ((d_direction - direction).Norm() > J_EPS)
            {
                return direction;
            }
        }
    }
    return direction;
}

Point AGCVMModel::CorrectDWall(Pedestrian *ped, Point d_direction, const Line& l) const
{
    Point direction = d_direction;
    Point pt = l.ShortestPoint(ped->GetPos());
    Point dist = pt - ped->GetPos(); // ped ---> wall
    double Distance = dist.Norm(); //Distance between the center of pedestrian and walls
    Point e_iw;
    if (Distance > J_EPS)
    {
        e_iw = dist.Normalized();
    }
    else
    {
        printf("ERROR: \tin AGCVMModel::ForceRepWall() eiw can not be calculated!!!\n");
        exit(EXIT_FAILURE);
    }
    double radius = ped->GetEllipse().GetBmax();
    double costheta = d_direction.ScalarProduct(e_iw);
    if (Distance - radius <= J_EPS && costheta > 0)
    {
        direction = d_direction - e_iw * costheta;
    }
    return direction.Normalized();
}

/*----------Functions not use now----------*/
bool AGCVMModel::ReArrange(const vector< Pedestrian* >& allPeds_ini, vector< Pedestrian* >& allPeds, Building* building)
{
    int nsize = allPeds_ini.size();
    vector<inf_pair> inf_ped = vector<inf_pair>();
    inf_ped.reserve(nsize);
    for (int p = 0; p < nsize; p++)
    {
        Pedestrian* ped = allPeds_ini[p];
        Room* room = building->GetRoom(ped->GetRoomID());
        const Point target = _direction->GetTarget(room, ped);
        double distance = (target - ped->GetPos()).Norm();
        double position = target._x - distance;
        inf_ped.push_back(inf_pair(p, position));
    }
    for (int p = 0; p < nsize; p++)
    {
        for (int q = 0; q < nsize - 1; q++)
        {
            if (std::get<1>(inf_ped[q]) < std::get<1>(inf_ped[q + 1]))
            {
                inf_ped[q].swap(inf_ped[q + 1]);
            }
        }
    }
    for (int p = 0; p < nsize; p++)
    {
        int number = std::get<0>(inf_ped[p]);
        allPeds.push_back(allPeds_ini[number]);
    }
    return true;
}

