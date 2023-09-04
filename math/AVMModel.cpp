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
    bool circle, double Dtime, double Dangle, double Dweight)
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

    //circle_antipode
    _Circle_Antipode = circle;
    _DAntiTime = Dtime;
    _DAngleStep = Dangle;
    _DWeight = Dweight;
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
        //CircleExperiment
        if (_Circle_Antipode)
        {
            Point target = Point(0, 0) - ped->GetPos();
            ped->SetAlwaysTarget(target);
            Point d = target - ped->GetPos();
            double dist = d.Norm();
            cosPhi = d._x / dist;
            sinPhi = d._y / dist;
            ped->InitV0(target);

            JEllipse E = ped->GetEllipse();
            E.SetCosPhi(cosPhi);
            E.SetSinPhi(sinPhi);
            ped->SetEllipse(E);
            continue;
        }
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

        // Just test for detour，it will be interesting --

        //Debug--------------------------------------------------------------------------------------------------------------------
        if (ped1->GetID() == -1)
            printf("\nCurrent Time: %f.\n", current);
        //Debug--------------------------------------------------------------------------------------------------------------------

        // Detour is here, I don't know if it is useful
        IniDirection = DetourDirection(ped1, room1, neighbours, periodic);
        //---------------------------------------------------------

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
                    //repPed += ForceRepPedGCVM(ped1, ped2, building, periodic); // here anticipation seems not work, so just have a try
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


        // Debug------------------------------------------------------
        if (ped1->GetID() == -1)
        {
            printf("The IniDirection of ped %d is (%.3f, %.3f), and repPed is (%.3f, %.3f), repWall is (%.3f, %.3f).\n",
                ped1->GetID(), IniDirection._x, IniDirection._y, repPed._x, repPed._y, repWall._x, repWall._y);
        }
        //---------------------------------------------------------------


        // test----------------------------------------------------
        //dDirection = IniDirection + repWall;
        //---------------------------------------------------------
        switch (GetModel())
        {
        case 0:
            direction = dDirection;
            break;
        case 1:
        case 2:
            Point aDirection = ped1->GetMoveDirection();
            double angleTau = GetTd();
            // The new method for turning direction-------------------------------------------------------
            dDirection = dDirection.Normalized();
            /*
            double Theta = 0;
            double cpTemp1 = aDirection.CrossProduct(IniDirection);
            double cpTemp2 = dDirection.CrossProduct(IniDirection);
            // actual direction and desired directions are on the same side of the target, so turning directly
            if (cpTemp1*cpTemp2 >= 0)
            {
                double tempScalar = dDirection.ScalarProduct(aDirection);
                tempScalar = tempScalar > 1 - J_EPS ? 1 - J_EPS : tempScalar;
                tempScalar = tempScalar < -1 + J_EPS ? -1 + J_EPS : tempScalar;
                Theta = acos(tempScalar);
            }
            else // antual direction and desired direction are not on the same side of the target, so turning to the direction of target
            {
                double tempScalar1 = dDirection.ScalarProduct(IniDirection);
                tempScalar1 = tempScalar1 > 1 - J_EPS ? 1 - J_EPS : tempScalar1;
                tempScalar1 = tempScalar1 < -1 + J_EPS ? -1 + J_EPS : tempScalar1;
                double tempScalar2 = aDirection.ScalarProduct(IniDirection);
                tempScalar2 = tempScalar2 > 1 - J_EPS ? 1 - J_EPS : tempScalar2;
                tempScalar2 = tempScalar2 < -1 + J_EPS ? -1 + J_EPS : tempScalar2;
                Theta = acos(tempScalar1) + acos(tempScalar2);
            }
            double ThetaTemp = Theta * deltaT / angleTau;
            if (aDirection.CrossProduct(dDirection) > 0)// dDireciton is on the left side of aDirection
            {
                direction = aDirection.Rotate(cos(ThetaTemp), sin(ThetaTemp)); //then turn left
            }
            else
            {
                direction = aDirection.Rotate(cos(ThetaTemp), sin(-ThetaTemp)); // then turn right
            }
            //------------------------------------------------------------------------------
            */
            // orignal method----------------------------------------------------------
            Point AccTu = Point(0, 0);
            AccTu = (dDirection.Normalized() - aDirection) / angleTau;
            direction = aDirection + AccTu * deltaT;
            /*
            if (IniDirection.ScalarProduct(dDirection)* IniDirection.ScalarProduct(aDirection) < 0 && aDirection.ScalarProduct(dDirection) < -0.9)
            {
                //direction = dDirection.Normalized();
            }
            */
            /*
            Point target = this->GetDirection()->GetTarget(room1, ped1);
            Point desiredDirection;
            if (_Circle_Antipode)
            {
                target = ped1->GetAlwaysTarget();
                double dist = (target - p1).Norm();
                if (dist < 0.5)
                {
                    direction = dDirection.Normalized();
                }
            }
            */
            //---------------------------------------------------------------------------------


            //Debug---------------------------------------------------------------------------------
            if (ped1->GetID() == -1)
            {
                printf("aDirection (%.3f,%.3f), dDirection (%.3f,%.3f),  direction (%.3f,%.3f)\n",
                    aDirection._x, aDirection._y, dDirection._x, dDirection._y, direction.Normalized()._x, direction.Normalized()._y);
            }
            //------------------------------------------------------------------------------------------------------------

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
        speed = ei * OptimalSpeed(ped1, spacing);
        resultAcc.push_back(speed);

        // Debug
        if (ped1->GetID() == -1)
            printf("spacing is %.3f, speed is %.3f.\n", spacing, speed.Norm());
        // Debug
    }

    //Update everything
    for (int p = start; p <= end; ++p)
    {
        Pedestrian* ped = allPeds[p];
        Point vNeu = resultAcc[p];
        Point dirNeu = resultDir[p];
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

    if (_Circle_Antipode)
    {
        target = ped->GetAlwaysTarget();
    }
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

// Check if pedestrian need to detour
Point AVMModel::DetourDirection(Pedestrian * ped, Room * room, vector<Pedestrian*> neighbours, int periodic) const
{
    // 
    int debug_id = -1;
    // Anticipation time for detour part
    double detourAntiT = GetDAntiTime();
    //printf("DAntiTime=%.2f\n", detourAntiT);
    /*
    // NOTUSE: Set the patience of pedestrians------------------------------------------------------------------------
    double maxPatience = 0.0; // patient time (the parameter is not used now)
    double currentPatience = ped->GetDetourPatient();
    ped->SetDetourPatient(currentPatience + 0.05);
    //printf("ped %d have %f second to wait.\n", ped->GetID(), ped->GetPatienceTime());
    //--------------------------------------------------------------------------------------------------------------------------------
    */

    // First get the target point----------------------------------------------------------------------------------------------
    Point target = this->GetDirection()->GetTarget(room, ped);
    Point desiredDirection = Point(0, 0);
    if (_Circle_Antipode)
    {
        target = ped->GetAlwaysTarget();
    }
    //----------------------------------------------------------------------------------------------------------------------------

    // Second get the original desired direction-----------------------------------------------------------------------
    const Point pos = ped->GetPos();
    double dist = (target - pos).Norm(); // this is the distance to the target
    if (dist > J_EPS_GOAL)
    {
        Point NewE0;
        NewE0 = (target - pos).Normalized();
        ped->SetLastE0(NewE0);
    }
    desiredDirection = ped->GetLastE0();
    /*
    if (dist <= 0.5) // already very close to the target, so don't need to detour
    {
        return desiredDirection;
    }
    */
    //---------------------------------------------------------------------------------------------------------------------------------

    /*
    // NOTUSE: Check if pedestrian already in detour and still have patience----------------------------------------
    if (ped->GetDetour() && currentPatience < maxPatience)
    {
        //printf("Ped %d is already in detour.\n", ped->GetID());
        //printf("Ped %d still have %f patience.\n", ped->GetID(), maxPatience - currentPatience);
        Point center = ped->GetDetourCenter();
        Point c2p = (pos - center).Normalized();
        Point c2t = (target - center).Normalized();
        if (c2p.CrossProduct(c2t) < 0) // c2p on the right of c2t
        {
            desiredDirection = c2p.Rotate(0, -1); // turn right
        }
        else
        {
            desiredDirection = c2p.Rotate(0, 1); // turn left
        }
        //printf("Ped %d detourAngle is %f.\n", ped->GetID(), ped->GetDetourAngle());
        //printf("Ped %d center is (%f, %f).\n", ped->GetID(), center._x, center._y);
        //printf("Ped %d c2p is (%f, %f).\n", ped->GetID(), c2p._x, c2p._y);
        //printf("Ped %d c2t is (%f, %f).\n", ped->GetID(), c2t._x, c2t._y);
        //printf("Ped %d desired direction is (%f, %f).\n", ped->GetID(), desiredDirection._x, desiredDirection._y);
        return desiredDirection;
    }
    // printf("ped %d current patience is %f.\n", ped->GetID(), currentPatience);
    // although pedestains are not in the detour, but still have patience
    if (currentPatience < maxPatience)
    {
        //printf("ped %d have %f second to wait.\n", ped->GetID(), ped->GetDetourPatient());
        return  desiredDirection;
    }
    // else, which pedestrian don't have patience, want to detour
    ped->SetDetour(false);
    ped->SetDetourPatient(0);
    */
    //-------------------------------------------------------------------------------------------------------------------------------------------------

    // Check if the current direction good enough-------------------------------------------------------------------------------------------------------------------
    // Debug
    if (ped->GetID() == debug_id)
        printf("The position of ped %d is (%f,%f), and the target is (%f,%f).\n", ped->GetID(), pos._x, pos._y, target._x, target._y);
    // Debug
    vector<double> pedOnway;
    Point UsefulDistance = Point(dist, 0); // x: free distance, y: congestion distance
    for (int i = 0; i < neighbours.size(); i++)
    {
        Pedestrian * ped2 = neighbours[i];
        // Judge if ped2 is on ped1's moving path
        Point pos2 = periodic ? GetPosPeriodic(ped, ped2) : ped2->GetPos();
        // if we need Anticipation ? we need to check it
        pos2 = pos2 + ped2->GetV() * detourAntiT;
        Point distp12 = pos2 - pos; //ped1 ---> ped2
        double Distance = distp12.Norm();
        Point ep12 = distp12.Normalized();
        double condition1 = desiredDirection.ScalarProduct(ep12);
        double l = ped->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax();
        // Rotate(0,1) turning left for 90 degree [Rotate(double ctheta, double stheta)]
        double condition2 = desiredDirection.Rotate(0, 1).ScalarProduct(ep12);
        condition2 = condition2 > 0 ? condition2 : -condition2;
        if ((condition1 > 0) && (condition2 < l / Distance) && (Distance < dist))
        {
            // Debug--------------------------------------------------------------------------------------
            //if (ped->GetID() == -1)
                //printf("Stupid ped %d is on the way, pos is (%.3f, %.3f).\n", ped2->GetID(), ped2->GetPos()._x, ped2->GetPos()._y);
            // Debug--------------------------------------------------------------------------------------
            pedOnway.push_back(ped2->GetV().ScalarProduct(desiredDirection)); // the projection of speed on the desired direction
            /*
            if (Distance < UsefulDistance._x)
            {
                UsefulDistance._x = Distance;
                UsefulDistance._y = dist - Distance;
                //Debug--------------------------------------------------------------------------------------------------------------------
                //printf("ped %d distance1 %f, distance2 %f.\n", ped->GetID(), UsefulDistance._x, UsefulDistance._y);
                //Debug--------------------------------------------------------------------------------------------------------------------
            }
            */
        }
    }
    // this current path is free
    double  average = ped->GetV0Norm();
    if (!pedOnway.empty())
    {
        double sum = 0;
        for (vector<double>::iterator it = pedOnway.begin(); it != pedOnway.end(); it++)
        {
            sum += *it;
        }
        average = sum / pedOnway.size(); // average speed
    }
    else
    {
        return desiredDirection;
    }

    // I have to say this function is very important, maybe there is a better way
    // double punish = 500000;
    // double heurist = average > 0 ? UsefulDistance._x / ped->GetV0Norm() + punish * UsefulDistance._y / average : FLT_MAX;
    // heurist = (pedOnway.size() + 1) / exp(average) / exp(average);
    double angel2Current = desiredDirection.ScalarProduct(ped->GetMoveDirection()) + 1;
    double heurist = (pedOnway.size() + 1)* dist / angel2Current / exp(average);
    // We propose a better function here 2023.06.02
    heurist = 1 + (angel2Current - 1) + average / ped->GetV0Norm() + exp(-1.0 * pedOnway.size());
    heurist = 1 + average / ped->GetV0Norm() + exp(-1.0 * pedOnway.size());
    double weight = GetDWeight();
    heurist = 1 + weight * average / ped->GetV0Norm();
    // 2023.06.27 a new exponentail function
    heurist = 1 + weight * exp(average - ped->GetV0Norm());
    // 2023.06.28 a double exponentail function
    heurist = 1 + weight * exp(average - ped->GetV0Norm());
    // 2023.06.30 a double exponentail function
    heurist = 1 + exp(weight*(average - ped->GetV0Norm()));
    //heurist = (pedOnway.size() + 1) / exp(average) / angelTemp;
    // 2023.07.31 considering something new
    heurist = weight * (ped->GetV0Norm() - average);

    //Debug--------------------------------------------------------------------------------------------------------------------
    if (ped->GetID() == debug_id)
    {
        printf("ped %d dist is %f, pedOnway is %d,  average is %f, angelTemp is %f,  heurist is %f.\n", ped->GetID(), dist, pedOnway.size(), average, angel2Current, heurist);
        printf("ped %d dist / newDist is %f, exp(-pedOnway.size()) is %f, average / ped->GetV0Norm() is %f, angel2Current is %f,  heurist is %f.\n", ped->GetID(), 1, exp(-1.0 * pedOnway.size()), average / ped->GetV0Norm(), angel2Current - 1, heurist);
    }
    //Debug--------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------------------

    // Find if there is a better detour router--------------------------------------------------------------------------------------------------------------
    double angleStep = GetDAngleStep(); //discrete angle step, also for test (it should be set in the infile, different discrete size)
    //printf("DAngleStep=%.2f\n", angleStep);
    vector<int> angles;
    for (int i = -90; i <= 90; i += angleStep)
    {
        if (i == 0)
            continue;
        angles.push_back(i);
        //printf("%d\n", i);
    }
    for (int i = 0; i < angles.size(); i++)
    {
        pedOnway.clear();
        // 1. Find the center of detour circle
        double angleTarget = atan2(target._y - pos._y, target._x - pos._x);
        double deviationAngle = angles[i] * M_PI / 180 + angleTarget;
        Point da = Point(cos(deviationAngle), sin(deviationAngle));
        double dist2Center = abs(dist / (2 * sin(angles[i] * M_PI / 180)));
        Point dir2Center = angles[i] > 0 ? Point(da._y, -da._x) : Point(-da._y, da._x);
        Point Center = pos + dir2Center * dist2Center;
        /*
        //Debug--------------------------------------------------------------------------------------------------------------------
        printf("deviation angle is %d.\n", angles[i]);
        printf("ped %d angleTarget is %f.\n", ped->GetID(), angleTarget);
        printf("ped %d deviation angle is (%f, %f).\n", ped->GetID(), da._x, da._y);
        printf("ped %d distance to center is %f.\n", ped->GetID(), dist2Center);
        printf("ped %d postion is (%f,%f).\n", ped->GetID(), ped->GetPos()._x, ped->GetPos()._y);
        printf("ped %d direction to center is (%f,%f).\n", ped->GetID(), dir2Center._x, dir2Center._y);
        printf("ped %d detour center is (%f,%f).\n", ped->GetID(), Center._x, Center._y);
        //Debug--------------------------------------------------------------------------------------------------------------------
        */
        //Debug--------------------------------------------------------------------------------------------------------------------
        if (ped->GetID() == debug_id)
            printf("deviation angle is %d.\n", angles[i]);
        //Debug--------------------------------------------------------------------------------------------------------------------

        // 2. Check is there any wall on the detour route
        SubRoom* subroom = room->GetSubRoom(ped->GetSubRoomID());
        bool intersectionWall = false;
        for (const auto & wall : subroom->GetAllWalls())
        {
            // extend the wall first to the size of pedestrian
            auto extendWall = wall.Enlarge(ped->GetEllipse().GetBmax());
            // then check the intersection points of wall and detour circles
            vector<Point> inters;
            extendWall.pointIntersectionWithCircle(inters, Center, dist2Center);
            // check if the intersection points are on the detour path (arc)
            for (auto interp : inters)
            {
                //printf("IMPORTANT: interpoint is (%f,%f).\n", interp._x, interp._y);
                double tempc1 = (interp - Center).CrossProduct(pos - Center);
                double tempc2 = (interp - Center).CrossProduct(target - Center);
                if (tempc1*tempc2 <= 0)
                {
                    intersectionWall = true;
                    //printf("ERROR: the intersection point is (%f,%f).\n", interp._x, interp._y);
                    break;
                }
            }
            if (intersectionWall == true)
            {
                break;
            }
        }
        if (intersectionWall == true)
        {
            //printf("ERROR: the detour direciton is bad.\n");
            continue;
        }
        // 3. Check the heurist function value of this direction
        double newDist = dist2Center * abs(2 * angles[i] * M_PI / 180); // the distace of detour path
        UsefulDistance = Point(newDist, 0);
        for (int j = 0; j < neighbours.size(); j++)
        {
            Pedestrian* ped2 = neighbours[j];
            Point pos2 = periodic ? GetPosPeriodic(ped, ped2) : ped2->GetPos();
            // Anticipation here, do we need it?
            pos2 = ped2->GetPos() + ped2->GetV() * detourAntiT;
            Point Ped2Center = pos2 - Center;
            bool condition1 = Ped2Center.Norm() > dist2Center - 2 * ped->GetEllipse().GetBmax() && Ped2Center.Norm() < dist2Center + 2 * ped->GetEllipse().GetBmax();
            Point midPoint = (target + pos) / 2;
            bool condition2 = Ped2Center.ScalarProduct(midPoint - Center) > (pos - Center).ScalarProduct(midPoint - Center);
            Point c2t = (target - Center).Normalized();
            Point tangleDirection = Point(0, 0);
            // The tangent direction of detour path
            if (Ped2Center.CrossProduct(c2t) < 0) // c2p on the right of c2t
            {
                tangleDirection = Ped2Center.Rotate(0, -1); // turn right
            }
            else
            {
                tangleDirection = Ped2Center.Rotate(0, 1); // turn left
            }
            if (condition1&&condition2)
            {
                /*
                double angleTemp = abs(acos(Ped2Center.ScalarProduct(c2t) / (Ped2Center.Norm()*c2t.Norm())));
                double distTemp = dist2Center * angleTemp;
                if (distTemp < UsefulDistance._x)
                {
                    UsefulDistance._x = distTemp;
                    UsefulDistance._y = newDist - distTemp;
                }
                */
                //  the projection of speed on the tagent direction of detour path
                // Debug------------------------------------------------------------------------------
                //if (ped->GetID() == -1)
                    //printf("Stupid ped %d is on the way, pos is (%.3f, %.3f).\n", ped2->GetID(), ped2->GetPos()._x, ped2->GetPos()._y);
                // Debug------------------------------------------------------------------------------
                pedOnway.push_back(ped2->GetV().ScalarProduct(tangleDirection.Normalized()));
            }
        }
        if (pedOnway.empty())
        {
            average = ped->GetV0Norm();
        }
        else
        {
            double sum = 0;
            //average speed;
            for (vector<double>::iterator it = pedOnway.begin(); it != pedOnway.end(); it++)
            {
                sum += *it;
            }
            average = sum / pedOnway.size();
            //average = average / pedOnway.size();
        }
        // This function has to be checked, very important
        //double newHeurist = average > 0 ? UsefulDistance._x / ped->GetV0Norm() + punish * UsefulDistance._y / average : FLT_MAX;
        //newHeurist = (pedOnway.size() + 1) / exp(average) / exp(average);
        angel2Current = da.ScalarProduct(ped->GetMoveDirection()) + 1;
        //newHeurist = (pedOnway.size() + 1)* newDist / angelTemp / exp(average);
        double newHeurist = (pedOnway.size() + 1)* newDist / exp(average) / angel2Current;
        // peopose a new heurist function here
        newHeurist = dist / newDist + (angel2Current - 1) + average / ped->GetV0Norm() + exp(-1.0*pedOnway.size());
        newHeurist = dist / newDist + average / ped->GetV0Norm() + exp(-1.0*pedOnway.size());
        newHeurist = dist / newDist + weight * average / ped->GetV0Norm();
        // 2023.06.27 a new exponentail function
        newHeurist = dist / newDist + weight * exp(average - ped->GetV0Norm());
        // 2023.06.28 a double exponentail function
        newHeurist = exp(dist - newDist) + weight * exp(average - ped->GetV0Norm());
        // 2023.06.29 a double exponentail function
        newHeurist = dist / newDist + exp(weight*(average - ped->GetV0Norm()));
        // 2023.07.31 considering something new
        newHeurist = (newDist - dist) + weight * (ped->GetV0Norm() - average);
        //Debug--------------------------------------------------------------------------------------------------------------------
        if (ped->GetID() == debug_id)
        {
            printf("ped %d newdist is %f, pedOnway is %d, newaverage is %f, angelTemp is %f,  newheurist is %f.\n", ped->GetID(), newDist, pedOnway.size(), average, angel2Current, newHeurist);
            printf("ped %d dist / newDist is %f, exp(-pedOnway.size()) is %f, average / ped->GetV0Norm() is %f, angel2Current is %f,  newheurist is %f.\n", ped->GetID(), dist / newDist, exp(-1.0*pedOnway.size()), average / ped->GetV0Norm(), angel2Current - 1, newHeurist);
        }
        //Debug--------------------------------------------------------------------------------------------------------------------
 // 4. Check if this path is better than others
        if (newHeurist < heurist)
        {
            ped->SetDetour(true);
            ped->SetDetourAngle(angles[i]);
            ped->SetDetourCenter(Center);
            heurist = newHeurist;
            desiredDirection = da;
        }
        else if (newHeurist == heurist) // for the same heuristic value, choose one randomly
        {
            int random = rand() % 10000;
            //Debug--------------------------------------------------------------------------------------------------------------------
            if (ped->GetID() == -1)
                printf("random value is %d.\n", random);
            //Debug--------------------------------------------------------------------------------------------------------------------
            if (random > 5000)
            {
                ped->SetDetour(true);
                ped->SetDetourAngle(angles[i]);
                ped->SetDetourCenter(Center);
                heurist = newHeurist;
                desiredDirection = da;
            }
        }
    }

    //Debug--------------------------------------------------------------------------------------------------------------------
    if (ped->GetID() == debug_id)
    {
        printf("ped %d detour angle is %f, detour center is (%f,%f).\n", ped->GetID(), ped->GetDetourAngle(), ped->GetDetourCenter()._x, ped->GetDetourCenter()._y);
        printf("ped %d desired direction is (%f,%f).\n", ped->GetID(), desiredDirection._x, desiredDirection._y);
        if (ped->GetDetour())
        {
            printf("ped %d is in detour jetzt.\n", ped->GetID());
        }
    }
    //Debug--------------------------------------------------------------------------------------------------------------------

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
    // for testing the detour model, pedestrians only influence others when they are too close
    /*
    if (Distance > 0.5)
    {
        //return FRep;
    }
    */
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
    // if pedestrian already arrive target, then don't move anymore
    if (_Circle_Antipode)
    {
        bool alreadyArriving = (ped->GetPos() - ped->GetAlwaysTarget()).Norm() < ped->GetEllipse().GetBmax();
        if (alreadyArriving)
        {
            ped->SetV(Point(0, 0));
            return;
        }
    }
    //--------------------------------------------------------------------------------
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

