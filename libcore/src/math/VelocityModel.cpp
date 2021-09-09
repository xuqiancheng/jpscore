/**
 * \file        VelocityModel.cpp
 * \date        Aug. 07, 2015
 * \version     v0.7
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
 * 3. Velocity Model: Tordeux2015
 *
 *
 **/
#include "VelocityModel.h"

#include "direction/walking/DirectionStrategy.h"
#include "general/Logger.h"
#include "general/OpenMP.h"
#include "geometry/SubRoom.h"
#include "geometry/Wall.h"
#include "neighborhood/NeighborhoodSearch.h"
#include "pedestrian/Pedestrian.h"

double xRight = 26.0;
double xLeft  = 0.0;
double cutoff = 2.0;


VelocityModel::VelocityModel(
    std::shared_ptr<DirectionManager> dir,
    double aped,
    double Dped,
    double awall,
    double Dwall,
    int covid,
    int fType,
    int maxNum,
    double checkTRate,
    double checkDRate,
    double SocialDist)
{
    _direction      = dir;
    _aPed           = aped;
    _DPed           = Dped;
    _aWall          = awall;
    _DWall          = Dwall;
    _isCovid        = covid;
    _fType          = fType;
    _maxNumber      = maxNum;
    _checkTimeRate  = checkTRate;
    _checkDisRate   = checkDRate;
    _socialDistance = SocialDist;
    _BestCheckout   = "checkout1";
    _checkoutFull   = false;
}


VelocityModel::~VelocityModel() {}

bool VelocityModel::Init(Building * building)
{
    _direction->Init(building);

    const std::vector<Pedestrian *> & allPeds = building->GetAllPedestrians();
    size_t peds_size                          = allPeds.size();
    for(unsigned int p = 0; p < peds_size; p++) {
        Pedestrian * ped = allPeds[p];
        double cosPhi, sinPhi;
        //a destination could not be found for that pedestrian
        int ped_is_waiting = 1; // quick and dirty fix
        // we should maybe differentiate between pedestrians who did not find
        // routs because of a bug in the router and these who simplyt just want
        // to wait in waiting areas
        int res = ped->FindRoute();
        if(!ped_is_waiting && res == -1) {
            std::cout << ped->GetID() << " has no route\n";
            LOG_ERROR(
                "VelocityModel::Init() cannot initialise route. ped {:d} is deleted in Room "
                "%d %d.\n",
                ped->GetID(),
                ped->GetRoomID(),
                ped->GetSubRoomID());
            building->DeletePedestrian(ped);
            // TODO KKZ track deleted peds
            p--;
            peds_size--;
            continue;
        }

        // TODO
        // HERE every ped should have a navline already
        //


        Point target = Point(0, 0);
        if(ped->GetExitLine())
            target = ped->GetExitLine()->ShortestPoint(ped->GetPos());
        else {
            std::cout << "Ped " << ped->GetID() << " has no exit line in INIT\n";
        }
        if(_isCovid == 1) {
            Room * room       = building->GetRoom(ped->GetRoomID());
            NavLine * NewExit = RandomExitLine(ped, room);
            ped->SetExitLine(NewExit);
            ped->SetExitIndex(NewExit->GetUniqueID());
            target = _direction->GetTarget(room, ped);
        }
        Point d     = target - ped->GetPos();
        double dist = d.Norm();
        if(dist != 0.0) {
            cosPhi = d._x / dist;
            sinPhi = d._y / dist;
        } else {
            LOG_ERROR("allPeds::Init() cannot initialise phi! dist to target is 0");
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

void VelocityModel::ComputeNextTimeStep(
    double current,
    double deltaT,
    Building * building,
    int periodic)
{
    // collect all pedestrians in the simulation.
    const std::vector<Pedestrian *> & allPeds = building->GetAllPedestrians();
    _BestCheckout                             = BestCheckout(building);
    std::vector<Pedestrian *> pedsToRemove;
    pedsToRemove.reserve(500);
    unsigned long nSize;
    nSize = allPeds.size();

    int nThreads = omp_get_max_threads();

    int partSize;
    partSize = ((int) nSize > nThreads) ? (int) (nSize / nThreads) : (int) nSize;
    if(partSize == (int) nSize)
        nThreads = 1; // not worthy to parallelize


//TODO richtig parallelisieren!
#pragma omp parallel default(shared) num_threads(nThreads)
    {
        std::vector<Point> result_acc = std::vector<Point>();
        result_acc.reserve(nSize);
        std::vector<my_pair> spacings = std::vector<my_pair>();
        spacings.reserve(nSize);             // larger than needed
        spacings.push_back(my_pair(100, 1)); // in case there are no neighbors
        const int threadID = omp_get_thread_num();

        int start = threadID * partSize;
        int end;
        end = (threadID < nThreads - 1) ? (threadID + 1) * partSize - 1 : (int) (nSize - 1);
        for(int p = start; p <= end; ++p) {
            Pedestrian * ped  = allPeds[p];
            Room * room       = building->GetRoom(ped->GetRoomID());
            SubRoom * subroom = room->GetSubRoom(ped->GetSubRoomID());
            Point repPed      = Point(0, 0);

            double virus   = 0;
            double contact = 0;

            //Count the time pedestrian in shop, if pedestrian not in the shop, set it as 0.
            if(InRightRoom(room, {"supermarket"})) {
                double timeInshop = ped->GetTimeInShop() + deltaT;
                ped->SetTimeInShop(timeInshop);
            } else {
                ped->SetTimeInShop(0);
            }

            //Pedestrian is the first one in the chekout, start couting checkout time
            if(InRightRoom(room, {"checkout1", "checkout2", "checkout3"}) &&
               ped->GetExitLine()->DistTo(ped->GetPos()) < 2 * ped->GetEllipse().GetAmin()) {
                double timeCheckout = ped->GetTimeCheckout() + deltaT;
                ped->SetTimeCheckout(timeCheckout);
            } else {
                ped->SetTimeCheckout(0);
            }

            std::vector<Pedestrian *> neighbours =
                building->GetNeighborhoodSearch().GetNeighbourhood(ped);

            int size = (int) neighbours.size();
            for(int i = 0; i < size; i++) {
                Pedestrian * ped1 = neighbours[i];
                if(ped1 == nullptr) {
                    std::cout << "Velocity Model debug: " << size << std::endl;
                }
                //if they are in the same subroom
                Point p1 = ped->GetPos();

                Point p2 = ped1->GetPos();
                //subrooms to consider when looking for neighbour for the 3d visibility
                std::vector<SubRoom *> emptyVector;
                emptyVector.push_back(subroom);
                emptyVector.push_back(
                    building->GetRoom(ped1->GetRoomID())->GetSubRoom(ped1->GetSubRoomID()));
                bool isVisible = building->IsVisible(p1, p2, emptyVector, false);
                if(!isVisible)
                    continue;
                if(ped->GetUniqueRoomID() == ped1->GetUniqueRoomID()) {
                    repPed += ForceRepPed(ped, ped1, periodic);
                    if(_isCovid == 5 &&
                       InRightRoom(room, {"checkout1", "checkout2", "checkout3", "supermarket"})) {
                        virus += VirusContactAmount(ped, ped1);
                        contact += ContactDegree(ped, ped1, _fType);
                    } else if(_isCovid == 1 || _isCovid == 2) {
                        virus += VirusContactAmount(ped, ped1);
                        contact += ContactDegree(ped, ped1, _fType);
                    }
                } else {
                    // or in neighbour subrooms
                    Room * room1  = building->GetRoom(ped1->GetRoomID());
                    SubRoom * sb2 = room1->GetSubRoom(ped1->GetSubRoomID());
                    if(subroom->IsDirectlyConnectedWith(sb2)) {
                        repPed += ForceRepPed(ped, ped1, periodic);
                        if(_isCovid == 5 &&
                           InRightRoom(
                               room, {"checkout1", "checkout2", "checkout3", "supermarket"}) &&
                           InRightRoom(
                               room1, {"checkout1", "checkout2", "checkout3", "supermarket"})) {
                            virus += VirusContactAmount(ped, ped1);
                            contact += ContactDegree(ped, ped1, _fType);
                        } else if(_isCovid == 1 || _isCovid == 2) {
                            virus += VirusContactAmount(ped, ped1);
                            contact += ContactDegree(ped, ped1, _fType);
                        }
                    }
                }
            } // for i
            if(_isCovid > 0) {
                //intergrate virus
                virus = ped->GetVirusContact() + virus * deltaT;
                ped->SetVirusContact(virus);
                ped->SetVirusGet(VirusGetAmount(ped));
                ped->SetProInfect(ProbInfect(ped));
                //intergrate contact
                contact = ped->GetContactDegree() + contact * deltaT;
                ped->SetContactDegree(contact);
            }

            // if the door is closed
            bool IfClose = false;
            if((InRightRoom(room, {"supermarket"}) &&
                ped->GetTimeInShop() < ped->GetMaxTimeInShop()) ||
               (InRightRoom(room, {"waiting"}) && IfMarketFull(building)) ||
               (InRightRoom(
                    room,
                    {
                        "checkout1",
                        "checkout2",
                        "checkout3",
                    }) &&
                ped->GetTimeCheckout() < ped->GetMaxTimeInShop() * _checkTimeRate)) {
                IfClose = true;
            }

            //repulsive forces to walls and closed transitions that are not my target
            Point repWall = ForceRepRoom(allPeds[p], subroom, IfClose);

            // calculate new direction ei according to (6)
            Point direction = e0(ped, room) + repPed + repWall;
            if(_isCovid == 5) {
                if(InRightRoom(room, {"checkout1", "checkout2", "checkout3"})) {
                    direction = e0(ped, room) + repWall + repPed * 0.5;
                }
                if(InRightRoom(room, {"waiting"})) {
                    direction = e0(ped, room) + repWall; //+ repPed;
                }
            }
            for(int i = 0; i < size; i++) {
                Pedestrian * ped1 = neighbours[i];
                // calculate spacing
                // my_pair spacing_winkel = GetSpacing(ped, ped1);
                if(ped->GetUniqueRoomID() == ped1->GetUniqueRoomID()) {
                    spacings.push_back(GetSpacingEllipse(ped, ped1, direction, room, periodic));
                } else {
                    // or in neighbour subrooms
                    SubRoom * sb2 =
                        building->GetRoom(ped1->GetRoomID())->GetSubRoom(ped1->GetSubRoomID());
                    if(subroom->IsDirectlyConnectedWith(sb2)) {
                        spacings.push_back(GetSpacingEllipse(ped, ped1, direction, room, periodic));
                    }
                }
            }
            //TODO get spacing to walls
            //TODO update direction every DT?
            // calculate min spacing
            std::sort(spacings.begin(), spacings.end(), sort_pred());
            double spacing = spacings.size() == 0 ? FLT_MAX : spacings[0].first;
            //============================================================
            // TODO: Hack for Head on situations: ped1 x ------> | <------- x ped2
            if(0 && direction.NormSquare() < 0.5) {
                double pi_half = 1.57079663;
                double alpha   = pi_half * exp(-spacing);
                direction      = e0(ped, room).Rotate(cos(alpha), sin(alpha));
                printf(
                    "\nRotate %f, %f, norm = %f alpha = %f, spacing = %f\n",
                    direction._x,
                    direction._y,
                    direction.NormSquare(),
                    alpha,
                    spacing);
                getc(stdin);
            }
            //============================================================
            Point speed = direction.Normalized() * OptimalSpeed(ped, spacing, room);
            result_acc.push_back(speed);
            /*
            if(ped->GetID() == 92 && spacing < 0.01) {
                direction = direction.Normalized();
                LOG_ERROR(
                    "ID:({:d}),pos=({:f},{:f}),direction({:f},{:f}),spacing({:f},size({:d}))",
                    ped->GetID(),
                    ped->GetPos()._x,
                    ped->GetPos()._y,
                    direction._x,
                    direction._y,
                    spacing,
                    spacings.size());
            }
			*/
            spacings.clear(); //clear for ped p
            /*
            // stuck peds get removed. Warning is thrown. low speed due to jam is omitted.
            if(ped->GetTimeInJam() > ped->GetPatienceTime() &&
               ped->GetGlobalTime() > 10000 + ped->GetPremovementTime() &&
               std::max(ped->GetMeanVelOverRecTime(), ped->GetV().Norm()) < 0.01 &&
               size == 0) // size length of peds neighbour vector
            {
                LOG_WARNING(
                    "ped {:d} with vmean {:f} has been deleted in room {:d}/{:d} after time "
                    "{:f}s (current={:f}",
                    ped->GetID(),
                    ped->GetMeanVelOverRecTime(),
                    ped->GetRoomID(),
                    ped->GetSubRoomID(),
                    ped->GetGlobalTime(),
                    current);
                //TODO KKZ track deleted peds
#pragma omp critical(VelocityModel_ComputeNextTimeStep_pedsToRemove)
                pedsToRemove.push_back(ped);
            }
			*/

        } // for p

#pragma omp barrier
        // update
        for(int p = start; p <= end; ++p) {
            Pedestrian * ped = allPeds[p];

            Point v_neu   = result_acc[p - start];
            Point pos_neu = ped->GetPos() + v_neu * deltaT;

            //Jam is based on the current velocity
            if(v_neu.Norm() >= ped->GetV0Norm() * 0.5) {
                ped->ResetTimeInJam();
            } else {
                ped->UpdateTimeInJam();
            }
            //only update the position if the velocity is above a threshold
            if(v_neu.Norm() >= 0) {
                ped->SetPhiPed();
            }
            ped->SetPos(pos_neu);
            if(periodic) {
                if(ped->GetPos()._x >= xRight) {
                    ped->SetPos(Point(ped->GetPos()._x - (xRight - xLeft), ped->GetPos()._y));
                    //ped->SetID( ped->GetID() + 1);
                }
                if(ped->GetPos()._x <= xLeft) {
                    ped->SetPos(Point(ped->GetPos()._x + (xRight - xLeft), ped->GetPos()._y));
                }
            }
            ped->SetV(v_neu);
        }
    } //end parallel

    // remove the pedestrians that have left the building
    for(unsigned int p = 0; p < pedsToRemove.size(); p++) {
        building->DeletePedestrian(pedsToRemove[p]);
    }
    pedsToRemove.clear();
}

Point VelocityModel::e0(Pedestrian * ped, Room * room) const
{
    Point target;

    if(_direction && ped->GetExitLine()) {
        // target is where the ped wants to be after the next timestep
        target = _direction->GetTarget(room, ped);
    } else { //@todo: we need a model for waiting pedestrians
        std::cout << ped->GetID() << " VelocityModel::e0 Ped has no navline.\n";
        // set random destination
        std::mt19937 mt(ped->GetBuilding()->GetConfig()->GetSeed());
        std::uniform_real_distribution<double> dist(0, 1.0);
        double random_x = dist(mt);
        double random_y = dist(mt);
        Point P1        = Point(ped->GetPos()._x - random_x, ped->GetPos()._y - random_y);
        Point P2        = Point(ped->GetPos()._x + random_x, ped->GetPos()._y + random_y);
        const NavLine L = Line(P1, P2);
        ped->SetExitLine((const NavLine *) &L);
        target = P1;
    }
    Point desired_direction;
    const Point pos = ped->GetPos();
    double dist     = 0.0;
    if(ped->GetExitLine())
        dist = ped->GetExitLine()->DistTo(pos);
    // check if the molified version works
    Point lastE0 = ped->GetLastE0();
    ped->SetLastE0(target - pos);
    SubRoom * subroom = room->GetSubRoom(ped->GetSubRoomID());

    if(_isCovid == 1 && dist < 1) {
        NavLine * NewExit = RandomExitLine(ped, room);
        ped->SetExitLine(NewExit);
        ped->SetExitIndex(NewExit->GetUniqueID());
        target = _direction->GetTarget(room, ped);
    } else if(
        _isCovid == 4 && subroom->GetCaption() == "shop" &&
        ped->GetTimeInShop() <= ped->GetMaxTimeInShop()) {
        if(!ped->GetCounterTogo() || ped->GetCounterTogo()->Contains(ped->GetPos())) {
            int counterSize = (int) subroom->GetAllCounters().size();
            int myseed      = ped->GetGlobalTime() * 100 + ped->GetID();
            srand(myseed);
            int randNum = rand() % counterSize;
            //LOG_ERROR("Test:(ID:{:d}, Value{:d}, time{:d})", ped->GetID(), randNum, myseed);
            Counter * counterTogo = subroom->GetAllCounters()[randNum];
            ped->SetCounterTogo(counterTogo);
        }
        target = ped->GetCounterTogo()->GetCentroid();
    } else if(_isCovid == 5) {
        const NavLine * NewExit = NewExitLineForMarket(ped, room);
        ped->SetExitLine(NewExit);
        ped->SetExitIndex(NewExit->GetUniqueID());
        target = _direction->GetTarget(room, ped);
    }

    if((dynamic_cast<DirectionLocalFloorfield *>(_direction->GetDirectionStrategy().get())) ||
       (dynamic_cast<DirectionSubLocalFloorfield *>(_direction->GetDirectionStrategy().get()))) {
        desired_direction = target - pos;
        if(desired_direction.NormSquare() < 0.25 && !ped->IsWaiting()) {
            desired_direction = lastE0;
            ped->SetLastE0(lastE0);
        }
    } else if(dist > J_EPS_GOAL) {
        desired_direction = ped->GetV0(target);
    } else {
        ped->SetSmoothTurning();
        desired_direction = ped->GetV0();
    }
    return desired_direction;
}


double VelocityModel::OptimalSpeed(Pedestrian * ped, double spacing, Room * room) const
{
    double v0    = ped->GetV0Norm();
    double T     = ped->GetT();
    double speed = (spacing) / T;
    speed        = (speed > 0) ? speed : 0;
    speed        = (speed < v0) ? speed : v0;
    //      (1-winkel)*speed;
    //todo use winkel
    return speed;
}

// return spacing and id of the nearest pedestrian
my_pair
VelocityModel::GetSpacing(Pedestrian * ped1, Pedestrian * ped2, Point ei, int periodic) const
{
    Point distp12 = ped2->GetPos() - ped1->GetPos(); // inversed sign
    if(periodic) {
        double x   = ped1->GetPos()._x;
        double x_j = ped2->GetPos()._x;

        if((xRight - x) + (x_j - xLeft) <= cutoff) {
            distp12._x = distp12._x + xRight - xLeft;
        }
        if((x - xLeft) + (xRight - x_j) <= cutoff) {
            distp12._x = xRight - xLeft - distp12._x;
        }
    }
    double Distance = distp12.Norm();
    double l        = 2 * ped1->GetEllipse().GetBmax();
    Point ep12;
    if(Distance >= J_EPS) {
        ep12 = distp12.Normalized();
    } else {
        LOG_WARNING(
            "VelocityModel::GetSPacing() ep12 can not be calculated! Pedestrians are to close to "
            "each other ({:f})",
            Distance);
        my_pair(FLT_MAX, ped2->GetID());
        exit(EXIT_FAILURE); //TODO
    }

    double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
    double condition2 =
        ei.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
    condition2 = (condition2 > 0) ? condition2 : -condition2; // abs

    if((condition1 >= 0) && (condition2 <= l / Distance))
        // return a pair <dist, condition1>. Then take the smallest dist. In case of equality the biggest condition1
        return my_pair(distp12.Norm(), ped2->GetID());
    else
        return my_pair(FLT_MAX, ped2->GetID());
}
Point VelocityModel::ForceRepPed(Pedestrian * ped1, Pedestrian * ped2, int periodic) const
{
    Point F_rep(0.0, 0.0);
    // x- and y-coordinate of the distance between p1 and p2
    Point distp12 = ped2->GetPos() - ped1->GetPos();

    if(periodic) {
        double x   = ped1->GetPos()._x;
        double x_j = ped2->GetPos()._x;
        if((xRight - x) + (x_j - xLeft) <= cutoff) {
            distp12._x = distp12._x + xRight - xLeft;
        }
        if((x - xLeft) + (xRight - x_j) <= cutoff) {
            distp12._x = xRight - xLeft - distp12._x;
        }
    }

    double Distance = distp12.Norm();
    Point ep12; // x- and y-coordinate of the normalized vector between p1 and p2
    double R_ij;
    //double l = 2 * ped1->GetEllipse().GetBmax();

    if(Distance >= J_EPS) {
        ep12 = distp12.Normalized();
    } else {
        LOG_ERROR(
            "VelocityModel::forcePedPed() ep12 can not be calculated! Pedestrians are too near to "
            "each other (dist={:f}). Adjust <a> value in force_ped to counter this. Affected "
            "pedestrians ped1 {:d} at ({:f},{:f}) and ped2 {:d} at ({:f}, {:f})",
            Distance,
            ped1->GetID(),
            ped1->GetPos()._x,
            ped1->GetPos()._y,
            ped2->GetID(),
            ped2->GetPos()._x,
            ped2->GetPos()._y);
        exit(EXIT_FAILURE); //TODO: quick and dirty fix for issue #158
                            // (sometimes sources create peds on the same location)
    }

    double dist;
    JEllipse Eped1 = ped1->GetEllipse();
    JEllipse Eped2 = ped2->GetEllipse();
    dist           = Eped1.EffectiveDistanceToEllipse(Eped2, &dist);
    //dist           = dist > 0 ? dist : 0;
    R_ij  = -_aPed * exp((-dist) / _DPed);
    F_rep = ep12 * R_ij;

    return F_rep;
} //END Velocity:ForceRepPed()

Point VelocityModel::ForceRepRoom(Pedestrian * ped, SubRoom * subroom, bool ifclose) const
{
    Point f(0., 0.);
    const Point & centroid = subroom->GetCentroid();
    bool inside            = subroom->IsInSubRoom(centroid);
    //first the walls
    for(const auto & wall : subroom->GetAllWalls()) {
        f += ForceRepWall(ped, wall, centroid, inside);
    }

    //then the obstacles
    for(const auto & obst : subroom->GetAllObstacles()) {
        if(obst->Contains(ped->GetPos())) {
            LOG_ERROR(
                "Agent {:d} is trapped in obstacle in room/subroom {:d}/{:d}",
                ped->GetID(),
                subroom->GetRoomID(),
                subroom->GetSubRoomID());
            exit(EXIT_FAILURE);
        } else
            for(const auto & wall : obst->GetAllWalls()) {
                f += ForceRepWall(ped, wall, centroid, inside);
            }
    }

    // and finally the closed doors
    for(const auto & trans : subroom->GetAllTransitions()) {
        if(!trans->IsOpen() || _isCovid == 1) {
            f += ForceRepWall(ped, *(static_cast<Line *>(trans)), centroid, inside);
        }
    }

    if(_isCovid == 3 && ped->GetGroup() != 0) {
        for(const auto & trans : subroom->GetAllTransitions()) {
            f += ForceRepWall(ped, *(static_cast<Line *>(trans)), centroid, inside);
        }
    } else if(_isCovid == 4 && ifclose) {
        for(const auto & trans : subroom->GetAllTransitions()) {
            f += ForceRepWall(ped, *(static_cast<Line *>(trans)), centroid, inside);
        }
    } else if(_isCovid == 5 && ifclose) {
        for(const auto & trans : subroom->GetAllTransitions()) {
            f += ForceRepWall(ped, *(static_cast<Line *>(trans)), centroid, inside);
        }
    }
    return f;
}

Point VelocityModel::ForceRepWall(
    Pedestrian * ped,
    const Line & w,
    const Point & centroid,
    bool inside) const
{
    Point F_wrep = Point(0.0, 0.0);
    Point pt     = w.ShortestPoint(ped->GetPos());

    Point dist       = pt - ped->GetPos(); // x- and y-coordinate of the distance between ped and p
    const double EPS = 0.000;              // molified see Koester2013
    double Distance  = dist.Norm() + EPS;  // distance between the centre of ped and point p
    Point e_iw; // x- and y-coordinate of the normalized vector between ped and pt
    double l = ped->GetEllipse().GetBmax();
    double R_iw;
    double min_distance_to_wall = 0.001; // 10 cm

    if(Distance > min_distance_to_wall) {
        e_iw = dist / Distance;
    } else {
        /*
        LOG_WARNING(
            "Velocity: forceRepWall() ped {:d} [{:f}, {:f}] is too near to the wall [{:f}, "
            "{:f}]-[{:f}, {:f}] (dist={:f})",
            ped->GetID(),
            ped->GetPos()._y,
            ped->GetPos()._y,
            w.GetPoint1()._x,
            w.GetPoint1()._y,
            w.GetPoint2()._x,
            w.GetPoint2()._y,
            Distance);
			*/
        Point new_dist = centroid - ped->GetPos();
        new_dist       = new_dist / new_dist.Norm();
        e_iw           = (inside ? new_dist : new_dist * -1);
    }
    //-------------------------

    const Point & pos = ped->GetPos();
    double distGoal   = 0.0;
    if(ped->GetExitLine())
        distGoal = ped->GetExitLine()->DistToSquare(pos);

    if(distGoal < J_EPS_GOAL * J_EPS_GOAL)
        return F_wrep;
    //-------------------------
    JEllipse Eped1 = ped->GetEllipse();
    double effdis  = Eped1.EffectiveDistanceToLine(w);
    //effdis         = effdis > 0 ? effdis : 0;
    R_iw = -_aWall * exp((-effdis) / _DWall);
    //R_iw   = -_aWall * exp((l - Distance) / _DWall);
    F_wrep = e_iw * R_iw;

    return F_wrep;
}

std::string VelocityModel::GetDescription()
{
    std::string rueck;
    char tmp[1024];

    sprintf(tmp, "\t\ta: \t\tPed: %f \tWall: %f\n", _aPed, _aWall);
    rueck.append(tmp);
    sprintf(tmp, "\t\tD: \t\tPed: %f \tWall: %f\n", _DPed, _DWall);
    rueck.append(tmp);
    return rueck;
}

double VelocityModel::GetaPed() const
{
    return _aPed;
}

double VelocityModel::GetDPed() const
{
    return _DPed;
}


double VelocityModel::GetaWall() const
{
    return _aWall;
}

double VelocityModel::GetDWall() const
{
    return _DWall;
}

NavLine * VelocityModel::RandomExitLine(Pedestrian * ped, Room * room) const
{
    SubRoom * subroom = room->GetSubRoom(ped->GetSubRoomID());
    subroom->GetAllTransitions();
    std::vector<Hline *> allGoals;
    const auto & crossings = subroom->GetAllCrossings();
    allGoals.insert(allGoals.end(), crossings.begin(), crossings.end());
    const auto & transitions = subroom->GetAllTransitions();
    allGoals.insert(allGoals.end(), transitions.begin(), transitions.end());
    const auto & hlines = subroom->GetAllHlines();
    allGoals.insert(allGoals.end(), hlines.begin(), hlines.end());
    int ExitSize         = allGoals.size();
    NavLine * RandomExit = allGoals[rand() % ExitSize];
    while(RandomExit->GetUniqueID() == ped->GetExitLine()->GetUniqueID()) {
        RandomExit = allGoals[rand() % ExitSize];
    }
    return RandomExit;
}

double VelocityModel::VirusContactAmount(Pedestrian * ped1, Pedestrian * ped2) const
{
    double Virus = 0.0;
    int status1  = ped1->GetInfection();
    int status2  = ped2->GetInfection();
    if(status1 == 0 && status2 == 1) {
        Point dist   = ped2->GetPos() - ped1->GetPos();
        double s     = dist.Norm();
        double alpha = ped2->GetCovidAlpha();
        double P     = ped2->GetCovidP();
        double r     = ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax();
        double D     = ped1->GetCovidD();
        Virus        = alpha * P * exp((r - s) / D);
    }
    return Virus;
}

double VelocityModel::VirusGetAmount(Pedestrian * ped) const
{
    double vc    = ped->GetVirusContact();
    double alpha = ped->GetCovidAlpha();
    double k     = ped->GetCovidK();
    double vg    = 1 - exp(-vc * alpha * k);
    return vg;
}

double VelocityModel::ProbInfect(Pedestrian * ped) const
{
    double vg   = ped->GetVirusget();
    double Q    = ped->GetCovidQ();
    double proI = vg * Q;
    return proI;
}

double VelocityModel::ContactDegree(Pedestrian * ped1, Pedestrian * ped2, int func) const
{
    double Contact = 0;
    Point dist     = ped2->GetPos() - ped1->GetPos();
    double s       = dist.Norm();
    double r       = ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax();
    if(func == 1) {
        int K   = 1;
        int D   = 1;
        Contact = K * exp((-s) / D);
    }
    if(func == 2) {
        int K   = 1;
        Contact = K / s;
    }
    if(func == 3) {
        int K     = 1;
        int safeR = 1;
        Contact   = s - safeR > 0 ? 0 : K;
    }
    return Contact;
}

bool VelocityModel::IfMarketFull(Building * building) const
{
    int size                                  = 0;
    const std::vector<Pedestrian *> & allPeds = building->GetAllPedestrians();
    for(auto ped : building->GetAllPedestrians()) {
        Room * room = building->GetRoom(ped->GetRoomID());
        if(InRightRoom(room, {"checkout1", "checkout2", "checkout3", "supermarket"})) {
            size++;
        }
    }
    bool full = false;
    if(size >= _maxNumber) {
        full = true;
    }
    return full;
}

const NavLine * VelocityModel::NewExitLineForMarket(Pedestrian * ped, Room * room) const
{
    const NavLine * oldgoal = ped->GetExitLine();
    Point pos               = ped->GetPos();
    SubRoom * subroom       = room->GetSubRoom(ped->GetSubRoomID());
    bool CanSeeExit         = subroom->IsVisible(pos, oldgoal->GetCentre(), false);
    if(InRightRoom(room, {"supermarket"})) {
        auto counters        = subroom->GetAllCounters();
        Counter * chooseArea = counters[0];
        std::vector<Hline *> allGoals;
        if(((ped->GetTimeInShop() <= ped->GetMaxTimeInShop() || _checkoutFull) &&
            oldgoal->DistTo(pos) < 2 * ped->GetEllipse().GetAmin()) ||
           ped->IsFeelingLikeInJam() || (!CanSeeExit)) {
            for(auto sub : room->GetAllSubRooms()) {
                const auto & crossings = sub.second->GetAllCrossings();
                allGoals.insert(allGoals.end(), crossings.begin(), crossings.end());
                const auto & hlines = sub.second->GetAllHlines();
                allGoals.insert(allGoals.end(), hlines.begin(), hlines.end());
            }
            bool SetNewGoal = false;
            do {
                //int myseed = ped->GetGlobalTime() * 100 + ped->GetID();
                //srand(myseed);
                int randNum           = rand() % (int(allGoals.size()));
                const Hline * newgoal = allGoals[randNum];
                if(newgoal == oldgoal) {
                    continue;
                }
                /*
                double dist = (newgoal->GetCentre() - oldgoal->GetCentre()).Norm();
                if(dist > 11) {
                    continue;
                }
				*/
                SetNewGoal = subroom->IsVisible(pos, newgoal->GetCentre(), false);
                if(SetNewGoal) {
                    //target = newgoal->GetCentre();
                    return newgoal;
                }
            } while(SetNewGoal != true);
        } else if(
            ped->GetTimeInShop() > ped->GetMaxTimeInShop() && chooseArea->Contains(pos) &&
            (!_checkoutFull)) {
            const auto & transitions = subroom->GetAllTransitions();
            for(auto trans : transitions) {
                if(trans->GetCaption() == _BestCheckout) {
                    return trans;
                }
            }
        }
    }
    return oldgoal;
}

std::string VelocityModel::BestCheckout(Building * building)
{
    const std::vector<Pedestrian *> allped = building->GetAllPedestrians();
    std::vector<bool> full(3);
    std::vector<int> num = {false, false, false};
    double end           = 5;
    for(Pedestrian * ped : allped) {
        int ID                  = ped->GetRoomID();
        double x                = ped->GetPos()._x;
        double v                = ped->GetV().Norm();
        double l                = 3 * ped->GetEllipse().GetAmin();
        std::string roomCaption = building->GetRoom(ID)->GetCaption();
        if(roomCaption == "checkout1") {
            num[0]++;
            if((end - x < l)) {
                full[0] = true;
            }
        } else if(roomCaption == "checkout2") {
            num[1]++;
            if((end - x < l)) {
                full[1] = true;
            }
        } else if(roomCaption == "checkout3") {
            num[2]++;
            if((end - x < l)) {
                full[2] = true;
            }
        }
    }
    auto f = std::find(full.begin(), full.end(), false);
    if(f == full.end()) {
        _checkoutFull = true;
    } else {
        _checkoutFull = false;
    }
    auto best = std::min_element(num.begin(), num.end());
    //LOG_ERROR("Test{:d}", *best);
    switch(std::distance(num.begin(), best)) {
        case 0:
            return "checkout1";
        case 1:
            return "checkout2";
        case 2:
            return "checkout3";
    }
}

bool VelocityModel::InRightRoom(Room * room, std::vector<std::string> captions) const
{
    std::string roomName = room->GetCaption();
    auto iter            = find(captions.begin(), captions.end(), roomName);
    if(iter == captions.end()) {
        return false;
    }
    return true;
}

my_pair VelocityModel::GetSpacingEllipse(
    Pedestrian * ped1,
    Pedestrian * ped2,
    Point ei,
    Room * room,
    int periodic) const
{
    double x1 = ped1->GetPos()._x;
    double y1 = ped1->GetPos()._y;

    Point distp12   = ped2->GetPos() - ped1->GetPos(); // inversed sign
    double Distance = distp12.Norm();
    Point ep12;
    if(Distance >= J_EPS) {
        ep12 = distp12.Normalized();
    } else {
        LOG_WARNING("In VelocityModel::GetSPacing() ep12 can not be calculated!!!");
        LOG_WARNING("Pedestrians are too near to each other {:f}.", Distance);
    }
    //calculate effective distance
    JEllipse eped1 = ped1->GetEllipse();
    JEllipse eped2 = ped2->GetEllipse();
    double dist;
    double eff_dist = eped1.EffectiveDistanceToEllipse(eped2, &dist);

    if(InRightRoom(room, {"checkout1", "checkout2", "checkout3"})) {
        double l = ped1->GetMaxTimeInShop() * _checkDisRate;
        l        = l < _socialDistance ? _socialDistance : l;
        eff_dist = distp12.Norm() <= l ? 0 : eff_dist;
    } else if(InRightRoom(room, {"waiting"})) {
        eff_dist = distp12.Norm() - 2 * 0.4;
    }
    double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
    if(condition1 < 0) {
        return my_pair(FLT_MAX, -1);
    }
    //When condition1==0,should no influence on velocity
    bool parallel = almostEqual(condition1, 0.0, J_EPS);
    if(parallel) {
        return my_pair(FLT_MAX, -1);
    }
    //Judge conllision
    //Obtain parameters
    double a1      = ped1->GetLargerAxis();
    Point v        = ped1->GetV();
    double b1      = ped1->GetSmallerAxis();
    b1             = ped1->GetEllipse().GetBmin();
    double a2      = ped2->GetLargerAxis();
    double b2      = ped2->GetSmallerAxis();
    double x2      = ped2->GetPos()._x;
    double y2      = ped2->GetPos()._y;
    double cosphi1 = ei.Normalized()._x;
    double sinphi1 = ei.Normalized()._y;
    double cosphi2 = ped2->GetEllipse().GetCosPhi();
    double sinphi2 = ped2->GetEllipse().GetSinPhi();
    //Judge the position of the center of ped2
    double d1 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) + b1;
    double d2 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) - b1;
    if(d1 * d2 <= 0) {
        //if the center between two lines, collision
        return my_pair(eff_dist, ped2->GetID());
    }
    //If the center not between two lines, Judge if ped2 contact with two lines
    Point D;
    D._x = cosphi1;
    D._y = sinphi1;
    Point De;
    De._x = cosphi1 * cosphi2 + sinphi1 * sinphi2;
    De._y = sinphi1 * cosphi2 - sinphi2 * cosphi1;
    Point Ne;
    Ne._x = -De._y;
    Ne._y = De._x;
    Point A1;
    A1._x = x1 + b1 * sinphi1;
    A1._y = y1 - b1 * cosphi1;
    Point A2;
    A2._x = x1 - b1 * sinphi1;
    A2._y = y1 + b1 * cosphi1;
    //Transfer A1 and A2 to ped2 coordinate
    //Point A1e = ((A1._x - x2)*cosphi2 + (A1._y - y2)*sinphi2, (A1._y - y2)*cosphi2 - (A1._x - x2)*sinphi2);
    Point A1e = A1.TransformToEllipseCoordinates(ped2->GetPos(), cosphi2, sinphi2);
    //Point A2e = ((A2._x - x2)*cosphi2 + (A2._y - y2)*sinphi2, (A2._y - y2)*cosphi2 - (A2._x - x2)*sinphi2);
    Point A2e = A2.TransformToEllipseCoordinates(ped2->GetPos(), cosphi2, sinphi2);
    // Judge if the direction of De is right (ellipse2 coordinate)
    double J1 = Ne.ScalarProduct(A1e);
    double J2 = Ne.ScalarProduct(A2e);
    Point De1 = (J1 >= 0) ? De * (-1, -1) : De;
    Point De2 = (J2 >= 0) ? De * (-1, -1) : De;
    //Calculate point R (ellipse2 coordinate)
    Point Ne1;
    Ne1._x = -De1._y;
    Ne1._y = De1._x;
    Point Ne2;
    Ne2._x = -De2._y;
    Ne2._y = De2._x;
    Point Te1;
    Te1._x = De1._y / b2;
    Te1._y = -De1._x / a2;
    Point Te2;
    Te2._x     = De2._y / b2;
    Te2._y     = -De2._x / a2;
    Point Te1n = Te1.Normalized();
    Point Te2n = Te2.Normalized();
    Point Re1;
    Re1._x = a2 * Te1n._x;
    Re1._y = b2 * Te1n._y;
    Point Re2;
    Re2._x = a2 * Te2n._x;
    Re2._y = b2 * Te2n._y;
    //Transfer R to global coordinate
    /*
	Point R1 = (Re1._x*cosphi2 - Re1._y*sinphi2 + x2, Re1._y*cosphi2 + Re1._x*sinphi2 + y2);
	Point R1 = Re1.TransformToCartesianCoordinates(ped2->GetPos(), cosphi2, sinphi2);
	Point R2 = (Re2._x*cosphi2 - Re2._y*sinphi2 + x2, Re2._y*cosphi2 + Re2._x*sinphi2 + y2);
	Point R2 = Re2.TransformToCartesianCoordinates(ped2->GetPos(), cosphi2, sinphi2);
	*/
    //Calculate distance between point R and line
    double Dis1 = Ne1.ScalarProduct(Re1) - Ne1.ScalarProduct(A1e);
    double Dis2 = Ne2.ScalarProduct(Re2) - Ne2.ScalarProduct(A2e);

    //Judge if the line contact with ellipse2
    if(Dis1 >= 0 && Dis2 >= 0)
        return my_pair(FLT_MAX, -1);
    else {
        return my_pair(eff_dist, ped2->GetID());
    }
}
