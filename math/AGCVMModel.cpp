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
int ifdrill = 1;
int ifcorrection = 1;
int vertialnoinf = 1;

AGCVMModel::AGCVMModel(std::shared_ptr<DirectionStrategy> dir, double aped, double Dped,
	double awall, double Dwall, double Ts, double Td, int GCVM,
	int Parallel, double waitingTime, double lb, double rb, double ub, double db, double co, 
	int Anticipation, int Cooperation, int AttracForce, double AntiT, double CoopT, int Push, double CoreSize)
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
	_GCVMUsing = GCVM; // If using GCVM

	// Clogging Parameter
	_Parallel = Parallel;
	_WaitingTime = waitingTime;

	// Boundary Case
	_left_boundary = lb;
	_right_boundary = rb;
	_up_boundary = ub;
	_down_boundary = db;
	_cutoff = co;

	_Anticipation = Anticipation;
	_Cooperation = Cooperation;
	_AttracForce = AttracForce;
	_AntiT = AntiT;
	_CoopT = CoopT;
	_Pushing = Push;
	_CoreSize = CoreSize;
}


AGCVMModel::~AGCVMModel()
{

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
				"ERROR:\tGCVMModel::Init() cannot initialise route. ped %d is deleted in Room %d %d.\n", ped->GetID(), ped->GetRoomID(), ped->GetSubRoomID());
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

	//Calculate direction for each pedestrians
	for (int p = start; p <= end; ++p)
	{
		Pedestrian* ped1 = allPeds[p];
		Point p1 = ped1->GetPos();
		Room* room = building->GetRoom(ped1->GetRoomID());
		SubRoom* subroom = room->GetSubRoom(ped1->GetSubRoomID());
		// v:Effect from neighbours, we consider two kinds effect here
		Point IniDirection = e0(ped1, room);//desired moving direction, direction3, and using core_size here.
		vector<Pedestrian*> neighbours;
		building->GetGrid()->GetNeighbourhood(ped1, neighbours);
		bool drill = Drill(ped1, neighbours, building, subroom, IniDirection,periodic);
		if (ShowInfo == 1)
		{
			printf("\nTime=%f, ID1=%d (%f, %f), e0=(%f, %f), drill=%d",
				current, ped1->GetID(), ped1->GetPos()._x, ped1->GetPos()._y,
				IniDirection._x, IniDirection._y, drill);
		}
		int size = (int)neighbours.size();
		Point repPed = Point(0, 0);
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
				Point ForcePush = ForceRepPed(ped1, ped2, IniDirection, periodic, true);//new method
				repPedPush += ForcePush;
			}
			if (1)// || subroom->IsDirectlyConnectedWith(subroom2))
			{
				Point Force = ForceRepPed(ped1, ped2, IniDirection, periodic,false);//new method
				repPed += Force;
			}
		} //for i
		
		Point repWall = ForceRepRoom(ped1, subroom, IniDirection);
		// the influence from next subroom should be considered
		/*
		for (const auto & subr : subroom->GetNeighbors())
		{
			repWall = repWall + ForceRepRoom(ped1, subr, IniDirection);
		}
		*/
		Point direction;
		Point a_direction=ped1->GetMoveDirection();
		Point d_direction = IniDirection + repPed + repWall+ repPedPush;
		if (ifdrill == 1)
		{
			if (drill == true)
			{
				d_direction = IniDirection.Normalized();
			}
		}
		if (ifcorrection == 1)
		{
			d_direction = CorrectD(ped1, d_direction, subroom);
		}
		Point AccTu = Point(0, 0);
		double angle_tau = GetTd();
		d_direction = d_direction.Normalized();
		AccTu = (d_direction.Normalized()-ped1->GetV()) / angle_tau;
		direction = ped1->GetV() + AccTu * deltaT;
		direction = direction.Normalized();
	
		//direction = IniDirection.Normalized();
		if (GetGCVMU() == 0)
		{
			direction = d_direction.Normalized();//original method
		}
		if (ped1->GetID()==-1)
		{
			Point old_d = (IniDirection + repPed + repWall + repPedPush).Normalized();
			printf("\nTime=%f, ID1=%d (%f, %f), a_direction=(%f, %f), new_d=(%f, %f), old_d=(%f, %f)\n",
				current, ped1->GetID(), ped1->GetPos()._x, ped1->GetPos()._y,
				a_direction._x, a_direction._y, d_direction._x, d_direction._y, old_d._x, old_d._y);
		}
		result_dir.push_back(direction);	
	}

	//update direction of each pedestrian
	for (int p = start; p <= end; ++p)
	{
		Pedestrian* ped = allPeds[p];
		Point direction = result_dir[p];
		ped->SetMoveDirection(direction);
	}

	// Calculate speed and minimal TTC
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

		// Saving all the information of ttcs
		vector< my_pair> ttcs = vector<my_pair>();
		ttcs.reserve(size);

		//Calculating spacing in front -------------------------------------------------------------------------------------------------------
		for (int i = 0; i < size; i++)
		{
			Pedestrian* ped2 = neighbours[i];
			SubRoom* subroom2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());
			if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom->IsDirectlyConnectedWith(subroom2)||1)
			{
				if (GetPushing() == 1)
				{
					spacings.push_back(GetSpacing(ped1, ped2, periodic, true));
				}
				else
				{
					spacings.push_back(GetSpacing(ped1, ped2, periodic, false));
				}
				ttcs.push_back(JudgeCollision(ped1, ped2,periodic));
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

		// Calculate min ttc and save
		std::sort(ttcs.begin(), ttcs.end(), sort_pred_agcvm());
		double mttc = ttcs.size() == 0 ? FLT_MAX : ttcs[0].first;
		if (mttc < GetCoopT())
		{
			ped1->SetMTTCP(ttcs[0].second);
		}
		else
		{
			ped1->SetMTTCP(-1);
		}
		/*
		printf("\nTime=%f, ID1=%d, ID2=%d, ttc=%f\n",
			current, ped1->GetID(), ped1->GetMTTCP(), ttcs[0].first);
		*/
		// Optimal speed function
		Point speed;
		Point ei = ped1->GetMoveDirection();
		speed = ei *OptimalSpeed(ped1, spacing);
		normal_acc.push_back(speed);
		if (ped1->GetID()==-1)
		{
			printf("\nTime=%f, ID1=%d (%f, %f), spacing=%f, speed=(%f, %f)\n",
				current, ped1->GetID(), ped1->GetPos()._x, ped1->GetPos()._y,
				spacing, speed._x, speed._y);
		}
	} // for p
	
	//Cooperation
	if (GetCooperation() == 1)
	{
		for (int p = start; p <= end; ++p)
		{
			Pedestrian* ped1 = allPeds[p];
			vector<Pedestrian*> neighbours;
			building->GetGrid()->GetNeighbourhood(ped1, neighbours);
			int size = (int)neighbours.size();
			if (ped1->GetMTTCP() == -1)
			{
				continue;
			}
			for (int i = 0; i < size; i++)
			{
				Pedestrian* ped2 = neighbours[i];
				if (ped2->GetID() == ped1->GetMTTCP())
				{
					double coop1 = ped1->GetCooperation();
					double coop2 = ped2->GetCooperation();
					double s1 = GetSpacing(ped1, ped2, periodic, false).first;
					double s2 = GetSpacing(ped2, ped1, periodic, false).first;
					if ((s1 < 100 && s2>100)|| (coop1 > coop2&&s1 > 100 && s2 > 100)|| (coop1 > coop2&&s1 < 100 && s2 < 100))
					{
						/*
						printf("\nTime=%f, ID1=%d (%f, %f) (coop=%f,v0=%f), ID2=%d(%f, %f) (coop=%f,v0=%f)\n",
							current, ped1->GetID(), ped1->GetPos()._x, ped1->GetPos()._y, ped1->GetCooperation(), ped1->GetV0Norm(),
							ped2->GetID(), ped2->GetPos()._x, ped2->GetPos()._y, ped2->GetCooperation(), ped2->GetV0Norm());
						*/	
						normal_acc[p] = Point(0, 0);
					}
					break;
				}
			}

		}
	}
	
	//Update everything
	for (int p = start; p <= end; ++p)
	{
		Pedestrian* ped = allPeds[p];
		Point v_neu = normal_acc[p];
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

Point AGCVMModel::e0(Pedestrian* ped, Room* room) const
{
	const Point target = _direction->GetTarget(room, ped); // target is where the ped wants to be after the next timestep
	Point desired_direction;
	const Point pos = ped->GetPos();
	double dist = (target - pos).Norm();
	if (dist > J_EPS_GOAL) // Adjust in the future
	{
		Point NewE0;
		NewE0=(target - pos).Normalized();
		ped->SetLastE0(NewE0);
	}
	desired_direction = ped->GetLastE0();
	return desired_direction;
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

my_pair AGCVMModel::GetSpacing(Pedestrian* ped1, Pedestrian* ped2, int periodic, bool collision) const
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
		printf("ERROR: \tin AGCVMModel::GetSpacing() ep12 can not be calculated!!!\n");
		exit(EXIT_FAILURE);
	}

	//Judge conllision
	double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
	double l = collision == true ? 2 * GetCoreSize() : (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());
	double condition2 = ei.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
	condition2 = (condition2 > 0) ? condition2 : -condition2; // abs
	if ((condition1 > 0) && (condition2 < l / Distance))
	{
		if (ped1->GetID()==-1)
		{
			printf("\nID1=%d (%f, %f), ID2=%d (%f, %f), Dis=%f, condition1=%f, ei=(%f,%f), ep12=(%f,%f)",
				ped1->GetID(), ped1->GetPos()._x, ped1->GetPos()._y,
				ped2->GetID(), ped2->GetPos()._x, ped2->GetPos()._y,
				Distance - l, condition1, ei._x, ei._y, ep12._x, ep12._y);
		}
		return  my_pair((Distance - l), ped2->GetID());
	}
	else
	{
		return  my_pair(FLT_MAX, -1);
	}
}

Point AGCVMModel::ForceRepPed(Pedestrian* ped1, Pedestrian* ped2, Point e0, int periodic, bool push) const
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
		printf("ERROR: \tin AGCVMModel::forcePedPed() ep12 can not be calculated!!!\n");
		printf("ped1 %d  ped2 %d\n", ped1->GetID(), ped2->GetID());
		printf("ped1 at (%f, %f), ped2 at (%f, %f)\n", p1._x, p1._y, p2._x, p2._y);
		exit(EXIT_FAILURE);
	}
	Point ei = ped1->GetMoveDirection();
	Point ei2 = ped2->GetMoveDirection();
	double dist = Distance - (ped1->GetEllipse().GetBmax() + ped2->GetEllipse().GetBmax());// Distance between boundaries
	if (push == true)
	{
		if (dist <= J_EPS)
		{
			double R_ij = _aPed * exp((-dist) / _DPed);
			F_rep = ep12 * (-R_ij);
			if (ShowInfo == -1)
			{
				printf("\nID1=%d (%f, %f), ID2=%d (%f, %f), Dis=%f, PushForce2->1=(%f, %f)",
					ped1->GetID(), p1._x, p1._y,
					ped2->GetID(), p2._x, p2._y,
					dist, F_rep._x, F_rep._y);
			}
		}
		return F_rep;
	}

	//Vision area
	double condition1 = e0.ScalarProduct(ep12);
	double condition2 = ei.ScalarProduct(ep12);

	//Anticipation
	double S_Gap = 0;
	double Dis_Gap = 0;
	double multi_e0 = ped1->GetV0().ScalarProduct(ped2->GetV0());
	double beta = (1 - multi_e0)/2;
	if (GetAnticipation() == 1)
	{
		double t_anti = GetAntiT();//Anticipation time
		//New S_Gap: Using desired speed instead of real speed
		S_Gap = (ei.ScalarProduct(ep12)*ped1->GetV0Norm() - ei2.ScalarProduct(ep12)*ped2->GetV0Norm());
		Dis_Gap=S_Gap * t_anti* beta;
	}
	
	double R_ij;
	if (GetGCVMU() == 0)
	{
		R_ij = _aPed * exp((-dist) / _DPed);
		F_rep = ep12 * (-R_ij);//ep12: from 1(i) to 2(j)
		return F_rep;
	}
	
	if (condition1>=-J_EPS||condition2>=-J_EPS)//rule:pedestrian's direction only influenced by pedestrian in vision area
	{
		Point infd = GetInfDirection(e0, ep12);
		double condition3 = e0.ScalarProduct(ei2);// ped2 move in the same direction of ped1's e0;		
		if ((GetAttracForce() == 1) && condition1>0&& condition2>0 && condition3>0 && S_Gap < 0 && dist>J_EPS)
		{
			double R_dist = dist + Dis_Gap;
			R_dist = R_dist < 0 ? 0 : R_dist;
			R_ij = _aPed * exp((-R_dist) / _DPed);
			infd = infd * -1;
			F_rep = infd * R_ij;
		}
		else
		{
			double R_dist = dist - Dis_Gap;
			R_dist = R_dist < 0 ? 0 : R_dist;
			R_ij = _aPed * exp((-R_dist ) / _DPed);
			F_rep = infd * R_ij;
		}
	}
	if (ped1->GetID()==-1)
	{
		printf("\nID1=%d (%f, %f), ID2=%d (%f, %f), Dis=%f, Force2->1=(%f, %f),condition1=%f, condition2=%f",
			ped1->GetID(), p1._x, p1._y,
			ped2->GetID(), p2._x, p2._y,
			dist, F_rep._x, F_rep._y,
			condition1,condition2);
	}
	if (vertialnoinf == 1)
	{
		F_rep = F_rep * abs(multi_e0);
	}
	return F_rep;
}

Point AGCVMModel::ForceRepRoom(Pedestrian* ped, SubRoom* subroom, Point e0) const
{
	Point f(0., 0.);
	const Point& centroid = subroom->GetCentroid();
	bool inside = subroom->IsInSubRoom(centroid);
	//first the walls
	for (const auto & wall : subroom->GetAllWalls())
	{
		f += ForceRepWall(ped, wall, centroid, inside, e0);
	}

	//then the obstacles
	for (const auto & obst : subroom->GetAllObstacles())
	{
		if (obst->Contains(ped->GetPos()))
		{
			Log->Write("ERROR:\t Agent [%d] is trapped in obstacle in room/subroom [%d/%d]", ped->GetID(), subroom->GetRoomID(), subroom->GetSubRoomID());
			exit(EXIT_FAILURE);
		}
		else
			for (const auto & wall : obst->GetAllWalls())
			{
				f += ForceRepWall(ped, wall, centroid, inside, e0);
			}
	}

	// and finally the closed doors
	for (const auto & goal : subroom->GetAllTransitions())
	{
		if (!goal->IsOpen())
		{
			f += ForceRepWall(ped, *(static_cast<Line*>(goal)), centroid, inside, e0);
		}
		/*
		//door is open, but it's not my door (has influence)
		int uid1 = goal->GetUniqueID();
		int uid2 = ped->GetExitIndex();
		if((uid1 != uid2) && (goal->IsOpen()==true ))
		{
			f +=  ForceRepWall(ped,*(static_cast<Line*>(goal)), centroid, inside, e0);
		}
		*/
	}
	return f;
}

Point AGCVMModel::ForceRepWall(Pedestrian* ped, const Line& w, const Point& centroid, bool inside, Point e0) const
{
	Point F_wrep = Point(0.0, 0.0);
	Point pt = w.ShortestPoint(ped->GetPos());
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

	//rule: wall in behind has no influence 
	Point pt1 = w.GetPoint1();
	Point pt2 = w.GetPoint2();
	Point dist_pt1 = pt1 - ped->GetPos();
	Point dist_pt2 = pt2 - ped->GetPos();
	Point e_iw1 = dist_pt1.Normalized();
	Point e_iw2 = dist_pt2.Normalized();
	Point ei = ped->GetMoveDirection();
	//vision area----------------------------------
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
	double R_iw = _aWall * exp((-effdis) / _DWall);
	if (GetGCVMU() == 1)
	{
		F_wrep = inf_direction * R_iw;//new method
	}
	else
	{
		F_wrep = e_iw * R_iw;//original method
	}
	if (ped->GetID()==-1)
	{
		printf("\nID1=%d (%f, %f), Dis_wall=%f, Force_wall=(%f, %f)",
			ped->GetID(), ped->GetPos()._x, ped->GetPos()._y,
			effdis, inf_direction._x, inf_direction._y);
	}
	return F_wrep;
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
	double b = ped->GetCore()-J_EPS;
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
	if (ped->GetID() == -1)
	{
		printf("\nID1=%d (%f, %f), Dis=%f, spacing=%f",
			ped->GetID(), ped->GetPos()._x, ped->GetPos()._y,spacing);
	}
	return spacing;
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

std::shared_ptr<DirectionStrategy> AGCVMModel::GetDirection() const
{
	return _direction;
}

double AGCVMModel::GetaPed() const
{
	return _aPed;
}

double AGCVMModel::GetDPed() const
{
	return _DPed;
}

double AGCVMModel::GetaWall() const
{
	return _aWall;
}

double AGCVMModel::GetDWall() const
{
	return _DWall;
}

double AGCVMModel::GetTs() const
{
	return _Ts;
}

double AGCVMModel::GetTd() const
{
	return _Td;
}

int AGCVMModel::GetGCVMU() const
{
	return _GCVMUsing;
}

int AGCVMModel::GetUpdate() const
{
	return _Parallel;
}

double AGCVMModel::GetWaitingTime() const
{
	return _WaitingTime;
}

double AGCVMModel::GetLeftBoundary() const
{
	return _left_boundary;
}

double AGCVMModel::GetRightBoundary() const
{
	return _right_boundary;
}

double AGCVMModel::GetUpBoundary() const
{
	return _up_boundary;
}

double AGCVMModel::GetDownBoundary() const
{
	return _down_boundary;
}

double AGCVMModel::GetCutoff() const
{
	return _cutoff;
}

void AGCVMModel::UpdatePed(Pedestrian* ped, Point speed, Point direction, double deltaT, int periodic)
{
	Point pos_neu = ped->GetPos() + speed * deltaT;
	Point e0 = ped->GetLastE0();
	Point MD = direction;
	Point MS = speed;
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

int AGCVMModel::GetAnticipation() const
{
	return _Anticipation;
}

int AGCVMModel::GetCooperation() const
{
	return _Cooperation;
}

int AGCVMModel::GetAttracForce() const
{
	return _AttracForce;
}

double AGCVMModel::GetAntiT() const
{
	return _AntiT;
}

double AGCVMModel::GetCoopT() const
{
	return _CoopT;
}

int AGCVMModel::GetPushing() const
{
	return _Pushing;
}

double AGCVMModel::GetCoreSize() const
{
	return _CoreSize;
}
my_pair AGCVMModel::JudgeCollision(Pedestrian* ped1, Pedestrian* ped2, int periodic) const
{
	double ttc = FLT_MAX;
	double At = GetCoopT();
	//At = 1;
	Point p1 = ped1->GetPos();
	Point p2 = periodic? GetPosPeriodic(ped1, ped2):ped2->GetPos();
	Point  e1 = ped1->GetMoveDirection();
	e1 = e1 * ped1->GetV0Norm();
	//e1 = ped1->GetLastE0()* ped2->GetV0Norm();
	Point  e2 = ped2->GetMoveDirection();
	e2 = e2 * ped2->GetV0Norm();
	//e2=ped2->GetLastE0()* ped2->GetV0Norm();
	Point distp12 = p2 - p1; //ped1 ---> ped2
	double distance = distp12.Norm();
	if (e1.ScalarProduct(distp12) < 0)
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
		/*
		if (ped1->GetID() == 2 && ped2->GetID() == 20)
		{
			printf("\nID1=%d (v0=%f), ID2=%d (v0=%f), delta=%f, t1=%f, t2=%f\n",
				ped1->GetID(), ped1->GetV0Norm(), ped2->GetID(), ped2->GetV0Norm(), delta, t1,t2);
		}
		*/
	}
	return my_pair(ttc,ped2->GetID());
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
		double x2_periodic = x2 + xR- xL;
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

	for (const auto & wall : subroom->GetAllWalls())
	{
		if (DrillWall(ped1, e0, wall) == false)
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
			if (DrillWall(ped1, e0, wall) == false)
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
			if (DrillWall(ped1, e0, *(static_cast<Line*>(goal))) == false)
			{
				drill = false;
				return drill;
			}
		}
	}
	return drill;
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
	d_direction= d_direction.Normalized();
	Point direction = d_direction;
	for (const auto & wall : subroom->GetAllWalls())
	{
		direction = CorrectDWall(ped, d_direction, wall);
		if ((d_direction- direction).Norm()>J_EPS)
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
	if (Distance-radius <= J_EPS && costheta>0)
	{
		direction = d_direction - e_iw * costheta;
	}
	return direction.Normalized();
}
