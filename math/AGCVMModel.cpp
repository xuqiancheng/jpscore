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

AGCVMModel::AGCVMModel(std::shared_ptr<DirectionStrategy> dir, double aped, double Dped,
	double awall, double Dwall, double Ts, double Td, int GCVM,
	int Parallel, double waitingTime, double lb, double rb, double ub, double db, double co, 
	int Anticipation, int ContactRep, int AttracForce, double AntiT)
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
	_ContactRep = ContactRep;
	_AttracForce = AttracForce;
	_AntiT = AntiT;
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
	vector< Point > coop_acc = vector<Point >();
	coop_acc.reserve(nSize);
	vector< Point > back_acc = vector<Point >();
	back_acc.reserve(nSize);
	vector< Point > defect_acc = vector<Point >();
	defect_acc.reserve(nSize);
	vector< Point > result_dir = vector<Point >();
	result_dir.reserve(nSize);

	int start = 0;
	int end = nSize - 1;

	//Calculate direction
	for (int p = start; p <= end; ++p)
	{
		Pedestrian* ped1 = allPeds[p];
		Room* room = building->GetRoom(ped1->GetRoomID());
		SubRoom* subroom = room->GetSubRoom(ped1->GetSubRoomID());

		Point repPed = Point(0, 0);
		vector<Pedestrian*> neighbours;
		building->GetGrid()->GetNeighbourhood(ped1, neighbours);
		Point inid_direction = e0(ped1, room);

		//Using a new method calculate the influence of pedestrian (the value of influence id decided by distance and the direction is vertical with desired direction) 
		Point inf_direction;// left side of pedestrian
		inf_direction._x = -inid_direction._y;
		inf_direction._y = inid_direction._x;
		inf_direction = inf_direction.Normalized();
		//--------------------------------------------------------------------------------------------------------------------------------------------------------------

		int size = (int)neighbours.size();
		for (int i = 0; i < size; i++)
		{
			Pedestrian* ped2 = neighbours[i];
			//if they are in the same subroom
			Point p1 = ped1->GetPos();
			Point p2 = ped2->GetPos();
			Point ep12 = p2 - p1;

			//Deciding the influence direction--------------------------------------------------------------------
			double result1 = inid_direction.CrossProduct(ep12);
			Point zero = Point(0, 0);
			if (bool equal = almostEqual(result1, 0, 0.00001))//if neighbour is in front or behind pedestrian
			{
				int random = rand() % 1000;//choose one direciton bu random
				if (random < 500)//when random is larger than 50, influence's direction is right, otherwise is left
				{
					inf_direction = zero - inf_direction;
				}
			}
			else
			{
				double result2 = inid_direction.CrossProduct(inf_direction);
				if (result1*result2 > 0)//is ep12 and inf_direction in the same side of inid_direction, change the direction of inf_direction
				{
					inf_direction = zero - inf_direction;
				}
			}

			//Check if two pedestrians can see each other
			SubRoom* sb2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());
			vector<SubRoom*> emptyVector;
			emptyVector.push_back(subroom);
			emptyVector.push_back(sb2);
			bool isVisible = building->IsVisible(p1, p2, emptyVector, false);
			if (!isVisible)
			{
				continue;
			}
			if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom->IsDirectlyConnectedWith(sb2))
			{
				repPed += ForceRepPed(ped1, ped2, inid_direction, inf_direction, periodic);//new method
			}
		} //for i

		// Calculating influence of walls-------------------------------------------------------------------------
		Point repWall = ForceRepRoom(ped1, subroom, inid_direction);

		//Caluculating desired direcition----------------------------------------------------------------------------------------------
		Point d_direction;
		d_direction = inid_direction + repPed + repWall;//new method

		//Calculating the actual direction of pedestrian at next timestep------------------------------------------------------------------
		Point a_direction;
		a_direction._x = ped1->GetEllipse().GetCosPhi();
		a_direction._y = ped1->GetEllipse().GetSinPhi();
		double angle_tau = GetTd();
		Point angle_v = (d_direction.Normalized() - a_direction) / angle_tau;
		Point direction = a_direction + angle_v * deltaT;

		
		/*
		//We dont use it now
		if (GetContactRep() == 1)// when using contact forece, pedestrian can move backward
		{
			//When the angle between actual moving direction and initial desired moving direction is bigger than 90 degree. turning to the desired moving direction directly.
			bool c_backward = (d_direction.ScalarProduct(inid_direction) < 0);//going move backward
			bool c_forward = (inid_direction.ScalarProduct(d_direction) > 0);// going move forward
			bool on_backward = (inid_direction.ScalarProduct(a_direction) < 0);// move backward now
			bool on_forward = (inid_direction.ScalarProduct(a_direction) > 0);// move forward now
			double condition_final = (inid_direction.ScalarProduct(d_direction)*inid_direction.ScalarProduct(a_direction));
			double condition_backward = d_direction.ScalarProduct(a_direction);
			if (condition_final < 0 && condition_backward < 0)
			{
				direction = d_direction.Normalized();
			}
		}
		*/
		

		if (GetGCVMU() == 0)
		{
			direction = d_direction.Normalized();//original method
		}
		direction = direction.Normalized();
		result_dir.push_back(direction);
	}

	//update direction of each pedestrian
	for (int p = start; p <= end; ++p)
	{
		Pedestrian* ped = allPeds[p];
		Point dir_neu = result_dir[p];
		double normdir = dir_neu.Norm();
		double cosPhi = dir_neu._x / normdir;
		double sinPhi = dir_neu._y / normdir;
		JEllipse e = ped->GetEllipse();
		e.SetCosPhi(cosPhi);
		e.SetSinPhi(sinPhi);
		ped->SetEllipse(e);
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
		vector< my_pair > spacings_back = vector<my_pair >();
		spacings_back.reserve(size);
		vector< my_pair > spacings_defect = vector<my_pair >();
		spacings_defect.reserve(size);
		Point e1 = result_dir[p];

		//Calculating spacing in front -------------------------------------------------------------------------------------------------------
		for (int i = 0; i < size; i++)
		{
			Pedestrian* ped2 = neighbours[i];
			SubRoom* sb2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());
			if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom->IsDirectlyConnectedWith(sb2))
			{
				spacings.push_back(GetSpacing(ped1, ped2, e1, periodic));
				spacings_defect.push_back(GetSpacing_collision(ped1, ped2, e1, periodic));
				spacings_back.push_back(GetSpacing_collision(ped1, ped2, e1*-1, periodic));
			}
		}//for i

		// Calculate min spacing
		std::sort(spacings.begin(), spacings.end(), sort_pred_agcvm());
		double spacing = spacings.size() == 0 ? 100 : spacings[0].first;
		double spacing_wall = GetSpacingRoom(ped1, subroom, e1);
		spacing = spacing < spacing_wall ? spacing : spacing_wall;

		std::sort(spacings_back.begin(), spacings_back.end(), sort_pred_agcvm());
		double spacing_back = spacings_back.size() == 0 ? 100 : spacings_back[0].first;
		double spacing_back_wall = GetSpacingRoom(ped1, subroom, e1*-1);
		spacing_back = spacing_back < spacing_back_wall ? spacing_back : spacing_back_wall;

		std::sort(spacings_defect.begin(), spacings_defect.end(), sort_pred_agcvm());
		double spacing_defect = spacings_defect.size() == 0 ? 100 : spacings_defect[0].first;
		spacing_defect = spacing_defect < spacing_wall ? spacing_defect : spacing_wall;

		// Optimal speed function
		Point speed;
		speed = e1.Normalized() *OptimalSpeed(ped1, spacing);
		Point back_speed;
		back_speed = e1.Normalized() *OptimalSpeed(ped1, spacing_back)*-1;
		Point defect_speed;
		defect_speed= e1.Normalized() *OptimalSpeed(ped1, spacing_defect);

		normal_acc.push_back(speed);
		back_acc.push_back(back_speed);
		defect_acc.push_back(defect_speed);
	} // for p
	

	//update speed
	for (int p = start; p <= end; ++p)
	{
		Pedestrian* ped = allPeds[p];
		Point v_neu = normal_acc[p];
		ped->SetV(v_neu);
		//if (ped->GetID() == 13 || ped->GetID() == 14)
		//printf("\nTest: t=%f, ID=%d, v_neu=(%f,%f), ped_v=(%f,%f).",current,ped->GetID(), v_neu._x, v_neu._y, ped->GetV()._x, ped->GetV()._y);
	}

	// Cooperation

		for (int p = start; p <= end; ++p)
		{
			Pedestrian* ped1 = allPeds[p];
			Room* room = building->GetRoom(ped1->GetRoomID());
			SubRoom* subroom = room->GetSubRoom(ped1->GetSubRoomID());

			vector<Pedestrian*> neighbours;
			building->GetGrid()->GetNeighbourhood(ped1, neighbours);
			int size = (int)neighbours.size();

			Point e1 = result_dir[p];

			int back_move = 0;
			int try_coop = 0;

			for (int i = 0; i < size; i++)
			{
				Pedestrian* ped2 = neighbours[i];
				SubRoom* sb2 = building->GetRoom(ped2->GetRoomID())->GetSubRoom(ped2->GetSubRoomID());

				Point p1 = ped1->GetPos();
				Point p2 = ped2->GetPos();
				Point ep12 = p2 - p1;
				ep12 = ep12.Normalized();

				double condition1 = e1.ScalarProduct(ep12);
				double cosphi = ped2->GetEllipse().GetCosPhi();
				double sinphi = ped2->GetEllipse().GetSinPhi();
				Point  e2 = Point(cosphi, sinphi);
				double condition2 = e2.ScalarProduct(ep12);

				if (ped1->GetUniqueRoomID() == ped2->GetUniqueRoomID() || subroom->IsDirectlyConnectedWith(sb2))
				{
					if (condition1 >= 0 && condition2 <= 0)
					{
						int collision = JudgeCollision(ped1, ped2);
						if (collision == 1)
						{
							try_coop = 1;

							double coop1 = ped1->GetCooperation();
							double coop2 = ped2->GetCooperation();
							Point d1 = ped1->GetV();
							Point d2 = ped2->GetV();
							double v1 = d1.Norm();
							double v2 = d2.Norm();

							if (v1 > 0.0001 && v2 <= 0.0001)
							{
								//if (ped1->GetID() == 8 || ped1->GetID() == 14)
								//printf("\nDefc,t=%f: (ID1=%d, v1=[%f,%f], cooperation1=%f, c1=%f), (ID2=%d, v2=[%f,%f], cooperation2=%f, c2=%f)", current, ped1->GetID(), d1._x, d1._y, coop1, c1, ped2->GetID(), d2._x, d2._y, coop2, c2);
							}
							else if (v1 <= 0.0001 && v2 > 0.0001)
							{
								//if (ped1->GetID()==14|| ped1->GetID()==8)
								//printf("\nCoop,t=%f: (ID1=%d, v1=[%f,%f], cooperation1=%f, c1=%f), (ID2=%d, v2=[%f,%f], cooperation2=%f, c2=%f)", current, ped1->GetID(), d1._x, d1._y, coop1, c1, ped2->GetID(), d2._x, d2._y, coop2, c2);
								back_move = 1;
							}
							else
							{
								if (coop1 > coop2)
								{
									//if (ped1->GetID()==14|| ped1->GetID()==8)
									//printf("\nCoop,t=%f: (ID1=%d, v1=[%f,%f], cooperation1=%f, c1=%f), (ID2=%d, v2=[%f,%f], cooperation2=%f, c2=%f)", current, ped1->GetID(), d1._x, d1._y, coop1, c1, ped2->GetID(), d2._x, d2._y, coop2, c2);
									back_move = 1;
								}
								else
								{
									//if (ped1->GetID() == 8 || ped1->GetID() == 14)
									//printf("\nDefc,t=%f: (ID1=%d, v1=[%f,%f], cooperation1=%f, c1=%f), (ID2=%d, v2=[%f,%f], cooperation2=%f, c2=%f)", current, ped1->GetID(), d1._x, d1._y, coop1, c1, ped2->GetID(), d2._x, d2._y, coop2, c2);
								}
							}
						}
					}
				}
			}//for i
			// Optimal speed function

			Point speed;
			if (0)
			{
				if (try_coop == 0)
				{
					speed = normal_acc[p];
				}
				else
				{
					if (back_move == 0)
					{
						speed = defect_acc[p];
					}
					else
					{
						speed = back_acc[p];
					}
				}
			}
			else
			{
				speed = normal_acc[p];
			}
			coop_acc.push_back(speed);
		} // for p
	
	//Update everything
	for (int p = start; p <= end; ++p)
	{
		Pedestrian* ped = allPeds[p];
		Point v_neu = coop_acc[p];
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
	double dist = ped->GetExitLine()->DistTo(pos);
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

my_pair AGCVMModel::GetSpacing(Pedestrian* ped1, Pedestrian* ped2, Point ei, int periodic) const
{
	double x1 = ped1->GetPos()._x;
	double y1 = ped1->GetPos()._y;
	double x2_real = ped2->GetPos()._x;
	double y2_real = ped2->GetPos()._y;
	Point ped2_current = ped2->GetPos();

	if (periodic) 
	{
		double xLeft_gcvm = GetLeftBoundary();
		double xRight_gcvm = GetRightBoundary();
		double yUp_gcvm = GetUpBoundary();
		double yDown_gcvm = GetDownBoundary();
		double cutoff_gcvm = GetCutoff();
		if ((xRight_gcvm - x1) + (x2_real - xLeft_gcvm) <= cutoff_gcvm) {
			double x2_periodic = x2_real + xRight_gcvm - xLeft_gcvm;
			ped2->SetPos(Point(x2_periodic, y2_real));
		}
		if ((x1 - xLeft_gcvm) + (xRight_gcvm - x2_real) <= cutoff_gcvm) {
			double x2_periodic = xLeft_gcvm - xRight_gcvm + x2_real;
			ped2->SetPos(Point(x2_periodic, y2_real));
		}
		if ((y1 - yDown_gcvm) + (yUp_gcvm - y2_real) <= cutoff_gcvm) {
			double y2_periodic = yDown_gcvm - yUp_gcvm + y2_real;
			ped2->SetPos(Point(x2_real, y2_periodic));
		}
		if ((y2_real - yDown_gcvm) + (yUp_gcvm - y1) <= cutoff_gcvm) {
			double y2_periodic = yUp_gcvm + y2_real - yDown_gcvm;
			ped2->SetPos(Point(x2_real, y2_periodic));
		}
	}

	Point distp12 = ped2->GetPos() - ped1->GetPos(); //ped1 ---> ped2
	double Distance = distp12.Norm();
	Point ep12;
	if (Distance >= J_EPS) {
		ep12 = distp12.Normalized();
	}
	else {
		printf("ERROR: \tin AGCVMModel::GetSpacing() ep12 can not be calculated!!!\n");
		Log->Write("WARNING: \tin  AGCVMModel::GetSpacing() ep12 can not be calculated!!!\n");
		Log->Write("\t\t Pedestrians are too near to each other (%f).", Distance);
		exit(EXIT_FAILURE);
	}

	//calculate effective distance
	JEllipse eped1 = ped1->GetEllipse();
	JEllipse eped2 = ped2->GetEllipse();
	double dist;
	double eff_dist = eped1.EffectiveDistanceToEllipse(eped2, &dist);

	// If ped2 has effect on ped1
	double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
	if (condition1 <= 0) 
	{
		ped2->SetPos(ped2_current);
		return  my_pair(FLT_MAX, -1);// ped2 is behind ped1, so not considered
	}

	//Judge conllision
	if (!ped1->GetEllipse().DoesStretch())
	{
		double l = ped1->GetLargerAxis() + ped2->GetLargerAxis();
		double condition2 = ei.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
		condition2 = (condition2 > 0) ? condition2 : -condition2; // abs
		ped2->SetPos(ped2_current);
		if ((condition1 >= 0) && (condition2 <= l / Distance))
			return  my_pair((distp12.Norm() - l), ped2->GetID());
		else
			return  my_pair(FLT_MAX, ped2->GetID());
	}
	else
	{
		double a1 = ped1->GetLargerAxis();
		Point v = ped1->GetV();
		//-------------------------------------------------------------
		// introduce disturbance here can eliminate clogging
		double b1 = ped1->GetSmallerAxis();
		//b1 = ped1->GetEllipse().GetBmin();
		b1 = ped1->GetEllipse().GetBmin();
		//-------------------------------------------------------------
		double a2 = ped2->GetLargerAxis();
		double b2 = ped2->GetSmallerAxis();
		b2 = ped2->GetEllipse().GetBmin();
		double x2 = ped2->GetPos()._x;
		double y2 = ped2->GetPos()._y;
		double cosphi1 = ei.Normalized()._x;
		double sinphi1 = ei.Normalized()._y;
		double cosphi2 = ped2->GetEllipse().GetCosPhi();
		double sinphi2 = ped2->GetEllipse().GetSinPhi();
		//Judge the position of the center of ped2
		double d1 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) + b1;
		double d2 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) - b1;
		if (d1*d2 <= 0)
		{
			//if the center between two lines, collision
			ped2->SetPos(ped2_current);
			return  my_pair(eff_dist, ped2->GetID());
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
		Point A1e = A1.TransformToEllipseCoordinates(ped2->GetPos(), cosphi2, sinphi2);
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
		Te2._x = De2._y / b2;
		Te2._y = -De2._x / a2;
		Point Te1n = Te1.Normalized();
		Point Te2n = Te2.Normalized();
		Point Re1;
		Re1._x = a2 * Te1n._x;
		Re1._y = b2 * Te1n._y;
		Point Re2;
		Re2._x = a2 * Te2n._x;
		Re2._y = b2 * Te2n._y;
		//Calculate distance between point R and line
		double Dis1 = Ne1.ScalarProduct(Re1) - Ne1.ScalarProduct(A1e);
		double Dis2 = Ne2.ScalarProduct(Re2) - Ne2.ScalarProduct(A2e);
		//Judge if the line contact with ellipse2
		ped2->SetPos(ped2_current);
		if (Dis1 >= 0 && Dis2 >= 0)
			return  my_pair(FLT_MAX, -1);
		else
		{
			return  my_pair(eff_dist, ped2->GetID());
		}
	}
}

Point AGCVMModel::ForceRepPed(Pedestrian* ped1, Pedestrian* ped2, Point e0, Point infd, int periodic) const
{
	Point F_rep(0.0, 0.0);
	double x_j = ped2->GetPos()._x;
	double y_j = ped2->GetPos()._y;

	if (periodic) {
		double xLeft_gcvm = GetLeftBoundary();
		double xRight_gcvm = GetRightBoundary();
		double yUp_gcvm = GetUpBoundary();
		double yDown_gcvm = GetDownBoundary();
		double cutoff_gcvm = GetCutoff();
		double x = ped1->GetPos()._x;
		double y = ped1->GetPos()._y;
		if ((xRight_gcvm - x) + (x_j - xLeft_gcvm) <= cutoff_gcvm) {
			double x2_periodic = x_j + xRight_gcvm - xLeft_gcvm;
			ped2->SetPos(Point(x2_periodic, y_j));
		}
		if ((x - xLeft_gcvm) + (xRight_gcvm - x_j) <= cutoff_gcvm) {
			double x2_periodic = xLeft_gcvm - xRight_gcvm + x_j;
			ped2->SetPos(Point(x2_periodic, y_j));
		}
		if ((y - yDown_gcvm) + (yUp_gcvm - y_j) <= cutoff_gcvm) {
			double y2_periodic = yDown_gcvm - yUp_gcvm + y_j;
			ped2->SetPos(Point(x_j, y2_periodic));
		}
		if ((y_j - yDown_gcvm) + (yUp_gcvm - y) <= cutoff_gcvm) {
			double y2_periodic = yUp_gcvm + y_j - yDown_gcvm;
			ped2->SetPos(Point(x_j, y2_periodic));
		}
	}

	Point distp12 = ped2->GetPos() - ped1->GetPos();
	double Distance = distp12.Norm();
	Point ep12;
	double R_ij;
	if (Distance >= J_EPS) 
	{
		ep12 = distp12.Normalized();
	}
	else 
	{
		printf("ERROR: \tin AGCVMModel::forcePedPed() ep12 can not be calculated!!!\n");
		Log->Write(KRED "\nWARNING: \tin AGCVMyModel::forcePedPed() ep12 can not be calculated!!!" RESET);
		Log->Write("\t\t Pedestrians are too near to each other (dist=%f).", Distance);
		Log->Write("\t\t Maybe the value of <a> in force_ped should be increased. Going to exit.\n");
		printf("ped1 %d  ped2 %d\n", ped1->GetID(), ped2->GetID());
		printf("ped1 at (%f, %f), ped2 at (%f, %f)\n", ped1->GetPos()._x, ped1->GetPos()._y, ped2->GetPos()._x, ped2->GetPos()._y);
		exit(EXIT_FAILURE);
	}
	Point ei;
	JEllipse Eped1 = ped1->GetEllipse();
	JEllipse Eped2 = ped2->GetEllipse();
	double dist;
	dist = Eped1.EffectiveDistanceToEllipse(Eped2, &dist);
	ei._x = Eped1.GetCosPhi();
	ei._y = Eped1.GetSinPhi();
	//Vision area-----------------------------------
	double condition1 = e0.ScalarProduct(ep12);
	double condition2 = ei.ScalarProduct(ep12);
	//-----------------------------------------------
	Point ei2;
	ei2._x = Eped2.GetCosPhi();
	ei2._y = Eped2.GetSinPhi();

	//Anticipation
	double S_Gap = 0;
	double Dis_Gap = 0;
	double multi_e0 = ped1->GetV0().ScalarProduct(ped2->GetV0());
	double beta = (3 - multi_e0)/2;
	if (GetAnticipation() == 1)
	{
		double t_anti = GetAntiT();//Anticipation time
		//New S_Gap: Using desired speed instead of real speed
		S_Gap = (ped1->GetV().Normalized().ScalarProduct(ep12)*ped1->GetV0Norm() - ped2->GetV().Normalized().ScalarProduct(ep12)*ped1->GetV0Norm());
		Dis_Gap=S_Gap * t_anti*beta;
	}

	if (GetGCVMU() == 0)
	{
		R_ij = _aPed * exp((-dist) / _DPed);
		F_rep = ep12 * (-R_ij);//ep12: from 1(i) to 2(j)
	}
	else if (condition1 > 0 || condition2 > 0)//rule:pedestrian's direction only influenced by pedestrian in vision area
	{
		double c1 = ped1->GetCooperation();
		double c2 = ped2->GetCooperation();
		int collision = JudgeCollision(ped1, ped2);
		double condition3 = e0.ScalarProduct(ped2->GetV().Normalized());// ped2 move in the same direction of ped1's e0;
		if ((dist < 0.001) && (GetContactRep() == 1))
		{
			double R_dist = dist - Dis_Gap;
			R_ij = _aPed * exp((-R_dist) / _DPed);
			F_rep = ep12 * (-R_ij);// Contact repulision force
		}
		else if ((GetAttracForce() == 1) && condition2>0 && condition3>0 && S_Gap < 0&& dist>0)
		{
			double R_dist = dist+Dis_Gap;
			R_ij = -1*_aPed * exp((-R_dist) / _DPed);
			//F_rep = ep12 * (R_ij);// Attractive force
			F_rep = infd * R_ij;
		}
		else
		{
			double R_dist = dist - Dis_Gap;
			R_ij = _aPed * exp((-R_dist) / _DPed);
			F_rep = infd * R_ij;// Normal repulsion force
		}

	}
	ped2->SetPos(Point(x_j, y_j));
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
	const double EPS = 0.000; 
	double Distance = dist.Norm() + EPS; 
	Point e_iw;
	double R_iw;
	double min_distance_to_wall = 0.00001;
	if (Distance > min_distance_to_wall) {
		e_iw = dist / Distance;
	}
	else {
		Log->Write("WARNING:\t AGCVMModel: forceRepWall() ped %d [%f, %f] is too near to the wall [%f, %f]-[%f, %f] (dist=%f)", ped->GetID(), ped->GetPos()._y, ped->GetPos()._y, w.GetPoint1()._x, w.GetPoint1()._y, w.GetPoint2()._x, w.GetPoint2()._y, Distance);
		Point new_dist = centroid - ped->GetPos(); // ped ---> center
		new_dist = new_dist / new_dist.Norm();
		//printf("new distance = (%f, %f) inside=%d\n", new_dist._x, new_dist._y, inside);
		e_iw = (inside ? new_dist * -1 : new_dist);
	}

	//rule: wall in behind has no influence 
	Point pt1 = w.GetPoint1();
	Point pt2 = w.GetPoint2();
	Point dist_pt1 = pt1 - ped->GetPos();
	Point dist_pt2 = pt2 - ped->GetPos();
	Point e_iw1 = dist_pt1 / dist_pt1.Norm();
	Point e_iw2 = dist_pt2 / dist_pt2.Norm();
	Point ei;
	JEllipse Eped1 = ped->GetEllipse();
	ei._x = Eped1.GetCosPhi();
	ei._y = Eped1.GetSinPhi();
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

	// using new method calculate influence direciton----------------------------------------
	Point inf_direction;// left side f pedestrian
	inf_direction._x = -e0._y;
	inf_direction._y = e0._x;
	inf_direction = inf_direction.Normalized();
	double result1 = e0.CrossProduct(e_iw);
	double result2 = e0.CrossProduct(inf_direction);
	if (bool equal = almostEqual(result1, 0, 0.00001))//Is there any possible that desired direction towards the wall?
	{
		Point zero = Point(0.0, 0.0);
		int random = rand() % 1000;//choose one direciton by random
		if (random > 500)//when random is larger than 50, influence's direction is right, otherwise is left
		{
			inf_direction = zero - inf_direction;
		}
	}
	else
	{
		double result2 = e0.CrossProduct(inf_direction);
		Point zero = Point(0.0, 0.0);
		if (result1*result2 < 0)
		{
			inf_direction = zero - inf_direction;
		}
	}
	//-----------------------------------------------------------------------------------------

	//use effective distance to calculate ForceRepWall
	double effdis = Eped1.EffectiveDistanceToLine(w);
	R_iw = -_aWall * exp((-effdis) / _DWall);
	if (GetGCVMU() == 1)
	{
		F_wrep = inf_direction * R_iw;//new method
	}
	else
	{
		F_wrep = e_iw * R_iw;//original method
	}
	return F_wrep;
}

double AGCVMModel::GetSpacingRoom(Pedestrian* ped, SubRoom* subroom, Point ei) const
{
	double spacing = FLT_MAX;
	double distance = FLT_MAX;
	//first the walls
	for (const auto & wall : subroom->GetAllWalls())
	{
		distance = GetSpacingWall(ped, wall, ei);
		spacing = spacing > distance ? distance : spacing;
	}
	//then the obstacles
	for (const auto & obst : subroom->GetAllObstacles())
	{
		for (const auto & wall : obst->GetAllWalls())
		{
			distance = GetSpacingWall(ped, wall, ei);
			spacing = spacing > distance ? distance : spacing;
		}
	}
	//and finally the closed doors
	for (const auto & goal : subroom->GetAllTransitions())
	{
		if (!goal->IsOpen())
		{
			distance = GetSpacingWall(ped, *(static_cast<Line*>(goal)), ei);
			spacing = spacing > distance ? distance : spacing;
		}
		/*
		//door is open, bur not my door
		int uid1 = goal->GetUniqueID();
		int uid2 = ped->GetExitIndex();
		if ((uid1 != uid2) && (goal->IsOpen() == true))
		{
			distance = GetSpacingWall(ped, *(static_cast<Line*>(goal)), ei);
			spacing = spacing > distance ? distance : spacing;
		}
		*/
	}
	return spacing;

}

double AGCVMModel::GetSpacingWall(Pedestrian* ped, const Line& l, Point ei) const
{
	double spacing = FLT_MAX;
	Point pp = ped->GetPos();
	Point pt = l.ShortestPoint(ped->GetPos());
	Point p1 = l.GetPoint1();
	Point p2 = l.GetPoint2();
	Point dist = pt - pp;
	Point ei_vertical;
	ei_vertical._x = -ei.Normalized()._y;
	ei_vertical._y = ei.Normalized()._x;
	//--------------------------------------------------------
	double b = ped->GetSmallerAxis();
	b = ped->GetEllipse().GetBmin();
	//b = ped->GetEllipse().GetBmax();
	b = 0.1;
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
	double effdis = ped->GetEllipse().EffectiveDistanceToLine(l);
	double cosangle = dist.ScalarProduct(ei) / (dist.Norm()*ei.Norm());
	if (cosangle < 0.00001)
	{
		return spacing;
	}
	spacing = effdis / cosangle;
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
	double normdir = direction.Norm();
	double cosPhi = direction._x / normdir;
	double sinPhi = direction._y / normdir;
	JEllipse e = ped->GetEllipse();
	e.SetCosPhi(cosPhi);
	e.SetSinPhi(sinPhi);
	ped->SetEllipse(e);
	//speed = speed.ScalarProduct(direction) > 0 ? speed : speed * -1;
	ped->SetV(speed);
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

bool AGCVMModel::RealClogging(Pedestrian* ped1, Pedestrian* ped2, Point ei, int periodic) const
{
	//-------------------------------------------------------------
	double limitation = 100;
	double x1 = ped1->GetPos()._x;
	double y1 = ped1->GetPos()._y;
	double x2_real = ped2->GetPos()._x;
	double y2_real = ped2->GetPos()._y;
	Point ped2_current = ped2->GetPos();

	if (periodic)
	{
		double xLeft_gcvm = GetLeftBoundary();
		double xRight_gcvm = GetRightBoundary();
		double yUp_gcvm = GetUpBoundary();
		double yDown_gcvm = GetDownBoundary();
		double cutoff_gcvm = GetCutoff();
		if ((xRight_gcvm - x1) + (x2_real - xLeft_gcvm) <= cutoff_gcvm) {
			double x2_periodic = x2_real + xRight_gcvm - xLeft_gcvm;
			ped2->SetPos(Point(x2_periodic, y2_real));
		}
		if ((x1 - xLeft_gcvm) + (xRight_gcvm - x2_real) <= cutoff_gcvm) {
			double x2_periodic = xLeft_gcvm - xRight_gcvm + x2_real;
			ped2->SetPos(Point(x2_periodic, y2_real));
		}
		if ((y1 - yDown_gcvm) + (yUp_gcvm - y2_real) <= cutoff_gcvm) {
			double y2_periodic = yDown_gcvm - yUp_gcvm + y2_real;
			ped2->SetPos(Point(x2_real, y2_periodic));
		}
		if ((y2_real - yDown_gcvm) + (yUp_gcvm - y1) <= cutoff_gcvm) {
			double y2_periodic = yUp_gcvm + y2_real - yDown_gcvm;
			ped2->SetPos(Point(x2_real, y2_periodic));
		}
	}

	Point distp12 = ped2->GetPos() - ped1->GetPos();
	double Distance = distp12.Norm();
	Point ep12;
	if (Distance >= J_EPS) {
		ep12 = distp12.Normalized();
	}
	else {
		Log->Write("WARNING: \tin AGCVMModel::GetSPacing() ep12 can not be calculated!!!\n");
		Log->Write("\t\t Pedestrians are too near to each other (%f).", Distance);
		exit(EXIT_FAILURE);
	}
	//calculate effective distance
	JEllipse eped1 = ped1->GetEllipse();
	JEllipse eped2 = ped2->GetEllipse();
	double dist;
	double eff_dist = eped1.EffectiveDistanceToEllipse(eped2, &dist);
	double condition1 = ei.ScalarProduct(ep12);
	if (condition1 <= 0) {
		ped2->SetPos(ped2_current);
		return true; // ped2---ped1---Target, may be a real clogging
	}
	if (!ped1->GetEllipse().DoesStretch())
	{
		double r = ped2->GetLargerAxis();
		double l = ped1->GetLargerAxis() + ped2->GetLargerAxis();
		double condition2 = ei.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
		condition2 = (condition2 > 0) ? condition2 : -condition2; // abs
		ped2->SetPos(ped2_current);
		if ((condition1 >= 0) && (condition2 <= r / Distance) && (distp12.Norm() - l > 0.01) && (distp12.Norm() - l < limitation))
			return  false; //1. in front 2. ped1---ped2---Target 3. l+0.01 <distance < l+limitation 
		else
			return  true;
	}
	else
	{
		//Judge conllision
		//Obtain parameters
		double a1 = ped1->GetLargerAxis();
		Point v = ped1->GetV();
		double b1 = 0.001;
		double a2 = ped2->GetLargerAxis();
		double b2 = ped2->GetSmallerAxis();
		double x2 = ped2->GetPos()._x;
		double y2 = ped2->GetPos()._y;
		double cosphi1 = ei.Normalized()._x;
		double sinphi1 = ei.Normalized()._y;
		double cosphi2 = ped2->GetEllipse().GetCosPhi();
		double sinphi2 = ped2->GetEllipse().GetSinPhi();
		//Judge the position of the center of ped2
		double d1 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) + b1;
		double d2 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) - b1;
		if (d1*d2 <= 0) {
			ped2->SetPos(ped2_current);
			if (eff_dist<0.01 || eff_dist>limitation)
			{
				return true;
			}
			return false;
		}
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
		Point A1e = A1.TransformToEllipseCoordinates(ped2->GetPos(), cosphi2, sinphi2);
		Point A2e = A2.TransformToEllipseCoordinates(ped2->GetPos(), cosphi2, sinphi2);
		double J1 = Ne.ScalarProduct(A1e);
		double J2 = Ne.ScalarProduct(A2e);
		Point De1 = (J1 >= 0) ? De * (-1, -1) : De;
		Point De2 = (J2 >= 0) ? De * (-1, -1) : De;
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
		Te2._x = De2._y / b2;
		Te2._y = -De2._x / a2;
		Point Te1n = Te1.Normalized();
		Point Te2n = Te2.Normalized();
		Point Re1;
		Re1._x = a2 * Te1n._x;
		Re1._y = b2 * Te1n._y;
		Point Re2;
		Re2._x = a2 * Te2n._x;
		Re2._y = b2 * Te2n._y;
		double Dis1 = Ne1.ScalarProduct(Re1) - Ne1.ScalarProduct(A1e);
		double Dis2 = Ne2.ScalarProduct(Re2) - Ne2.ScalarProduct(A2e);
		ped2->SetPos(ped2_current);
		if (Dis1 >= 0 && Dis2 >= 0)
			return  true;
		else
		{
			if (eff_dist<0.01 || eff_dist>limitation)
			{
				return true;
			}
			return false;
		}
	}

}

int AGCVMModel::GetAnticipation() const
{
	return _Anticipation;
}

int AGCVMModel::GetContactRep() const
{
	return _ContactRep;
}

int AGCVMModel::GetAttracForce() const
{
	return _AttracForce;
}

double AGCVMModel::GetAntiT() const
{
	return _AntiT;
}

int AGCVMModel::JudgeCollision(Pedestrian* ped1, Pedestrian* ped2) const
{
	int collision=0;
	double At = GetAntiT();
	//At = 0.1;
	Point p1 = ped1->GetPos();
	Point p2 = ped2->GetPos();
	double cosphi1 = ped1->GetEllipse().GetCosPhi();
	double sinphi1 = ped1->GetEllipse().GetSinPhi();
	Point  e1 = Point(cosphi1, sinphi1);
	e1 = e1 * ped2->GetV0Norm();
	//e1 = ped1->GetV();
	double cosphi2 = ped2->GetEllipse().GetCosPhi();
	double sinphi2 = ped2->GetEllipse().GetSinPhi();
	Point  e2 = Point(cosphi2, sinphi2);
	e2 = e2 * ped2->GetV0Norm();
	//e2 = ped2->GetV();
	Point distp12 = p2 - p1; //ped1 ---> ped2
	double distance = distp12.Norm();
	JEllipse Eped1 = ped1->GetEllipse();
	JEllipse Eped2 = ped2->GetEllipse();
	double dist;
	dist = Eped1.EffectiveDistanceToEllipse(Eped2, &dist);
	double r = ped1->GetEllipse().GetBmax();

	double a = (e1._x - e2._x)*(e1._x - e2._x) + (e1._y - e2._y)*(e1._y - e2._y);
	double b = 2 * ((p1._x - p2._x)*(e1._x - e2._x) + (p1._y - p2._y)*(e1._y - e2._y));
	double c = distance * distance - 4 * r * r;
	double delta = b * b - 4 * a * c;

	// Distance is zero, collision
	if (dist < 0.00001)
	{
		collision = 1;
		return collision;
	}

	if (delta == 0)
	{
		double t = -b / (2 * a);
		if (t > 0 && t < At)
		{
			collision = 1;
		}
	}
	else if (delta > 0)
	{
		double sd = sqrt(delta);
		double t1 = (-b - sd) / (2 * a);
		double t2 = (-b + sd) / (2 * a);
		if ((t1 >= 0 && t1 < At) || (t1 < 0 && t2 >= 0 && t2 < At))
		{
			collision = 1;
		}
	}
	return collision;
}

my_pair AGCVMModel::GetSpacing_collision(Pedestrian* ped1, Pedestrian* ped2, Point ei, int periodic) const
{
	double x1 = ped1->GetPos()._x;
	double y1 = ped1->GetPos()._y;
	double x2_real = ped2->GetPos()._x;
	double y2_real = ped2->GetPos()._y;
	Point ped2_current = ped2->GetPos();

	if (periodic)
	{
		double xLeft_gcvm = GetLeftBoundary();
		double xRight_gcvm = GetRightBoundary();
		double yUp_gcvm = GetUpBoundary();
		double yDown_gcvm = GetDownBoundary();
		double cutoff_gcvm = GetCutoff();
		if ((xRight_gcvm - x1) + (x2_real - xLeft_gcvm) <= cutoff_gcvm) {
			double x2_periodic = x2_real + xRight_gcvm - xLeft_gcvm;
			ped2->SetPos(Point(x2_periodic, y2_real));
		}
		if ((x1 - xLeft_gcvm) + (xRight_gcvm - x2_real) <= cutoff_gcvm) {
			double x2_periodic = xLeft_gcvm - xRight_gcvm + x2_real;
			ped2->SetPos(Point(x2_periodic, y2_real));
		}
		if ((y1 - yDown_gcvm) + (yUp_gcvm - y2_real) <= cutoff_gcvm) {
			double y2_periodic = yDown_gcvm - yUp_gcvm + y2_real;
			ped2->SetPos(Point(x2_real, y2_periodic));
		}
		if ((y2_real - yDown_gcvm) + (yUp_gcvm - y1) <= cutoff_gcvm) {
			double y2_periodic = yUp_gcvm + y2_real - yDown_gcvm;
			ped2->SetPos(Point(x2_real, y2_periodic));
		}
	}

	Point distp12 = ped2->GetPos() - ped1->GetPos(); //ped1 ---> ped2
	double Distance = distp12.Norm();
	Point ep12;
	if (Distance >= J_EPS) {
		ep12 = distp12.Normalized();
	}
	else {
		printf("ERROR: \tin AGCVMModel::GetSpacing() ep12 can not be calculated!!!\n");
		Log->Write("WARNING: \tin  AGCVMModel::GetSpacing() ep12 can not be calculated!!!\n");
		Log->Write("\t\t Pedestrians are too near to each other (%f).", Distance);
		exit(EXIT_FAILURE);
	}

	//calculate effective distance
	JEllipse eped1 = ped1->GetEllipse();
	JEllipse eped2 = ped2->GetEllipse();
	double dist;
	double eff_dist = eped1.EffectiveDistanceToEllipse(eped2, &dist);

	// If ped2 has effect on ped1
	double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
	if (condition1 <= 0)
	{
		ped2->SetPos(ped2_current);
		return  my_pair(FLT_MAX, -1);// ped2 is behind ped1, so not considered
	}

	//Judge conllision
	if (!ped1->GetEllipse().DoesStretch())
	{
		double l = 0.15;// specif case in collision
		double condition2 = ei.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
		condition2 = (condition2 > 0) ? condition2 : -condition2; // abs
		ped2->SetPos(ped2_current);
		if ((condition1 >= 0) && (condition2 <= l / Distance))
			return  my_pair((distp12.Norm() - l), ped2->GetID());
		else
			return  my_pair(FLT_MAX, ped2->GetID());
	}
	else
	{
		double a1 = 0.05;
		Point v = ped1->GetV();
		//-------------------------------------------------------------
		// introduce disturbance here can eliminate clogging
		double b1 = ped1->GetSmallerAxis();
		b1 = 0.05;
		//-------------------------------------------------------------
		double a2 = 0.05;
		double b2 = 0.05;
		double x2 = ped2->GetPos()._x;
		double y2 = ped2->GetPos()._y;
		double cosphi1 = ei.Normalized()._x;
		double sinphi1 = ei.Normalized()._y;
		double cosphi2 = ped2->GetEllipse().GetCosPhi();
		double sinphi2 = ped2->GetEllipse().GetSinPhi();
		//Judge the position of the center of ped2
		double d1 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) + b1;
		double d2 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) - b1;
		if (d1*d2 <= 0)
		{
			//if the center between two lines, collision
			ped2->SetPos(ped2_current);
			return  my_pair(eff_dist, ped2->GetID());
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
		Point A1e = A1.TransformToEllipseCoordinates(ped2->GetPos(), cosphi2, sinphi2);
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
		Te2._x = De2._y / b2;
		Te2._y = -De2._x / a2;
		Point Te1n = Te1.Normalized();
		Point Te2n = Te2.Normalized();
		Point Re1;
		Re1._x = a2 * Te1n._x;
		Re1._y = b2 * Te1n._y;
		Point Re2;
		Re2._x = a2 * Te2n._x;
		Re2._y = b2 * Te2n._y;
		//Calculate distance between point R and line
		double Dis1 = Ne1.ScalarProduct(Re1) - Ne1.ScalarProduct(A1e);
		double Dis2 = Ne2.ScalarProduct(Re2) - Ne2.ScalarProduct(A2e);
		//Judge if the line contact with ellipse2
		ped2->SetPos(ped2_current);
		if (Dis1 >= 0 && Dis2 >= 0)
			return  my_pair(FLT_MAX, -1);
		else
		{
			return  my_pair(eff_dist, ped2->GetID());
		}
	}
}