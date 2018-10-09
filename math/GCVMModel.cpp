/**
* \file       GCVM.cpp
* \date        Oct 8, 2018
* \version     v0.7
* \copyright   <2009-2015> Forschungszentrum J¨¹lich GmbH. All rights reserved.
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
* Generalized collision-free velocity model: Qiancheng (3)
*
*
**/

# define NOMINMAX
#include "../pedestrian/Pedestrian.h"
//#include "../routing/DirectionStrategy.h"
#include "../mpi/LCGrid.h"
#include "../geometry/Wall.h"
#include "../geometry/SubRoom.h"


#include "GCVMModel.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads()  1
#endif

double xRight_gcvm = 26.0;
double xLeft_gcvm = 0.0;
double cutoff_gcvm = 2.0;

//--------------------------------
int Vision_area = 1;//Vision_Area=1:open,Vision_Area=0:close,Vision_Area=2:ei_half,Vision=3:e0_half
int Direction_smooth = 1;//Direction_smooth=1:using tau,Direction_smooth=0:not using tau
int Vertical_influence = 1;//Vertical_influence=1:vertical influence,Vertical_influence=0:original influence
int Velocity_influence = 0;//Velocity-influence=1:velocity wll be considered when calculate influence direction
int bf_use = 0;//when calculate direction 1:use b, 0:effective distance
int bv_use = 1;//when calculate velocity 1: use b, 0: effecitve distance
int bmin_use = 1;//width of area in front 1: bmin, 0: b
int real_distance = 0;//when calculate spacing between two pedestrian
					  //-------------------------------

using std::vector;
using std::string;

GCVMModel::GCVMModel(std::shared_ptr<DirectionStrategy> dir, double aped, double Dped,
	double awall, double Dwall, double Ts, double Td)
{
	_direction = dir;
	// Force_rep_PED Parameter
	_aPed = aped;
	_DPed = Dped;
	// Force_rep_WALL Parameter
	_aWall = awall;
	_DWall = Dwall;
	_Ts = Ts;
	_Td = Td;
}


GCVMModel::~GCVMModel()
{

}

bool GCVMModel::Init(Building* building)
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

void GCVMModel::ComputeNextTimeStep(double current, double deltaT, Building* building, int periodic)
{
	// collect all pedestrians in the simulation.
	const vector< Pedestrian* >& allPeds = building->GetAllPedestrians();
	vector<Pedestrian*> pedsToRemove;
	pedsToRemove.reserve(500);
	unsigned long nSize;
	nSize = allPeds.size();
	int nThreads = omp_get_max_threads();
	nThreads = 1; //debug only, I have to use other pedestrian's information
	int partSize;
	partSize = ((int)nSize > nThreads) ? (int)(nSize / nThreads) : (int)nSize;
	if (partSize == (int)nSize)
		nThreads = 1; // not worthy to parallelize
	double v_k = 0;//anticipation parameter, x_future=x_current+v_current*v_k*timesteps
#pragma omp parallel default(shared) num_threads(nThreads)
	{
		vector< Point > result_acc = vector<Point >();
		result_acc.reserve(nSize);
		vector< Point > result_dir = vector<Point >();
		result_dir.reserve(nSize);
		vector< my_pair > spacings = vector<my_pair >();
		spacings.reserve(nSize); // larger than needed
		vector< vector<int> > relations = vector< vector<int> >();
		relations.reserve(nSize);
		vector<int> relation = vector<int>();
		relation.reserve(nSize);
		vector< Point > second_acc = vector<Point>();
		second_acc.reserve(nSize);
		vector<Point> f_pos = vector<Point>();
		f_pos.reserve(nSize);
		const int threadID = omp_get_thread_num();

		int start = threadID * partSize;
		int end;
		end = (threadID < nThreads - 1) ? (threadID + 1) * partSize - 1 : (int)(nSize - 1);
		for (int p = start; p <= end; ++p) {
			Pedestrian* ped = allPeds[p];
			Room* room = building->GetRoom(ped->GetRoomID());
			SubRoom* subroom = room->GetSubRoom(ped->GetSubRoomID());
			Point repPed = Point(0, 0);
			vector<Pedestrian*> neighbours;
			building->GetGrid()->GetNeighbourhood(ped, neighbours);
			Point inid_direction = e0(ped, room);
			double alpha = 1;//may be useful in the future

			//using a new method calculate the influence of pedestrian (the value of influence id decided by distance and the direction is vertical with desired direction) 
			Point inf_direction;// left side of pedestrian
			inf_direction._x = -inid_direction._y;
			inf_direction._y = inid_direction._x;
			inf_direction = inf_direction.Normalized();
			//--------------------------------------------------------------------------------------------------------------------------------------------------------------

			int size = (int)neighbours.size();
			for (int i = 0; i < size; i++)
			{
				Pedestrian* ped1 = neighbours[i];
				//if they are in the same subroom
				Point p1 = ped->GetPos();
				Point p2 = ped1->GetPos();
				Point ep12 = p2 - p1;

				//Todo:Anticipation(will be used in the future,when v_k = 0 the anticipation is useless)--------------
				Point ped1_current = p2;
				Point ped1_future = p2 + ped1->GetV()*deltaT*v_k;
				ped1->SetPos(ped1_future);
				//----------------------------------------------------------------------------------------------------

				//Deciding the influence direction--------------------------------------------------------------------
				double result1 = inid_direction.CrossProduct(ep12);
				Point zero = Point(0, 0);
				if (bool equal = almostEqual(result1, 0, 0.001))//if neighbour is in front or behind pedestrian
				{
					int random = rand() % 100;//choose one direciton bu random
					if (random > 200)//when random is larger than 50, influence's direction is right, otherwise is left
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
				//test code
				/*
				if (ped->GetID() == 14)
				{

				fprintf(stderr, "\npedID=%d,ped.x=%f,ped.y=%f\n",ped->GetID(), ped->GetPos()._x, ped->GetPos()._y);
				fprintf(stderr, "ped1ID=%d,ped1.x=%f,ped.y=%f\n", ped1->GetID(), ped1->GetPos()._x, ped1->GetPos()._y);
				fprintf(stderr, "inf_direction._x=%f,inf_direction._y=%f,inid_direction._x=%f,inid_direction._y=%f\n",
				inf_direction._x, inf_direction._y, inid_direction._x, inid_direction._y);
				}
				*/
				//----------------------------------------------------------------------------------------------------

				//subrooms to consider when looking for neighbour for the 3d visibility
				vector<SubRoom*> emptyVector;
				emptyVector.push_back(subroom);
				emptyVector.push_back(building->GetRoom(ped1->GetRoomID())->GetSubRoom(ped1->GetSubRoomID()));
				bool isVisible = building->IsVisible(p1, p2, emptyVector, false);
				if (!isVisible)
				{
					ped1->SetPos(ped1_current);
					continue;
				}
				if (ped->GetUniqueRoomID() == ped1->GetUniqueRoomID()) {
					if (Vertical_influence == 1)
					{
						repPed += inf_direction * ForceRepPed(ped, ped1, inid_direction, periodic).Norm();//new method
					}
					else
					{
						repPed += ForceRepPed(ped, ped1, inid_direction, periodic);//original method
					}
				}
				else {
					// or in neighbour subrooms
					SubRoom* sb2 = building->GetRoom(ped1->GetRoomID())->GetSubRoom(ped1->GetSubRoomID());
					if (subroom->IsDirectlyConnectedWith(sb2)) {
						if (Vertical_influence == 1)
						{
							repPed += inf_direction * ForceRepPed(ped, ped1, inid_direction, periodic).Norm();//new method
						}
						else
						{
							repPed += ForceRepPed(ped, ped1, inid_direction, periodic);//original method
						}
					}
				}
				ped1->SetPos(ped1_current);
			} //for i

			  // Calculating influence of walls-------------------------------------------------------------------------
			Point pdesire_direction = inid_direction + repPed;// no use now, may be useful in the future
			Point repWall = ForceRepRoom(allPeds[p], subroom, inid_direction, pdesire_direction);
			//-----------------------------------------------------------------------------------------------------------

			//Caluculating desired direcition----------------------------------------------------------------------------------------------
			Point d_direction;
			if (Direction_smooth == 1 && Vertical_influence == 3)
			{
				//original method, alpha smaller than 1 can avoid move back
				Point repSum = repPed + repWall;
				double abs_repSum = repSum.Norm();
				repSum = abs_repSum > alpha ? repSum * (alpha / abs_repSum) : repSum;
				d_direction = inid_direction + repSum;
			}
			else
			{
				d_direction = inid_direction + repPed + repWall;//new method
			}
			//------------------------------------------------------------------------------------------------------------------------------

			//Calculating the actual direction of pedestrian at next timestep------------------------------------------------------------------
			Point a_direction;
			a_direction._x = ped->GetEllipse().GetCosPhi();
			a_direction._y = ped->GetEllipse().GetSinPhi();
			double angle_tau = _Td;
			//Two ways to calculate new direction (vector and angle)
			/*
			//case 1 (vector)
			Point angle_v = (d_direction.Normalized()-a_direction)/ angle_tau;
			Point direction =a_direction+angle_v*deltaT ;
			*/
			//case 2 (angle)
			double angle_dd = atan2(d_direction._y, d_direction._x);
			double angle_ad = atan2(a_direction._y, a_direction._x);
			double angle_diff = angle_dd - angle_ad;
			//angle_diff = angle_diff > 3.1415926 ? angle_diff - 3.1415926 * 2 : angle_diff;
			//angle_diff = angle_diff < -3.1414926 ? angle_diff + 3.1415926 * 2 : angle_diff;
			double angle_velocity = angle_diff / angle_tau;
			double angle_d = angle_ad + angle_velocity * deltaT;
			if (angle_diff > 3.1415926 || angle_diff < -3.1415926)
			{
				angle_d = angle_dd;
			}
			Point direction;
			direction._x = cos(angle_d);
			direction._y = sin(angle_d);

			/*
			//Todo:when v=0, peddestrian change direction directly.
			if (ped->GetV().Norm() < 0.01)
			{
			direction = d_direction;
			}
			*/

			/*
			//test code
			if (ped->GetID() == 11)
			{
			fprintf(stderr, "\n pedID=%d,postion=(%f,%f),e0=(%f,%f)\n",ped->GetID(),ped->GetPos()._x,ped->GetPos()._y,inid_direction._x,inid_direction._y);
			fprintf(stderr, "a_direction=(%f,%f),d_direction=(%f,%f),direction=(%f,%f)\n", a_direction._x, a_direction._y, d_direction._x, d_direction._y, direction._x, direction._y);
			fprintf(stderr, "reped=(%f,%f),repwall=(%f,%f)\n", repPed._x, repPed._y, repWall._x, repWall._y);
			}
			*/
			if (Direction_smooth == 0)
			{
				direction = d_direction;//original method
			}
			result_dir.push_back(direction);
			//------------------------------------------------------------------------------------------------------------------------------------


			//Calculating spacing in front -------------------------------------------------------------------------------------------------------
			for (int i = 0; i < size; i++) {
				Pedestrian* ped1 = neighbours[i];
				// calculate spacing
				if (ped->GetUniqueRoomID() == ped1->GetUniqueRoomID()) {
					spacings.push_back(GetSpacing(ped, ped1, direction, v_k*deltaT, periodic));
				}
				else {
					// or in neighbour subrooms
					SubRoom* sb2 = building->GetRoom(ped1->GetRoomID())->GetSubRoom(ped1->GetSubRoomID());
					if (subroom->IsDirectlyConnectedWith(sb2)) {
						spacings.push_back(GetSpacing(ped, ped1, direction, v_k*deltaT, periodic));
					}
				}
			}//for i
			 //------------------------------------------------------------------------------------------------------------------------------------

			 // Calculate min spacing, and save all spacing=0 to vector relation--------------------------------------------------------------------
			std::sort(spacings.begin(), spacings.end(), sort_pred_gcvm());
			double spacing = spacings.size() == 0 ? 100 : spacings[0].first;
			double first_ID = spacings.size() == 0 ? -1 : spacings[0].second;
			double sec_spacing = spacings.size() <= 1 ? 100 : spacings[1].first;
			relation.push_back(ped->GetID());
			relation.push_back(first_ID);
			if (spacings.size() > 1)
			{
				for (int i = 1; i < spacings.size(); i++)
				{
					if (spacings[i].first < 0.001)
					{
						relation.push_back(spacings[i].second);
					}
				}
			}
			// add this part to avoid pedestrian cross the wall directly
			spacing = spacing > GetSpacingRoom(ped, subroom, direction) ? GetSpacingRoom(ped, subroom, direction) : spacing;
			relations.push_back(relation);
			relation.clear();
			Point speed;
			speed = direction.NormalizedMolified() *OptimalSpeed(ped, spacing);
			result_acc.push_back(speed);
			Point sec_speed;
			sec_speed = direction.NormalizedMolified()*OptimalSpeed(ped, sec_spacing);
			second_acc.push_back(sec_speed);
			spacings.clear(); //clear for ped p
							  //-----------------------------------------------------------------------------------------------------------------------------------------
		} // for p

		  //Coopertion to solve block---------------------------------------------------------------------------------------------------------------------
		  /*
		  //TODO:cooperation according to future position
		  for (int p = start; p <= end; ++p)
		  {
		  int k = 2;
		  Pedestrian* fped = allPeds[p];
		  Point f_v = result_acc[p - start];
		  Point f_dir = result_dir[p - start];
		  f_pos.push_back(fped->GetPos() + f_v * deltaT*k);
		  }
		  for (int p = start; p <= end; ++p)
		  {
		  Pedestrian *ped = allPeds[p];
		  Point p_pos = f_pos[p];
		  for (int i = start; i <= end; ++i)
		  {
		  if (i != p)
		  {
		  Pedestrian *iped = allPeds[i];
		  if ((p_pos - f_pos[i]).Norm() < 0.5)
		  {
		  if (result_acc[p].Norm() < result_acc[i].Norm())
		  {
		  Point zero = Point(0, 0);
		  result_acc[p] = zero;
		  break;

		  }
		  if (result_acc[p].Norm() <= result_acc[i].Norm())
		  {
		  if (ped->GetID() > iped->GetID())//Whose ID is bigger has to wait
		  {
		  Point zero = Point(0, 0);
		  result_acc[p] = zero;
		  break;
		  }
		  }
		  }
		  }
		  }
		  }
		  */
		  /*
		  //solve the block(cooperation with all pedestrians who has influence)
		  for (int p = start; p <= end; ++p)
		  {
		  Pedestrian* ped = allPeds[p];
		  Point v1 = result_acc[p - start];
		  int ID1 = relations[p - start][0];
		  for (int i = 0; i < nSize; i++)
		  {
		  int ID3 = relations[i][0];
		  for (int j = 1; j < relations[p - start].size(); j++)
		  {
		  int ID2 = relations[p - start][j];
		  if (ID2 == ID3)
		  {
		  for (int k = 1; k < relations[i].size(); k++)
		  {
		  int ID4 = relations[i][k];
		  Point v2 = result_acc[i];
		  Point v1_neu = second_acc[p - start];
		  Point v2_neu = second_acc[i];
		  if (ID1 == ID4)// two peddestrian influence eachother
		  {
		  //ped->SetSpotlight(true);// when conflict happen, turn red
		  if (v1_neu.Norm()> v2_neu.Norm())// who's second velocity is bigger go first, others need to wait.
		  {
		  result_acc[p - start] = second_acc[p - start];
		  }
		  else
		  {
		  Point zero = Point(0, 0);
		  result_acc[p - start] = zero;
		  }
		  break;
		  }
		  }

		  }

		  }
		  }
		  }
		  */
		  //--------------------------------------------------------------------------------------------------------------------------------------------------

#pragma omp barrier
		  // update
		for (int p = start; p <= end; ++p) {
			Pedestrian* ped = allPeds[p];

			Point v_neu = result_acc[p - start];
			Point dir_neu = result_dir[p - start];
			Point pos_neu = ped->GetPos() + v_neu * deltaT;
			//Jam is based on the current velocity
			if (v_neu.Norm() >= ped->GetV0Norm()*0.5) {
				ped->ResetTimeInJam();
			}
			else {
				ped->UpdateTimeInJam();
			}
			//calculate ellipse orientation
			double normdir = dir_neu.Norm();
			double cosPhi = dir_neu._x / normdir;
			double sinPhi = dir_neu._y / normdir;
			JEllipse e = ped->GetEllipse();
			e.SetCosPhi(cosPhi);
			e.SetSinPhi(sinPhi);
			ped->SetEllipse(e);

			ped->SetV(v_neu);
			ped->SetPos(pos_neu);
			if (periodic) {
				if (ped->GetPos()._x >= xRight_gcvm) {
					ped->SetPos(Point(ped->GetPos()._x - (xRight_gcvm - xLeft_gcvm), ped->GetPos()._y));
					//ped->SetID( ped->GetID() + 1);
				}
			}
		}
	}//end parallel

	 // remove the pedestrians that have left the building
	for (unsigned int p = 0; p<pedsToRemove.size(); p++) {
		building->DeletePedestrian(pedsToRemove[p]);
	}
	pedsToRemove.clear();
}

Point GCVMModel::e0(Pedestrian* ped, Room* room) const
{
	const Point target = _direction->GetTarget(room, ped); // target is where the ped wants to be after the next timestep
	Point desired_direction;
	const Point pos = ped->GetPos();
	double dist = ped->GetExitLine()->DistTo(pos);
	// check if the molified version works
	Point lastE0 = ped->GetLastE0();
	ped->SetLastE0(target - pos);

	if ((dynamic_cast<DirectionFloorfield*>(_direction.get())) ||
		(dynamic_cast<DirectionLocalFloorfield*>(_direction.get())) ||
		(dynamic_cast<DirectionSubLocalFloorfield*>(_direction.get()))) {
		if (dist > 50 * J_EPS_GOAL) {
			desired_direction = target - pos; //ped->GetV0(target);
		}
		else {
			desired_direction = lastE0;
			ped->SetLastE0(lastE0); //keep old vector (revert set operation done 9 lines above)
		}
	}
	else if (dist > J_EPS_GOAL) {
		desired_direction = ped->GetV0(target);
	}
	else {
		ped->SetSmoothTurning();
		desired_direction = ped->GetV0();
	}
	//test result
	//fprintf(stderr,"\n %d %f %f %f %f",ped->GetID(),ped->GetPos()._x,ped->GetPos()._y,target._x,target._y);
	//test result
	return desired_direction;
}


double GCVMModel::OptimalSpeed(Pedestrian* ped, double spacing) const
{
	double v0 = ped->GetV0Norm();
	double T = _Ts;
	//spacing is approximate here
	double speed = (spacing) / T;
	speed = (speed>0) ? speed : 0;
	speed = (speed<v0) ? speed : v0;
	//      (1-winkel)*speed;
	//todo use winkel
	return speed;
}

// return spacing and id of the nearest pedestrian
my_pair GCVMModel::GetSpacing(Pedestrian* ped1, Pedestrian* ped2, Point ei, double k_deltaT, int periodic) const
{
	//--------------------------------------------------------------
	//Anticipation
	Point ped2_pos = ped2->GetPos();
	Point ped2_current = ped2_pos;
	Point ped2_future = ped2_pos + ped2->GetV()*k_deltaT;
	ped2->SetPos(ped2_future);
	//---------------------------------------------------------------

	double x1 = ped1->GetPos()._x;
	double x2_real = ped2->GetPos()._x;
	double y2 = ped2->GetPos()._y;

	if (periodic) {
		if ((xRight_gcvm - x1) + (x2_real - xLeft_gcvm) <= cutoff_gcvm) {
			double x2_periodic = x2_real + xRight_gcvm - xLeft_gcvm;
			ped2->SetPos(Point(x2_periodic, y2));
		}
	}

	Point distp12 = ped2->GetPos() - ped1->GetPos(); // inversed sign
	double Distance = distp12.Norm();
	Point ep12;
	if (Distance >= J_EPS) {
		ep12 = distp12.Normalized();
	}
	else {
		//printf("ERROR: \tin VelocityModel::forcePedPed() ep12 can not be calculated!!!\n");
		Log->Write("WARNING: \tin VelocityModel::GetSPacing() ep12 can not be calculated!!!\n");
		Log->Write("\t\t Pedestrians are too near to each other (%f).", Distance);
		exit(EXIT_FAILURE);
	}
	//calculate effective distance
	JEllipse eped1 = ped1->GetEllipse();
	JEllipse eped2 = ped2->GetEllipse();
	double dist;
	double eff_dist = eped1.EffectiveDistanceToEllipse(eped2, &dist);
	if (bv_use == 1)
	{
		//eff_dist = Distance - 2 * eped1.GetBmax();
		eff_dist = Distance - eped1.GetEB() - eped2.GetEB();
	}
	//using real distance
	if (real_distance == 1)
	{
		double cosangle = ei.ScalarProduct(distp12) / (ei.Norm()*distp12.Norm());
		if (cosangle <= 0.0000001)
		{
			return my_pair(FLT_MAX, -1);
		}
		eff_dist = eff_dist / cosangle;
	}

	double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
												/*
												//test code
												if (ped1->GetID() == 15 && ped2->GetID() == 11)
												{
												fprintf(stderr, "\nped1ID=%d,ped1.x=%f,ei.x=%f,ei.y=%f,ped2ID=%d,ped2.x=%f,ep12.x=%f,ep12.y=%f,condition1=%f", ped1->GetID(), ped1->GetPos()._x, ei._x,ei._y, ped2->GetID(), ped2->GetPos()._x,ep12._x,ep12._y,condition1);
												}
												*/
	if (condition1 < 0) {
		ped2->SetPos(ped2_current);
		return  my_pair(FLT_MAX, -1);
	}
	//When condition1==0,should no influence on velocity
	bool parallel = almostEqual(condition1, 0.0, J_EPS);
	if (parallel)
	{
		ped2->SetPos(ped2_current);
		return  my_pair(FLT_MAX, -1);
	}
	//Judge conllision
	//Obtain parameters
	double a1 = ped1->GetLargerAxis();
	Point v = ped1->GetV();
	double b1 = ped1->GetSmallerAxis();
	if (bmin_use == 1)
	{
		b1 = ped1->GetEllipse().GetBmin();
	}
	//Avoid block (drill,small tricks)-----------------------------------------------------
	/*
	if (fabs(v._x) < J_EPS && fabs(v._y) < J_EPS) // v==0
	{
	const Point& pos = ped1->GetPos();
	double distGoal = ped1->GetExitLine()->DistToSquare(pos);
	if (distGoal < 0.5)
	{
	b1 = 0.1;
	//b1 = ped1->GetEllipse().GetBmin();

	}
	}
	*/
	//-----------------------------------------------------------------------------------

	double a2 = ped2->GetLargerAxis();
	double b2 = ped2->GetSmallerAxis();
	double x2 = ped2->GetPos()._x;
	double y1 = ped1->GetPos()._y;
	double cosphi1 = ei.Normalized()._x;
	double sinphi1 = ei.Normalized()._y;
	double cosphi2 = ped2->GetEllipse().GetCosPhi();
	double sinphi2 = ped2->GetEllipse().GetSinPhi();
	//Judge the position of the center of ped2
	double d1 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) + b1;
	double d2 = -sinphi1 * (x2 - x1) + cosphi1 * (y2 - y1) - b1;
	if (d1*d2 <= 0) {
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

	/*
	//test code
	printf("\nID1=%d,ID2=%d,ped1x=%f,ped1y=%f,ped2x=%f,ped2y=%f\n", ped1->GetID(), ped2->GetID(), ped1->GetPos()._x, ped1->GetPos()._y, ped2->GetPos()._x, ped2->GetPos()._y);
	printf("\na1=%f,b1=%f,a2=%f,b2=%f,cosphi1=%f,sinphi1=%f,cosphi2=%f,sinphi2=%f\n", a1, b1, a2, b2, cosphi1, sinphi1, cosphi2, sinphi2);
	printf("\nDx=%f,Dy=%f\n", D._x, D._y);
	printf("\nDex=%f,Dey=%f\n", De._x, De._y);
	printf("\nNex=%f,Ney=%f\n", Ne._x, Ne._y);
	printf("\nA1x=%f,A1y=%f\n", A1._x, A1._y);
	printf("\nA2x=%f,A2y=%f\n", A2._x, A2._y);
	printf("\nA1ex=%f,A1ey=%f\n", A1e._x, A1e._y);
	printf("\nA2ex=%f,A2ey=%f\n", A2e._x, A2e._y);
	printf("\J1=%f,J2=%f\n", J1, J2);
	printf("\nDe1x=%f,De1y=%f\n", De1._x, De1._y);
	printf("\nDe2x=%f,De2y=%f\n", De2._x, De2._y);
	printf("\nNe1x=%f,Ne1y=%f\n", Ne1._x, Ne1._y);
	printf("\nNe2x=%f,Ne2y=%f\n", Ne2._x, Ne2._y);
	printf("\nTe1x=%f,Te1y=%f\n", Te1._x, Te1._y);
	printf("\nTe2x=%f,Te2y=%f\n", Te2._x, Te2._y);
	printf("\nTe1nx=%f,Te1ny=%f\n", Te1n._x, Te1n._y);
	printf("\nTe2nx=%f,Te2ny=%f\n", Te2n._x, Te2n._y);
	printf("\nRe1x=%f,Re1y=%f\n", Re1._x, Re1._y);
	printf("\nRe2x=%f,Re2y=%f\n", Re2._x, Re2._y);
	printf("\nDis1=%f,Dis2=%f\n", Dis1, Dis2);
	*/

	//Judge if the line contact with ellipse2
	ped2->SetPos(ped2_current);
	if (Dis1 >= 0 && Dis2 >= 0)
		return  my_pair(FLT_MAX, -1);
	else
	{
		return  my_pair(eff_dist, ped2->GetID());
	}
}
Point GCVMModel::ForceRepPed(Pedestrian* ped1, Pedestrian* ped2, Point e0, int periodic) const
{
	Point F_rep(0.0, 0.0);
	// x- and y-coordinate of the distance between p1 and p2

	double x_j = ped2->GetPos()._x;
	double y_j = ped2->GetPos()._y;

	if (periodic) {
		double x = ped1->GetPos()._x;
		if ((xRight_gcvm - x) + (x_j - xLeft_gcvm) <= cutoff_gcvm) {
			double x2_periodic = x_j + xRight_gcvm - xLeft_gcvm;
			ped2->SetPos(Point(x2_periodic, y_j));
		}
	}

	Point distp12 = ped2->GetPos() - ped1->GetPos();
	double Distance = distp12.Norm();
	Point ep12; // x- and y-coordinate of the normalized vector between p1 and p2
	double R_ij;

	if (Distance >= J_EPS) {
		ep12 = distp12.Normalized();
	}
	else {
		//printf("ERROR: \tin VelocityModel::forcePedPed() ep12 can not be calculated!!!\n");
		Log->Write(KRED "\nWARNING: \tin VelocityModel::forcePedPed() ep12 can not be calculated!!!" RESET);
		Log->Write("\t\t Pedestrians are too near to each other (dist=%f).", Distance);
		Log->Write("\t\t Maybe the value of <a> in force_ped should be increased. Going to exit.\n");
		printf("ped1 %d  ped2 %d\n", ped1->GetID(), ped2->GetID());
		printf("ped1 at (%f, %f), ped2 at (%f, %f)\n", ped1->GetPos()._x, ped1->GetPos()._y, ped2->GetPos()._x, ped2->GetPos()._y);
		exit(EXIT_FAILURE);
	}
	Point ei;
	JEllipse Eped1 = ped1->GetEllipse();
	ei._x = Eped1.GetCosPhi();
	ei._y = Eped1.GetSinPhi();

	//Version area-----------------------------------
	double condition1 = e0.ScalarProduct(ep12); // < e_i , e_ij > should be positive
	double condition2 = ei.ScalarProduct(ep12);
	if (Vision_area == 0)//all
	{
		condition1 = 1;
		condition2 = 1;
	}
	else if (Vision_area == 2)//ei_half
	{
		condition1 = -1;
	}
	else if (Vision_area == 3)//e0_half
	{
		condition2 = -1;
	}
	//-----------------------------------------------

	//rule:pedestrian's direction only influenced by pedestrian in version area
	if (condition1 > 0 || condition2 > 0)
	{
		double dist;
		JEllipse Eped2 = ped2->GetEllipse();
		dist = Eped1.EffectiveDistanceToEllipse(Eped2, &dist);
		if (bf_use == 1)
		{
			//dist = Distance-2*ped1->GetEllipse().GetBmax();
			dist = Distance - Eped1.GetEB() - Eped2.GetEB();
		}
		//---------------------
		if (Velocity_influence == 1)
		{
			Point vped2 = ped2->GetV();
			Point vped1 = ped1->GetV();
			Point v1_influence = vped1.ScalarProduct(ep12) / ep12.Norm();
			Point v2_influence = vped2.ScalarProduct(ep12) / ep12.Norm();
			dist = dist - v1_influence.Norm() + v2_influence.Norm();
		}
		//---------------------
		// dist is approximate here
		dist = dist > 0 ? dist : 0;
		R_ij = -_aPed * exp((-dist) / _DPed);
		F_rep = ep12 * R_ij;
		/*
		//test code
		if (ped1->GetID() == 40)
		{
		fprintf(stderr, "\npedID=%d,ped2ID=%d,pedpos=(%f,%f),effdis=%f,aped=%f,Dped=%f,Rij=%f\n", ped1->GetID(),ped2->GetID(), ped1->GetPos()._x, ped1->GetPos()._y, dist, _aPed, _DPed, R_ij);
		//fprintf(stderr, "pedID=%d,ped.x=%f,ped.y=%f,ep_12._x=%f,ep_12._y=%f \n", ped1->GetID(), ped1->GetPos()._x, ped1->GetPos()._y, ep12._x, ep12._y);
		}
		*/
	}
	ped2->SetPos(Point(x_j, y_j));
	return F_rep;
}//END Velocity:ForceRepPed()

Point GCVMModel::ForceRepRoom(Pedestrian* ped, SubRoom* subroom, Point e0, Point pdesire) const
{
	Point f(0., 0.);
	const Point& centroid = subroom->GetCentroid();
	bool inside = subroom->IsInSubRoom(centroid);
	//first the walls
	for (const auto & wall : subroom->GetAllWalls())
	{
		f += ForceRepWall(ped, wall, centroid, inside, e0, pdesire);
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
				f += ForceRepWall(ped, wall, centroid, inside, e0, pdesire);
			}
	}

	// and finally the closed doors
	for (const auto & goal : subroom->GetAllTransitions())
	{
		if (!goal->IsOpen())
		{
			f += ForceRepWall(ped, *(static_cast<Line*>(goal)), centroid, inside, e0, pdesire);
		}
		int uid1 = goal->GetUniqueID();
		int uid2 = ped->GetExitIndex();
		// ignore my transition consider closed doors
		//closed doors are considered as wall
		//door is open, but it's not my door (has influence)
		/*
		if((uid1 != uid2) && (goal->IsOpen()==true ))
		{
		f +=  ForceRepWall(ped,*(static_cast<Line*>(goal)), centroid, inside, e0, pdesire);
		}
		*/
	}
	return f;
}

Point GCVMModel::ForceRepWall(Pedestrian* ped, const Line& w, const Point& centroid, bool inside, Point e0, Point pdesire) const
{
	Point F_wrep = Point(0.0, 0.0);
	Point pt = w.ShortestPoint(ped->GetPos());
	Point dist = pt - ped->GetPos(); // x- and y-coordinate of the distance between ped and p
	const double EPS = 0.000; // molified see Koester2013
	double Distance = dist.Norm() + EPS; // distance between the centre of ped and point p
										 //double vn = w.NormalComp(ped->GetV()); //normal component of the velocity on the wall
	Point e_iw; // x- and y-coordinate of the normalized vector between ped and pt
				//double K_iw;
	double R_iw;
	double min_distance_to_wall = 0.001; // 10 cm
	if (Distance > min_distance_to_wall) {
		e_iw = dist / Distance;
	}
	else {
		Log->Write("WARNING:\t Velocity: forceRepWall() ped %d [%f, %f] is too near to the wall [%f, %f]-[%f, %f] (dist=%f)", ped->GetID(), ped->GetPos()._y, ped->GetPos()._y, w.GetPoint1()._x, w.GetPoint1()._y, w.GetPoint2()._x, w.GetPoint2()._y, Distance);
		Point new_dist = centroid - ped->GetPos();
		new_dist = new_dist / new_dist.Norm();
		//printf("new distance = (%f, %f) inside=%d\n", new_dist._x, new_dist._y, inside);
		e_iw = (inside ? new_dist * -1 : new_dist);
	}

	//rules--------------------------------------
	/*
	//rule1: when pedestrian very close to the exit, the wall is no influence
	const Point& pos = ped->GetPos();
	double distGoal = ped->GetExitLine()->DistToSquare(pos);
	if(distGoal < J_EPS_GOAL*J_EPS_GOAL)
	return F_wrep;
	*/

	/*
	//rule2: wall's influence is 0 when pedestrian out the range of the wall
	const Point& Point1 = w.GetPoint1();
	const Point& Point2 = w.GetPoint2();
	Point l_direction = Point1 - Point2;
	if (fabs(l_direction.ScalarProduct(ped->GetPos() - pt)) > J_EPS) {
	Point v_direction;
	v_direction._x = -l_direction._y;
	v_direction._y = l_direction._x;
	Point lPoint1 = Point1 - v_direction * 10;
	Point rPoint1 = Point1 + v_direction * 10;
	Point lPoint2 = Point2 - v_direction * 10;
	Point rPoint2 = Point2 + v_direction * 10;
	Line line1 = Line(lPoint1, rPoint1);
	Line line2 = Line(lPoint2, rPoint2);
	JEllipse Eped = ped->GetEllipse();
	double effdis1 = Eped.EffectiveDistanceToLine(line1);
	double effdis2 = Eped.EffectiveDistanceToLine(line2);
	if (effdis1*effdis2 >= 0)
	{
	return F_wrep;
	}
	}
	*/

	/*
	// rule3:pedestrian will not influence by walls behind and walls parallel
	double tmp;
	double judge;
	JEllipse Eped = ped->GetEllipse();
	Point vd;
	vd._x = Eped.GetCosPhi();
	vd._y = Eped.GetSinPhi();
	tmp = pdesire.ScalarProduct(e_iw); // < v_i , e_ij >;
	judge = (tmp + fabs(tmp));
	if (!judge)
	return F_wrep;
	*/

	/*
	//rule4:when velocity is very small, the influence of wall is 0
	Point v = ped->GetV();
	if (fabs(v._x) < J_EPS && fabs(v._y) < J_EPS) // v==0
	return F_wrep;
	*/

	//rule5: wall in behind has no influence 
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
	//version area----------------------------------
	double result_e01 = e0.ScalarProduct(e_iw1);
	double result_e02 = e0.ScalarProduct(e_iw2);
	double result_ei1 = ei.ScalarProduct(e_iw1);
	double result_ei2 = ei.ScalarProduct(e_iw2);
	//----------------------------------------------
	if (Vision_area == 0)//all
	{
		result_e01 = 1;
		result_e02 = 1;
		result_ei1 = 1;
		result_ei2 = 1;
	}
	else if (Vision_area == 2)//ei_half
	{
		result_e01 = -1;
		result_e02 = -1;
	}
	else if (Vision_area == 3)//e0_half
	{
		result_ei1 = -1;
		result_ei2 = -1;
	}
	if (result_e01 < 0 && result_e02 < 0 && result_ei1 < 0 && result_ei1 < 0)
	{
		return F_wrep;
	}

	//rules end------------------------------------------------------------------------------

	//------------------------------------------------------

	// using new method calculate influence direciton----------------------------------------
	Point inf_direction;// left side f pedestrian
	inf_direction._x = -e0._y;
	inf_direction._y = e0._x;
	inf_direction = inf_direction.Normalized();
	double result1 = e0.CrossProduct(e_iw);
	double result2 = e0.CrossProduct(inf_direction);
	if (bool equal = almostEqual(result1, 0, 0.001))//Is there any possible that desired direction towards the wall?
	{
		//double alpha = ped->GetAlpha();
		Point zero = Point(0.0, 0.0);
		int random = rand() % 100;//choose one direciton bu random
		if (random > 200)//when random is larger than 50, influence's direction is right, otherwise is left
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

	/*
	// TODO:The influence direction of wall is parellel with wall
	Point e_w = pt1 - pt2;
	Point inf_direction = e_w.Normalized();
	if (e_w.ScalarProduct(e0) > 0)
	{
	Point zero = Point(0.0, 0.0);
	inf_direction = zero - inf_direction;
	}
	*/

	//use effective distance to calculate ForceRepWall
	double effdis = Eped1.EffectiveDistanceToLine(w);
	if (bf_use == 1)
	{
		//effdis = Distance - Eped1.GetBmax();
		effdis = Distance - Eped1.GetEB();
	}
	// effdist is approximate here
	effdis = effdis > 0 ? effdis : 0;
	R_iw = -_aWall * exp((-effdis) / _DWall);
	if (Vertical_influence == 1)
	{
		F_wrep = inf_direction * R_iw;//new method
	}
	else
	{
		F_wrep = e_iw * R_iw;//original method
	}

	//test code
	/*
	if (ped->GetID() == 11)
	{
	fprintf(stderr, "\pedID=%d,ped.x=%f,ped.y=%f,effdis=%f,awall=%f,Dwall=%f,Riw=%f\n", ped->GetID(), ped->GetPos()._x, ped->GetPos()._y, effdis, _aWall, _DWall, R_iw);
	fprintf(stderr, "\pedID=%d,ped.x=%f,ped.y=%f,e_iw._x=%f,e_iw._y=%f \n\n", ped->GetID(), ped->GetPos()._x, ped->GetPos()._y, e_iw._x, e_iw._y);
	}
	*/
	return F_wrep;
}

double GCVMModel::GetSpacingRoom(Pedestrian* ped, SubRoom* subroom, Point ei) const
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
		int uid1 = goal->GetUniqueID();
		int uid2 = ped->GetExitIndex();
		//door is open, bur not my door
		if ((uid1 != uid2) && (goal->IsOpen() == true))
		{
			distance = GetSpacingWall(ped, *(static_cast<Line*>(goal)), ei);
			spacing = spacing > distance ? distance : spacing;
		}
	}
	return spacing;

}

double GCVMModel::GetSpacingWall(Pedestrian* ped, const Line& l, Point ei) const
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
	double b = ped->GetSmallerAxis();
	if (bmin_use == 1)
	{
		b = ped->GetEllipse().GetBmin();
	}
	//b = 0;
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
	if (bv_use == 1)
	{
		//effdis = dist.Norm() - ped->GetEllipse().GetBmax();
		effdis = dist.Norm() - ped->GetEllipse().GetEB();
	}
	double cosangle = dist.ScalarProduct(ei) / (dist.Norm()*ei.Norm());
	if (cosangle < 0.00001)
	{
		return spacing;
	}
	spacing = effdis / cosangle;
	/*
	double sinangle_sp = dist.CrossProduct(p12) / (dist.Norm()*p12.Norm());
	double sinangle_ei = ei.CrossProduct(p12) / (ei.Norm()*p12.Norm());
	sinangle_sp = sinangle_sp > 0 ? sinangle_sp : -sinangle_sp;
	if (sinangle_sp <= 0.01)
	{
	return spacing;
	}
	sinangle_ei = sinangle_ei > 0 ? sinangle_ei : -sinangle_ei;
	spacing = effdis * sinangle_sp / sinangle_ei;
	*/
	/*
	//test code
	if (ped->GetID() == 17)
	{
	printf("\npedID=%d,position=(%f,%f),pt=(%f,%f),effdis=%f,spacing=%f\n", ped->GetID(), pp._x, pp._y, pt._x, pt._y, effdis, spacing);
	}
	*/
	return spacing;
}

string GCVMModel::GetDescription()
{
	string rueck;
	char tmp[CLENGTH];

	sprintf(tmp, "\t\ta: \t\tPed: %f \tWall: %f\n", _aPed, _aWall);
	rueck.append(tmp);
	sprintf(tmp, "\t\tD: \t\tPed: %f \tWall: %f\n", _DPed, _DWall);
	rueck.append(tmp);
	return rueck;
}

std::shared_ptr<DirectionStrategy> GCVMModel::GetDirection() const
{
	return _direction;
}


double GCVMModel::GetaPed() const
{
	return _aPed;
}

double GCVMModel::GetDPed() const
{
	return _DPed;
}


double GCVMModel::GetaWall() const
{
	return _aWall;
}

double GCVMModel::GetDWall() const
{
	return _DWall;
}
