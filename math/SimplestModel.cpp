/**
* \file       SimplestModel.cpp
* \date       Jan 9, 2019
* \version    v0.8
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
* Generalized collision-free velocity model: Qiancheng (6)
*
*
**/

# define NOMINMAX
#include "../pedestrian/Pedestrian.h"
//#include "../routing/DirectionStrategy.h"
#include "../mpi/LCGrid.h"
#include "../geometry/Wall.h"
#include "../geometry/SubRoom.h"

#include "SimplestModel.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads()  1
#endif

double xRight_simplest = 26.0;
double xLeft_simplest = 0.0;
double cutoff_simplest = 2.0;

//--------------------------------
/*
int Vision_area = 1;//Vision_Area=1:open,Vision_Area=0:close,Vision_Area=2:ei_half,Vision=3:e0_half
int Direction_smooth = 1;//Direction_smooth=1:using tau,Direction_smooth=0:not using tau
int Vertical_influence = 1;//Vertical_influence=1:vertical influence,Vertical_influence=0:original influence
int Velocity_influence = 0;//Velocity-influence=1:velocity wll be considered when calculate influence direction
int bf_use = 0;//when calculate direction 1:use b, 0:effective distance
int bv_use = 0;//when calculate velocity 1: use b, 0: effecitve distance
int bmin_use = 1;//width of area in front 1: bmin, 0: b
int real_distance = 0;//when calculate spacing between two pedestrian
*/
using std::vector;
using std::string;

SimplestModel::SimplestModel(std::shared_ptr<DirectionStrategy> dir, double aped, double Dped,
	double awall, double Dwall, double Ts, double Td, int Parallel, double waitingTime, double areasize, int sDirection, int sSpeed, int GCVMU)
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
	_Parallel = Parallel;
	_WaitingTime = waitingTime;
	_SubmodelDirection = sDirection;
	_SubmodelSpeed = sSpeed;
	_GCVMUsing = GCVMU;
	_AreaSize = areasize;
}


SimplestModel::~SimplestModel()
{

}

bool SimplestModel::Init(Building* building)
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
				"ERROR:\tSimplestModel::Init() cannot initialise route. ped %d is deleted in Room %d %d.\n", ped->GetID(), ped->GetRoomID(), ped->GetSubRoomID());
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

void SimplestModel::ComputeNextTimeStep(double current, double deltaT, Building* building, int periodic)
{
	// collect all pedestrians in the simulation.
	const vector< Pedestrian* >& allPeds_ini = building->GetAllPedestrians();
	vector<Pedestrian*> pedsToRemove;
	pedsToRemove.reserve(500);

	unsigned long nSize;
	nSize = allPeds_ini.size();

	//Rearrange pedestrians according to the distance to exit (Using new function ReArrange)
	vector<Pedestrian*> allPeds = vector<Pedestrian*>();
	allPeds.reserve(nSize);
	ReArrange(allPeds_ini, allPeds, building);

	vector< Point > result_acc = vector<Point >();
	result_acc.reserve(nSize);

	vector< Point > result_dir = vector<Point >();
	result_dir.reserve(nSize);

	vector< my_pair > spacings = vector<my_pair >();
	spacings.reserve(nSize); // larger than needed

	vector< ID_pair > relations = vector<ID_pair>();
	relations.reserve(nSize);

	int start = 0;
	int end = nSize - 1;
	for (int p = start; p <= end; ++p) {
		Pedestrian* ped = allPeds[p];
		Room* room = building->GetRoom(ped->GetRoomID());
		SubRoom* subroom = room->GetSubRoom(ped->GetSubRoomID());
		Point repPed = Point(0, 0);
		vector<Pedestrian*> neighbours;
		building->GetGrid()->GetNeighbourhood(ped, neighbours);

		// Get desired moving direction
		Point inid_direction = e0(ped, room);

		int size = (int)neighbours.size();
		Point direction;

		int UDirection = GetSDirection(); // If using direction part?
		if (!UDirection)
		{
			direction = inid_direction;
		}
		else
		{
			//using a new method calculate the influence of pedestrian (the value of influence id decided by distance and the direction is vertical with desired direction)
			Point inf_direction;// left side of pedestrian

			// If using GCVM or CVM?
			if (GetGCVMU())
			{
				inf_direction._x = -inid_direction._y;
				inf_direction._y = inid_direction._x;
				inf_direction = inf_direction.Normalized();
			}
			//--------------------------------------------------------------------------------------------------------------------------------------------------------------

			//Calculating influence of pedestrians---------------------------------------------------------------------
			for (int i = 0; i < size; i++) {
				Pedestrian* ped1 = neighbours[i];
				//if they are in the same subroom
				Point p1 = ped->GetPos();
				Point p2 = ped1->GetPos();
				Point ep12 = p2 - p1;

				//Deciding the influence direction--------------------------------------------------------------------
				if (GetGCVMU())
				{
					double result1 = inid_direction.CrossProduct(ep12);
					Point zero = Point(0, 0);
					if (bool equal = almostEqual(result1, 0, 0.00001))
					{
						int random = rand() % 1000;//choose one direciton by random
						if (random < 500)
						{
							inf_direction = zero - inf_direction;
						}
					}
					else
					{
						double result2 = inid_direction.CrossProduct(inf_direction);
						if (result1*result2 > 0)
						{
							inf_direction = zero - inf_direction;
						}
					}
				}

				//----------------------------------------------------------------------------------------------------
				//subrooms to consider when looking for neighbour for the 3d visibility
				vector<SubRoom*> emptyVector;
				emptyVector.push_back(subroom);
				emptyVector.push_back(building->GetRoom(ped1->GetRoomID())->GetSubRoom(ped1->GetSubRoomID()));

				//This part is not my work, may be have to read it. But in my simple case, may be no influence.
				bool isVisible = building->IsVisible(p1, p2, emptyVector, false);
				if (!isVisible) {
					continue;
				}
				//-------------------------------------------------------

				if (ped->GetUniqueRoomID() == ped1->GetUniqueRoomID()) {
					if (GetGCVMU())
					{
						repPed += inf_direction * ForceRepPed(ped, ped1, inid_direction, periodic).Norm();//GCVM
					}
					else
					{
						repPed += ForceRepPed(ped, ped1, inid_direction, periodic);//CVM
					}
				}
				else {
					SubRoom* sb2 = building->GetRoom(ped1->GetRoomID())->GetSubRoom(ped1->GetSubRoomID());
					if (subroom->IsDirectlyConnectedWith(sb2)) {
						if (GetGCVMU())
						{
							repPed += inf_direction * ForceRepPed(ped, ped1, inid_direction, periodic).Norm();//GCVM
						}
						else
						{
							repPed += ForceRepPed(ped, ped1, inid_direction, periodic);//CVM
						}
					}
				}
			}
			//-----------------------------------------------------------------------------------------------------------

			// Calculating influence of walls-------------------------------------------------------------------------
			Point repWall = ForceRepRoom(ped, subroom, inid_direction);
			if (ped->GetExitLine()->DistTo(ped->GetPos()) < 0.2)
			{
				std::vector<SubRoom*> Nsubrooms = subroom->GetNeighbors();
				for (int i = 0; i < Nsubrooms.size(); i++)
					repWall += ForceRepRoom(ped, Nsubrooms[i], inid_direction);
			}
			//-----------------------------------------------------------------------------------------------------------

			//Caluculating desired direcition----------------------------------------------------------------------------------------------
			Point d_direction;
			d_direction = inid_direction + repPed + repWall;
			//------------------------------------------------------------------------------------------------------------------------------
			direction = d_direction;

			if (GetGCVMU())
			{
				//Calculating the actual direction of pedestrian at next timestep-----
				Point a_direction;
				a_direction._x = ped->GetEllipse().GetCosPhi();
				a_direction._y = ped->GetEllipse().GetSinPhi();
				double angle_tau = _Td;
				Point angle_v = (d_direction.Normalized() - a_direction) / angle_tau;
				direction = a_direction + angle_v * deltaT;

			}
		}
		//Update direction
		direction = direction.Normalized();
		result_dir.push_back(direction);
		//------------------------------------------------------------------------------------------------------------------------------------


		//Calculating spacing in front -------------------------------------------------------------------------------------------------------
		for (int i = 0; i < size; i++) {
			Pedestrian* ped1 = neighbours[i];
			// calculate spacing
			if (ped->GetUniqueRoomID() == ped1->GetUniqueRoomID()) {
				spacings.push_back(GetSpacing(ped, ped1, direction, periodic));
			}
			else {
				// or in neighbour subrooms
				SubRoom* sb2 = building->GetRoom(ped1->GetRoomID())->GetSubRoom(ped1->GetSubRoomID());
				if (subroom->IsDirectlyConnectedWith(sb2)) {
					spacings.push_back(GetSpacing(ped, ped1, direction, periodic));
				}
			}
		}
		//------------------------------------------------------------------------------------------------------------------------------------

		// Calculate min spacing, and save all spacing=0 to vector relation--------------------------------------------------------------------
		std::sort(spacings.begin(), spacings.end(), sort_pred_Simplest());
		double spacing = spacings.size() == 0 ? 100 : spacings[0].first;
		double first_ID = spacings.size() == 0 ? -1 : spacings[0].second;

		// add this part to avoid pedestrian cross the wall directly
		// some pedestrian are blocked by wall
		double spacing_wall = GetSpacingRoom(ped, subroom, direction);
		spacing = spacing > spacing_wall ? spacing_wall : spacing;

		//optimal speed function
		Point speed;
		speed = direction.NormalizedMolified() *OptimalSpeed(ped, spacing);
		result_acc.push_back(speed);
		spacings.clear(); //clear for ped p
		//unparallel update
		if (_Parallel != 1) {
			Point pos_neu = ped->GetPos() + speed * deltaT;
			//calculate ellipse orientation
			double normdir = direction.Norm();
			double cosPhi = direction._x / normdir;
			double sinPhi = direction._y / normdir;
			JEllipse e = ped->GetEllipse();
			e.SetCosPhi(cosPhi);
			e.SetSinPhi(sinPhi);
			ped->SetEllipse(e);
			ped->SetV(speed);

			if (periodic) {
				if ((ped->GetPos()._x < xRight_simplest) && (pos_neu._x >= xRight_simplest)) {
					ped->SetPos(Point(pos_neu._x - (xRight_simplest - xLeft_simplest), pos_neu._y));
				}
				else if ((ped->GetPos()._x > xLeft_simplest) && (pos_neu._x <= xLeft_simplest)) {
					ped->SetPos(Point(pos_neu._x + (xRight_simplest - xLeft_simplest), pos_neu._y));
				}
				else {
					ped->SetPos(pos_neu);
				}
			}
			else {
				ped->SetPos(pos_neu);
			}
		}

	}
	// parallel update
	if (_Parallel == 1) {
		for (int p = start; p <= end; ++p) {
			Pedestrian* ped = allPeds[p];
			Point v_neu = result_acc[p - start];
			Point dir_neu = result_dir[p - start];
			Point pos_neu = ped->GetPos() + v_neu * deltaT;
			//calculate ellipse orientation
			double normdir = dir_neu.Norm();
			double cosPhi = dir_neu._x / normdir;
			double sinPhi = dir_neu._y / normdir;
			JEllipse e = ped->GetEllipse();
			e.SetCosPhi(cosPhi);
			e.SetSinPhi(sinPhi);
			ped->SetEllipse(e);
			ped->SetV(v_neu);

			if (periodic) {
				if ((ped->GetPos()._x < xRight_simplest) && (pos_neu._x >= xRight_simplest)) {
					ped->SetPos(Point(pos_neu._x - (xRight_simplest - xLeft_simplest), pos_neu._y));
				}
				else if ((ped->GetPos()._x > xLeft_simplest) && (pos_neu._x <= xLeft_simplest)) {
					ped->SetPos(Point(pos_neu._x + (xRight_simplest - xLeft_simplest), pos_neu._y));
				}
				else {
					ped->SetPos(pos_neu);
				}
			}
			else {
				ped->SetPos(pos_neu);
			}
		}
	}

	// Save all the couples in clogging into relations and update blockage time
	int blockage = 1;
	for (int p = start; p <= end; ++p)
	{
		Pedestrian* ped1 = allPeds[p];
		Room* room = building->GetRoom(ped1->GetRoomID());
		SubRoom* subroom = room->GetSubRoom(ped1->GetSubRoomID());
		vector<Pedestrian*> neighbours;
		building->GetGrid()->GetNeighbourhood(ped1, neighbours);
		int size = (int)neighbours.size();
		int updateTime = 0;
		for (int i = 0; i < size; i++)
		{
			Pedestrian* ped2 = neighbours[i];
			bool InClogging = NewClogging(ped1, ped2);
			if (InClogging)
			{
				relations.push_back(ID_pair(ped1->GetID(), ped2->GetID()));
				updateTime = 1;
			}
		}
		if (updateTime)
		{
			double InCloggingTime = ped1->GetInCloggingTime() + deltaT;
			ped1->SetInCloggingTime(InCloggingTime);
		}
		else
		{
			ped1->SetInCloggingTime(0);
		}

		//if blockage happens
		double xpos = ped1->GetPos()._x;
		double xpos_before = xpos - ped1->GetV()._x*deltaT;
		double xlimit = 10;
		if (xpos > xlimit&&xpos_before < xlimit)
		{
			blockage = 0;
		}
	}
	if (blockage)
	{
		_blocktime = _blocktime + deltaT;
		_cool_clock = _cool_clock + deltaT;
	}
	else
	{
		_cool_clock = 0;
		_blocktime = 0;
	}

	//if clogging should be solved, then find the closest stable clogging and remove one pedestrian from it.
	if (_cool_clock >= _WaitingTime)
	{
		// find the closest pedestrian who satisfy
		int ID_closest = -1;
		double dis_closest = FLT_MAX;
		for (int p = start; p <= end; ++p)
		{
			Pedestrian* ped = allPeds[p];
			int ID = ped->GetID();
			for (vector<ID_pair>::iterator iter = relations.begin(); iter < relations.end(); ++iter)
			{
				if (iter->first == ID || iter->second == ID)
				{
					double distGoal = (ped->GetExitLine()->GetCentre() - ped->GetPos()).Norm();
					if (ped->GetPos()._x > 10)
					{
						distGoal = 10 - ped->GetPos()._x;
					}
					if (distGoal < dis_closest && ped->GetPos()._x >= 0)
					{
						ID_closest = ID;
						dis_closest = distGoal;
					}
					break;
				}
			}
		}

		//remove it
		for (int p = start; p <= end; ++p)
		{
			Pedestrian* ped = allPeds[p];
			if (ped->GetID() == ID_closest)
			{
				//printf("\nTest:ID=%d time=%f", ID_closest, current);
				if (_blocktime < 2 * _WaitingTime && abs(_blocktime - current)>_WaitingTime)
				{
					_clogging_times++;
				}
				std::ofstream ofile;
				string ProjectFileName = building->GetProjectFilename();
				int startN = ProjectFileName.find_last_of("\\");
				startN = startN == -1 ? ProjectFileName.find_last_of("/") : startN;
				int endN = ProjectFileName.find(".xml");
				string InifileName = ProjectFileName.substr(startN + 1, endN - startN - 1);
				if (_clogging_times == 1) {
					ofile.open(building->GetProjectRootDir() + "CloggingLog_" + InifileName + ".txt", std::ofstream::trunc);
					ofile << "#inifile: " << building->GetProjectFilename() << "\n";
					ofile << "#Commit date: " << GIT_COMMIT_DATE << "\n";
					ofile << "#Branch: " << GIT_BRANCH << "\n";
					ofile << "#Timestep: " << deltaT << " (s)\n";
					ofile << "#Waiting time: " << _WaitingTime << " (s)\n";
					ofile << "#Parallel: " << _Parallel << " (1:parallel,0:unparallel)\n";
					ofile << "#Direction: " << _SubmodelDirection << " (1:Using direction submodel,0:Not using direction submodel)\n";
					ofile << "#Speed: " << _SubmodelSpeed << " (1:Using speed submodel,0:Not using speed submodel)\n";
					ofile << "#GCVM: " << _GCVMUsing << " (1:Using GCVM instead of simplest model,0:Using simplest model)\n";
					ofile << "#ID\ttime(s)\tamount\tposition_x\tpostion_y\n";
				}
				else {
					ofile.open(building->GetProjectRootDir() + "CloggingLog_" + InifileName + ".txt", std::ofstream::app);
				}
				ofile << ped->GetID() << "\t" << current << "\t" << _clogging_times << "\t" << ped->GetPos()._x << "\t" << ped->GetPos()._y << "\n";
				ofile.close();
				Point position = ped->GetPos();
				double position_wx = position._x > -8 ? position._x - 18 : position._x - 2;
				Point position_w(position_wx, position._y);
				ped->SetPos(position_w, true);
				ped->SetmoveManually(true);

				for (int i = start; i <= end; ++i)
				{
					Pedestrian* ped1 = allPeds[i];
					ped1->SetInCloggingTime(0);
				}

				_cool_clock = 0;
				break;
			}
		}
	}
	// remove the pedestrians that have left the building
	for (unsigned int p = 0; p < pedsToRemove.size(); p++) {
		building->DeletePedestrian(pedsToRemove[p]);
	}
	pedsToRemove.clear();
}

Point SimplestModel::e0(Pedestrian* ped, Room* room) const
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


double SimplestModel::OptimalSpeed(Pedestrian* ped, double spacing) const
{
	double v0 = ped->GetV0Norm();
	int USpeed = GetSSpeed();
	double speed = 0;
	if (USpeed == 1)
	{
		//optimal speed function
		double T = _Ts;
		//spacing is approximate here
		speed = (spacing) / T;
		speed = (speed > 0) ? speed : 0;
		speed = (speed < v0) ? speed : v0;
	}
	else
	{
		//simplest speed function
		speed = spacing > 0 ? v0 : 0;
		//speed = spacing < v0 ? spacing : v0;
	}
	return speed;
}

// return spacing and id of the nearest pedestrian
my_pair SimplestModel::GetSpacing(Pedestrian* ped1, Pedestrian* ped2, Point ei, int periodic) const
{
	//-------------------------------------------------------------
	double x1 = ped1->GetPos()._x;
	double x2_real = ped2->GetPos()._x;
	double y2 = ped2->GetPos()._y;
	Point ped2_current = ped2->GetPos();

	if (periodic) {
		if ((xRight_simplest - x1) + (x2_real - xLeft_simplest) <= cutoff_simplest) {
			double x2_periodic = x2_real + xRight_simplest - xLeft_simplest;
			ped2->SetPos(Point(x2_periodic, y2));
		}
		if ((x1 - xLeft_simplest) + (xRight_simplest - x2_real) <= cutoff_simplest) {
			double x2_periodic = xLeft_simplest - xRight_simplest + x2_real;
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
		Log->Write("WARNING: \tin SimplestModel::GetSPacing() ep12 can not be calculated!!!\n");
		Log->Write("\t\t Pedestrians are too near to each other (%f).", Distance);
		//exit(EXIT_FAILURE);
	}
	//calculate effective distance
	JEllipse eped1 = ped1->GetEllipse();
	JEllipse eped2 = ped2->GetEllipse();
	double dist;
	double eff_dist = eped1.EffectiveDistanceToEllipse(eped2, &dist);
	double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
	if (condition1 <= 0) {
		ped2->SetPos(ped2_current);
		return  my_pair(FLT_MAX, ped2->GetID());
	}
	if (!ped1->GetEllipse().DoesStretch())
	{
		double l = ped1->GetLargerAxis() + ped2->GetLargerAxis();
		double condition2 = ei.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
		condition2 = (condition2 > 0) ? condition2 : -condition2; // abs
		ped2->SetPos(ped2_current);
		if ((condition1 >= 0) && (condition2 <= l / Distance))
			// return a pair <dist, condition1>. Then take the smallest dist. In case of equality the biggest condition1
			return  my_pair(distp12.Norm() - l, ped2->GetID());
		else
			return  my_pair(FLT_MAX, ped2->GetID());
	}
	else
	{
		//Judge conllision
		//Obtain parameters
		double a1 = ped1->GetLargerAxis();
		Point v = ped1->GetV();
		//double b1 = ped1->GetSmallerAxis();
		double b1 = ped1->GetEllipse().GetBmin();
		/*
		//Avoid block (drill,small tricks)
		if (fabs(v._x) < J_EPS && fabs(v._y) < J_EPS) // v==0
		{
			const Point& pos = ped1->GetPos();
			double distGoal = ped1->GetExitLine()->DistToSquare(pos);
			if (distGoal < 0.2)
			{
				b1 = 0.1;
				//b1 = ped1->GetEllipse().GetBmin();
			}
		}
		*/
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
		//Judge if the line contact with ellipse2
		ped2->SetPos(ped2_current);
		if (Dis1 >= 0 && Dis2 >= 0)
			return  my_pair(FLT_MAX, ped2->GetID());
		else
		{
			return  my_pair(eff_dist, ped2->GetID());
		}
	}

}
Point SimplestModel::ForceRepPed(Pedestrian* ped1, Pedestrian* ped2, Point e0, int periodic) const
{
	Point F_rep(0.0, 0.0);
	// x- and y-coordinate of the distance between p1 and p2

	double x_j = ped2->GetPos()._x;
	double y_j = ped2->GetPos()._y;

	if (periodic) {
		double x = ped1->GetPos()._x;
		if ((xRight_simplest - x) + (x_j - xLeft_simplest) <= cutoff_simplest) {
			double x2_periodic = x_j + xRight_simplest - xLeft_simplest;
			ped2->SetPos(Point(x2_periodic, y_j));
		}
		if ((x - xLeft_simplest) + (xRight_simplest - x_j) <= cutoff_simplest) {
			double x2_periodic = xLeft_simplest - xRight_simplest + x_j;
			ped2->SetPos(Point(x2_periodic, y_j));
		}
	}

	Point distp12 = ped2->GetPos() - ped1->GetPos();
	double Distance = distp12.Norm();
	Point ep12; // x- and y-coordinate of the normalized vector between p1 and p2
	double R_ij;

	//This part is not my work, may be read it later.
	if (Distance >= J_EPS) {
		ep12 = distp12.Normalized();
	}
	else {
		printf("ERROR: \tin VelocityModel::forcePedPed() ep12 can not be calculated!!!\n");
		Log->Write(KRED "\nWARNING: \tin SimplestModel::forcePedPed() ep12 can not be calculated!!!" RESET);
		Log->Write("\t\t Pedestrians are too near to each other (dist=%f).", Distance);
		Log->Write("\t\t Maybe the value of <a> in force_ped should be increased. Going to exit.\n");
		printf("ped1 %d  ped2 %d\n", ped1->GetID(), ped2->GetID());
		printf("ped1 at (%f, %f), ped2 at (%f, %f)\n", ped1->GetPos()._x, ped1->GetPos()._y, ped2->GetPos()._x, ped2->GetPos()._y);
		//exit(EXIT_FAILURE);
	}
	//--------------------------------------------------

	Point ei;
	JEllipse Eped1 = ped1->GetEllipse();
	ei._x = Eped1.GetCosPhi();
	ei._y = Eped1.GetSinPhi();

	//Vision area-----------------------------------
	double condition1 = e0.ScalarProduct(ep12); // < e_i , e_ij > should be positive
	double condition2 = ei.ScalarProduct(ep12);
	if (!GetGCVMU())
	{
		condition1 = 1;
	}
	//-----------------------------------------------
	//rule:pedestrian's direction only influenced by pedestrian in version area
	if (condition1 > 0 || condition2 > 0)
	{
		double dist;
		JEllipse Eped2 = ped2->GetEllipse();
		dist = Eped1.EffectiveDistanceToEllipse(Eped2, &dist);
		// dist is approximate here
		dist = dist > 0 ? dist : 0;
		R_ij = -_aPed * exp((-dist) / _DPed);
		F_rep = ep12 * R_ij;
	}
	ped2->SetPos(Point(x_j, y_j));
	return F_rep;
}//END Velocity:ForceRepPed()

Point SimplestModel::ForceRepRoom(Pedestrian* ped, SubRoom* subroom, Point e0) const
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
		int uid1 = goal->GetUniqueID();
		int uid2 = ped->GetExitIndex();
		// ignore my transition consider closed doors
		//closed doors are considered as wall
		//door is open, but it's not my door (has influence)
		/*
		if((uid1 != uid2) && (goal->IsOpen()==true ))
		{
			f +=  ForceRepWall(ped,*(static_cast<Line*>(goal)), centroid, inside, e0);
		}
		*/

	}
	return f;
}

Point SimplestModel::ForceRepWall(Pedestrian* ped, const Line& w, const Point& centroid, bool inside, Point e0) const
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
	if (GetGCVMU())
	{
		//version area----------------------------------
		double result_e01 = e0.ScalarProduct(e_iw1);
		double result_e02 = e0.ScalarProduct(e_iw2);
		double result_ei1 = ei.ScalarProduct(e_iw1);
		double result_ei2 = ei.ScalarProduct(e_iw2);
		//----------------------------------------------
		if (result_e01 < 0 && result_e02 < 0 && result_ei1 < 0 && result_ei1 < 0)
		{
			return F_wrep;
		}
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
	if (bool equal = almostEqual(result1, 0, 0.00001))//Is there any possible that desired direction towards the wall?
	{
		//double alpha = ped->GetAlpha();
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
	// effdist is approximate here
	effdis = effdis > 0 ? effdis : 0;
	R_iw = -_aWall * exp((-effdis) / _DWall);
	if (GetGCVMU())
	{
		F_wrep = inf_direction * R_iw;//GCVM
	}
	else
	{
		F_wrep = e_iw * R_iw;//Simplest model
	}
	return F_wrep;
}

double SimplestModel::GetSpacingRoom(Pedestrian* ped, SubRoom* subroom, Point ei) const
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
		/*
		if ((uid1 != uid2) && (goal->IsOpen() == true))
		{
			distance = GetSpacingWall(ped, *(static_cast<Line*>(goal)), ei);
			spacing = spacing > distance ? distance : spacing;
		}
		*/
	}
	return spacing;

}

double SimplestModel::GetSpacingWall(Pedestrian* ped, const Line& l, Point ei) const
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
	//double b = ped->GetSmallerAxis();
	double b = ped->GetEllipse().GetBmin();
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

string SimplestModel::GetDescription()
{
	string rueck;
	char tmp[CLENGTH];

	sprintf(tmp, "\t\ta: \t\tPed: %f \tWall: %f\n", _aPed, _aWall);
	rueck.append(tmp);
	sprintf(tmp, "\t\tD: \t\tPed: %f \tWall: %f\n", _DPed, _DWall);
	rueck.append(tmp);
	return rueck;
}

std::shared_ptr<DirectionStrategy> SimplestModel::GetDirection() const
{
	return _direction;
}


double SimplestModel::GetaPed() const
{
	return _aPed;
}

double SimplestModel::GetDPed() const
{
	return _DPed;
}


double SimplestModel::GetaWall() const
{
	return _aWall;
}

double SimplestModel::GetDWall() const
{
	return _DWall;
}

int SimplestModel::GetUpdate() const
{
	return _Parallel;
}

double SimplestModel::GetWaitingTime() const
{
	return _WaitingTime;
}

int SimplestModel::GetSDirection() const
{
	return _SubmodelDirection;
}

int SimplestModel::GetSSpeed() const
{
	return _SubmodelSpeed;
}

int SimplestModel::GetGCVMU() const
{
	return _GCVMUsing;
}

double SimplestModel::GetAreaSize() const
{
	return _AreaSize;
}

bool SimplestModel::ReArrange(const vector< Pedestrian* >& allPeds_ini, vector< Pedestrian* >& allPeds, Building* building)
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


bool SimplestModel::RealClogging(Pedestrian* ped1, Pedestrian* ped2, Point ei, int periodic) const
{
	//-------------------------------------------------------------
	double limitation = 100;
	double x1 = ped1->GetPos()._x;
	double x2_real = ped2->GetPos()._x;
	double y2 = ped2->GetPos()._y;
	Point ped2_current = ped2->GetPos();

	if (periodic) {
		if ((xRight_simplest - x1) + (x2_real - xLeft_simplest) <= cutoff_simplest) {
			double x2_periodic = x2_real + xRight_simplest - xLeft_simplest;
			ped2->SetPos(Point(x2_periodic, y2));
		}
		if ((x1 - xLeft_simplest) + (xRight_simplest - x2_real) <= cutoff_simplest) {
			double x2_periodic = xLeft_simplest - xRight_simplest + x2_real;
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
		Log->Write("WARNING: \tin SimplestModel::GetSPacing() ep12 can not be calculated!!!\n");
		Log->Write("\t\t Pedestrians are too near to each other (%f).", Distance);
		//exit(EXIT_FAILURE);
	}
	//calculate effective distance
	JEllipse eped1 = ped1->GetEllipse();
	JEllipse eped2 = ped2->GetEllipse();
	double dist;
	double eff_dist = eped1.EffectiveDistanceToEllipse(eped2, &dist);
	double condition1 = ei.ScalarProduct(ep12); // < e_i , e_ij > should be positive
	if (condition1 <= 0) {
		ped2->SetPos(ped2_current);
		return true;
	}
	if (!ped1->GetEllipse().DoesStretch())
	{
		double r = ped2->GetLargerAxis();
		double l = ped1->GetLargerAxis() + ped2->GetLargerAxis();
		double condition2 = ei.Rotate(0, 1).ScalarProduct(ep12); // theta = pi/2. condition2 should <= than l/Distance
		condition2 = (condition2 > 0) ? condition2 : -condition2; // abs
		ped2->SetPos(ped2_current);
		if ((condition1 >= 0) && (condition2 <= r / Distance) && (distp12.Norm() - l > 0.01) && (distp12.Norm() - l < limitation))
			// return a pair <dist, condition1>. Then take the smallest dist. In case of equality the biggest condition1
			return  false;
		else
			return  true;
	}
	else
	{
		//Judge conllision
		//Obtain parameters
		double a1 = ped1->GetLargerAxis();
		Point v = ped1->GetV();
		//double b1 = ped1->GetSmallerAxis();
		double b1 = 0.001;
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
			if (eff_dist<0.01 || eff_dist>limitation)
			{
				return true;
			}
			return false;
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
		//Judge if the line contact with ellipse2
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

bool SimplestModel::NewClogging(Pedestrian* ped1, Pedestrian* ped2)
{
	JEllipse eped1 = ped1->GetEllipse();
	JEllipse eped2 = ped2->GetEllipse();
	double dist;
	double eff_dist = eped1.EffectiveDistanceToEllipse(eped2, &dist);

	Point e1;
	e1._x = eped1.GetCosPhi();
	e1._y = eped1.GetSinPhi();
	Point p1 = ped1->GetPos();

	Point e2;
	e2._x = eped2.GetCosPhi();
	e2._y = eped2.GetSinPhi();
	Point p2 = ped2->GetPos();
	Point ep12 = p2 - p1;

	//condition1 : s<0
	double d12 = _AreaSize;
	int condition1 = eff_dist < d12 ? 1 : 0;

	//condition2
	double res1 = e1.ScalarProduct(ep12);
	double res2 = e2.ScalarProduct(ep12);

	//condition3
	Point speed1 = ped1->GetV();
	Point speed2 = ped2->GetV();
	double v01 = ped1->GetV0Norm();
	double v02 = ped2->GetV0Norm();
	double epsi = (v01 + v02) / 100;
	int condition3 = (speed1.Norm() + speed2.Norm()) < epsi ? 1 : 0;

	bool clogging = false;
	if (condition1 == 1 && condition3 == 1 && res1 > 0 && res2 < 0)
	{
		clogging = true;
	}
	return clogging;
}
