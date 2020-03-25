/**
* \file        AGCVM.h
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


#ifndef AGCVMMODEL_H_
#define AGCVMMODEL_H_

#include <vector>
#include <math.h>
#include "../geometry/Building.h"
#include "OperationalModel.h"

typedef std::pair<double, double> my_pair;
typedef std::pair<int, int> ID_pair;
typedef std::tuple<int, double> inf_pair;
// sort with respect to first element (ascending).
// In case of equality sort with respect to second element (descending)

struct sort_pred_agcvm
{
	bool operator () (const my_pair& left, const my_pair& right)
	{
		return (left.first == right.first) ?
			(left.second > right.second) :
			(left.first < right.first);
	}
};



//forward declaration
class Pedestrian;
class DirectionStrategy;


class AGCVMModel : public OperationalModel {
private:

	// Modellparameter (CVM)
	double _aPed;
	double _DPed;
	double _aWall;
	double _DWall;

	// GCVM
	double _Ts;
	double _Td;
	int _GCVMUsing=1;// Keep it for incase

	// Clogging
	int _Parallel=1;
	double _WaitingTime = 2;
	int _clogging_times = 0;
	// Boundary case
	double _left_boundary = -100;
	double _right_boundary = 100;
	double _up_boundary = 100;
	double _down_boundary = -100;
	double _cutoff = 2;

	// Switch
	int _Anticipation = 1;
	int _Cooperation = 1;
	int _AttracForce = 1;
	double _AntiT=0;
	double _CoopT = 0;

	// Functions
	double OptimalSpeed(Pedestrian* ped, double spacing) const;
	Point e0(Pedestrian *ped, Room* room) const;
	my_pair GetSpacing(Pedestrian* ped1, Pedestrian* ped2, int periodic, bool collision) const;
	Point ForceRepPed(Pedestrian* ped1, Pedestrian* ped2, Point e0, int periodic, bool push) const;
	Point ForceRepRoom(Pedestrian* ped, SubRoom* subroom, Point e0) const;
	Point ForceRepWall(Pedestrian* ped, const Line& l, const Point& centroid, bool inside, Point e0) const;
	double GetSpacingRoom(Pedestrian* ped, SubRoom* subroom) const;
	double GetSpacingWall(Pedestrian* ped, const Line& l) const;
	void UpdatePed(Pedestrian* ped, Point speed, Point direction, double deltaT, int periodic);
	bool ReArrange(const vector< Pedestrian* >& allPeds_ini, vector< Pedestrian* >& allPeds, Building* building);
	int GetGCVMU() const;
	int GetUpdate() const;
	double GetWaitingTime() const;
	double GetLeftBoundary() const;
	double GetRightBoundary() const;
	double GetUpBoundary() const;
	double GetDownBoundary() const;
	double GetCutoff() const;

	int GetAnticipation() const;
	int GetCooperation() const;
	int GetAttracForce() const;
	double GetAntiT() const;
	double GetCoopT() const;

	my_pair JudgeCollision(Pedestrian* ped1, Pedestrian* ped2) const;
	Point GetInfDirection(Point e0, Point ep12) const;
	Point GetPosPeriodic(Pedestrian* ped1, Pedestrian* ped2) const;
public:

	AGCVMModel(std::shared_ptr<DirectionStrategy> dir, double aped, double Dped,
		double awall, double Dwall, double Ts, double Td, int GCVM, 
		int Parallel, double waitingTime, double lb, double rb, double ub, double db, double co,
		int Anticipation, int Cooperation, int AttracForce, double AntiT, double CoopT);
	virtual ~AGCVMModel(void);


	std::shared_ptr<DirectionStrategy> GetDirection() const;
	double GetaPed() const;
	double GetDPed() const;
	double GetaWall() const;
	double GetDWall() const;
	double GetTs() const;
	double GetTd() const;
	virtual std::string GetDescription();
	virtual bool Init(Building* building);
	virtual void ComputeNextTimeStep(double current, double deltaT, Building* building, int periodic);
};


#endif 
