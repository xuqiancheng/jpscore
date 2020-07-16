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

	// Model Parameter (CVM)
	double _aPed;
	double _DPed;
	double _aWall;
	double _DWall;

	// GCVM
	double _Ts;
	double _Td;
	int _GCVM;

	// Boundary case
	double _LeftBoundary;
	double _RightBoundary;
	double _UpBoundary;
	double _DownBoundary;
	double _CutOff;

	// Anticipation
	int _Anticipation;
	int _Cooperation;
	int _AttractiveForce;
	int _PushingForce;
	double _AntiTime;
	double _CoopTime;
	double _CoreSize;

	// Functions 
	Point DesireDirection(Pedestrian *ped, Room* room) const;

	Point ForceRepPed(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const;
	Point ForceRepPedPush(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const;
	Point ForceRepRoom(Pedestrian* ped, SubRoom* subroom) const;
	Point ForceRepWall(Pedestrian* ped, const Line& l, const Point& centroid, bool inside) const;

	my_pair GetSpacing(Pedestrian* ped1, Pedestrian* ped2, int periodic) const;
	double GetSpacingRoom(Pedestrian* ped, SubRoom* subroom) const;
	double GetSpacingWall(Pedestrian* ped, const Line& l) const;

	double OptimalSpeed(Pedestrian* ped, double spacing) const;

	// Functions helpful
	Point GetPosPeriodic(Pedestrian* ped1, Pedestrian* ped2) const;//Get the periodic position of ped2 for ped1
	Point GetInfDirection(Point e0, Point ep12) const;
	void UpdatePed(Pedestrian* ped, Point speed, Point direction, double deltaT, int periodic);

	// Function may helpful
	my_pair JudgeCollision(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const;
	bool Drill(Pedestrian* ped, vector<Pedestrian*> neighbours, Building* building, SubRoom* subroom, Point e0, int periodic) const;
	bool DrillRoom(Pedestrian* ped, SubRoom* subroom, Point e0) const;
	bool DrillWall(Pedestrian* ped, Point e0, const Line& l) const;
	Point CorrectD(Pedestrian *ped, Point d_direction, SubRoom* subroom) const;
	Point CorrectDWall(Pedestrian *ped, Point d_direction, const Line& l) const;

	// Function not use not
	bool ReArrange(const vector< Pedestrian* >& allPeds_ini, vector< Pedestrian* >& allPeds, Building* building);



public:
	AGCVMModel(std::shared_ptr<DirectionStrategy> dir,
		double aped, double Dped, double awall, double Dwall,
		double Ts, double Td, int GCVM,
		double lb, double rb, double ub, double db, double co,
		int Anticipation, int Cooperation, int AttracForce, int Push,
		double AntiT, double CoopT, double CoreSize);
	~AGCVMModel(void) override;

	std::string GetDescription() override;
	bool Init(Building* building) override;
	void ComputeNextTimeStep(double current, double deltaT, Building* building, int periodic) override;

	inline std::shared_ptr<DirectionStrategy> GetDirection() const { return _direction; };
	inline double GetaPed() const { return _aPed; };
	inline double GetDPed() const { return _DPed; };
	inline double GetaWall() const { return _aWall; };
	inline double GetDWall() const { return _DWall; };

	inline double GetTs() const { return _Ts; };
	inline double GetTd() const { return _Td; };
	inline int GetGCVMU() const { return _GCVM; };

	inline double GetLeftBoundary() const { return _LeftBoundary; };
	inline double GetRightBoundary() const { return _RightBoundary; };
	inline double GetUpBoundary() const { return _UpBoundary; };
	inline double GetDownBoundary() const { return _DownBoundary; };
	inline double GetCutoff() const { return _CutOff; };

	inline int GetAnticipation() const { return _Anticipation; };
	inline int GetCooperation() const { return _Cooperation; };
	inline int GetAttracForce() const { return _AttractiveForce; };
	inline int GetPushing() const { return _PushingForce; };

	inline double GetAntiT() const { return _AntiTime; };
	inline double GetCoopT() const { return _CoopTime; };
	inline double GetCoreSize() const { return _CoreSize; };

};
#endif 
