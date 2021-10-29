/**
* \file        AVM.h
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


#ifndef AVMMODEL_H_
#define AVMMODEL_H_

#include <vector>
#include <math.h>
#include "../geometry/Building.h"
#include "OperationalModel.h"

typedef std::pair<double, double> my_pair;

struct sort_pred_agcvm
{
    bool operator () (const my_pair& left, const my_pair& right)
    {
        return (left.first == right.first) ?
            (left.second > right.second) :
            (left.first < right.first);
    }
};

class Pedestrian;
class DirectionStrategy;

class AVMModel : public OperationalModel {
private:
    //The model used (0:CSM, 1:GCVM, 2:AVM)
    int _Model;

    //Model Parameters (CVM)
    double _aPed;
    double _DPed;
    double _aWall;
    double _DWall;

    //Model Parameters (GCVM)
    double _Ts;
    double _Td;

    //Model Parameters (AVM)
    double _AntiTime;
    bool _ConstantAlpha;

    //Boundary case
    double _LeftBoundary;
    double _RightBoundary;
    double _UpBoundary;
    double _DownBoundary;
    double _CutOff;

    //Circle experiment
    bool _Circle_Antipode;


    // Functions 
    Point DesireDirection(Pedestrian *ped, Room* room) const;

    Point ForceRepPedGCVM(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const;
    Point ForceRepPedCSM(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const;
    Point ForceRepPedAVM(Pedestrian* ped1, Pedestrian* ped2, Building* building, int periodic) const;

    Point ForceRepRoom(Pedestrian* ped, SubRoom* subroom) const;
    Point ForceRepWall(Pedestrian* ped, const Line& l, const Point& centroid, bool inside) const;

    my_pair GetSpacing(Pedestrian* ped1, Pedestrian* ped2, int periodic) const;
    double GetSpacingRoom(Pedestrian* ped, SubRoom* subroom) const;
    double GetSpacingWall(Pedestrian* ped, const Line& l) const;
    double OptimalSpeed(Pedestrian* ped, double spacing) const;

    Point GetPosPeriodic(Pedestrian* ped1, Pedestrian* ped2) const;//Get the periodic position of ped2 for ped1
    Point GetInfDirection(Point e0, Point ep12) const;
    void UpdatePed(Pedestrian* ped, Point speed, Point direction, double deltaT, int periodic);

public:
    AVMModel(std::shared_ptr<DirectionStrategy> dir, int model,
        double aped, double Dped, double awall, double Dwall,
        double Ts, double Td,
        double AntiT, bool calpha,
        double lb, double rb, double ub, double db, double co,
        bool circle);
    ~AVMModel(void);

    std::string GetDescription() override;
    bool Init(Building* building) override;
    void ComputeNextTimeStep(double current, double deltaT, Building* building, int periodic) override;

    std::shared_ptr<DirectionStrategy> GetDirection() const { return _direction; };
    int GetModel() const { return _Model; }
    double GetaPed() const { return _aPed; };
    double GetDPed() const { return _DPed; };
    double GetaWall() const { return _aWall; };
    double GetDWall() const { return _DWall; };
    double GetTs() const { return _Ts; };
    double GetTd() const { return _Td; };
    double GetAntiT() const { return _AntiTime; };
    bool GetConstantAlpha() const { return _ConstantAlpha; };

    double GetLeftBoundary() const { return _LeftBoundary; };
    double GetRightBoundary() const { return _RightBoundary; };
    double GetUpBoundary() const { return _UpBoundary; };
    double GetDownBoundary() const { return _DownBoundary; };
    double GetCutoff() const { return _CutOff; };
};
#endif 
