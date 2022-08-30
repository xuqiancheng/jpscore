/**
 * \file        Pedestrian.h
 * \date        Sep 30, 2010
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
 *
 *
 **/
#pragma once

#include "Ellipse.hpp"
#include "Journey.hpp"
#include "Line.hpp"
#include "Macros.hpp"
#include "UniqueID.hpp"

#include <memory>

class Pedestrian
{
public:
    using UID = jps::UniqueID<Pedestrian>;

    // This is evaluated by the "strategic level"
    std::unique_ptr<Behaviour> behaviour{};

    // This is evaluated by the "operational level"
    Point destination{};
    Point waypoint{};

private:
    UID _uid{};

    // gcfm specific parameters
    double _mass = 1; // Mass: 1
    double _tau = 0.5; // Reaction time: 0.5
    double _t = 1.0; // OV function

    double _deltaT = 0.01; // step size
    JEllipse _ellipse{}; // the shape of this pedestrian
    Point _v0 = Point(0, 0); // vector V0

    int _newOrientationDelay = 0;

public:
    Pedestrian() = default;
    ~Pedestrian() = default;

    void SetTau(double tau);
    void SetEllipse(const JEllipse& e);

    double GetT() const;
    void SetT(double T);
    void SetDeltaT(double dt);
    void SetPos(const Point& pos);
    void SetV(const Point& v);
    void SetV0Norm(double v0);
    void SetV0(const Point& p) { _v0 = p; }
    void SetSmoothTurning();
    void IncrementOrientationDelay();
    void SetPhiPed();

    UID GetUID() const;
    double GetMass() const;
    double GetTau() const;
    const JEllipse& GetEllipse() const;
    const Point& GetPos() const;
    const Point& GetV() const;
    const Point& GetV0() const;
    Point GetV0(const Point& target) const;
    void InitV0(const Point& target);

    /**
     * the desired speed is the projection of the speed on the horizontal plane.
     * @return the norm of the desired speed.
     */
    double GetV0Norm() const;
};

std::ostream& operator<<(std::ostream& out, const Pedestrian& pedestrian);
