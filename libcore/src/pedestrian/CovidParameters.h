/**
 * \file        CovidParameters.h
 * \date        May 3, 2020
 * \version     v0.8
 * \copyright   <2016-2022> Forschungszentrum J¨¹lich GmbH. All rights reserved.
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
 * This class contains the different covid related parameters for the agents.
 * They are used to defined different population group, for instance infective and susceptible.
 *
 **/
#pragma once

#include <random>

class CovidParameters
{
public:
    /**
	 * Constructor
	 */
    CovidParameters(int id, int seed = 1234);

    /**
	 * Destructor
	 */
    virtual ~CovidParameters();

    /**
	 * @return the ID of the covid parameters sets.
	 */
    int GetID();

    /**
	 * Set the status of infection
	 * @param infection, the status of infection
	 */
    void SetInfection(int infection);

    /**
	 * Set the value of parameter k
	 */
    void Setk(double k);

    /**
	 * Set the value of parameter D
	 */
    void SetD(double D);

    /**
	 * Initialize the probability of releasing virus
	 * @param mean, mean value
	 * @param stv, standard deviation
	 */
    void InitP(double mean, double stdv);

    /**
	 * Initialize the probability of infection after the virus entering the body
	 * @param mean, mean value
	 * @param stv, standard deviation
	 */
    void InitQ(double mean, double stdv);

    /**
	 * Initialize the protecting measurement taken by the pedestrian
	 * @param mean, mean value
	 * @param stv, standard deviation
	 */
    void InitAlpha(double mean, double stdv);

    /**
	 * @return the status of infection
	 */
    int GetInfection();

    /**
	 * @return the value of parameter k
	 */
    double Getk();

    /**
	 * @return the value of parameter D
	 */
    double GetD();


    /**
	 * @return a random number following the distribution
	 */
    double GetP();

    /**
	 * @return a random number following the distribution
	 */
    double GetQ();

    /**
	 * @return a random number following the distribution
	 */
    double GetAlpha();

private:
    int _id;
    int _infective    = 0;
    double _fitting_k = 1;
    double _fitting_D = 1;
    std::default_random_engine _generator;
    std::normal_distribution<double> _P;
    std::normal_distribution<double> _Q;
    std::normal_distribution<double> _alpha;
    const double _judge = 10000;
};
