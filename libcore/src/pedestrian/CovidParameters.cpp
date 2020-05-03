/**
 * \file        CovidParameters.cpp
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

#include "CovidParameters.h"

#include "general/Macros.h"

CovidParameters::CovidParameters(int id, int seed)
{
	_id = id;
	_generator = std::default_random_engine(seed);
}

CovidParameters::~CovidParameters() {}

int CovidParameters::GetID()
{
	return _id;
}

void CovidParameters::SetInfection(int infection)
{
	_infective = infection;
}

void CovidParameters::Setk(double k)
{
	_fitting_k = k;
}

void CovidParameters::SetD(double D)
{
	_fitting_D = D;
}

int CovidParameters::GetInfection()
{
	return _infective;
}

void CovidParameters::InitP(double mean, double stdv)
{
	stdv = stdv == 0 ? _judge : stdv;
	_P = std::normal_distribution<double>(mean, stdv);
}

void CovidParameters::InitQ(double mean, double stdv)
{
	stdv = stdv == 0 ? _judge : stdv;
	_Q = std::normal_distribution<double>(mean, stdv);
}

void CovidParameters::InitAlpha(double mean, double stdv)
{
	stdv = stdv == 0 ? _judge : stdv;
	_alpha = std::normal_distribution<double>(mean, stdv);
}

double CovidParameters::GetP()
{
	return _P.stddev() == _judge ? _P.mean() : _P(_generator);
}

double CovidParameters::GetQ()
{
	return _Q.stddev() == _judge ? _Q.mean() : _Q(_generator);
}

double CovidParameters::GetAlpha()
{
	return _alpha.stddev() == _judge ? _alpha.mean() : _alpha(_generator);
}

double CovidParameters::Getk()
{
	return _fitting_k;
}

double CovidParameters::GetD()
{
	return _fitting_D;
}


