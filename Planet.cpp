/*
 * Planet.cpp
 *
 *  Created on: Nov 8 2012
 *      Author: dh4gan
 */

#include "Planet.h"

Planet::Planet() :
    Body()
    {
    }
Planet::Planet(string &namestring, double &m, double &rad, Vector3D &pos,
	Vector3D &vel) :
    Body(namestring, m, rad, pos, vel)
    {
    }

//

Planet::Planet(string &namestring, double &m, double &rad, double semimaj,
	double ecc, double inc, double trueAnom, double longascend,
	double argper, double G, double totalMass) :
    Body(namestring, m, rad,semimaj, ecc, inc, trueAnom, longascend, argper, G,totalMass)
{
    }

Planet::~Planet()
    {
    }

