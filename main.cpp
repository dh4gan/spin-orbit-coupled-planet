/*
 * CLIMATES ON SPIN-ORBIT RESONANT PLANETS
 *
 * main.cpp
 *
 * Written by dh4gan 7/8/13
 *
 * Program creates a Star and Planet Object, and measures the flux as a
 * function of longitude and latitude given the spin-orbit coupling
 *
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Star.h"
#include "Planet.h"
#include "Constants.h"

#include <fstream>
#include <sstream>
using namespace std;

int main()
    {

    string par;
    string line;
    string fileString, numString, prefixString;
    char inputFile[100],outputFile[100];

    double mEarthToMSol = 3.0034e-6;

    int nTime, nLongitude, nLatitude, nlambda;
    string starname, planetname;
    double mstar, radstar, lstar, totalMass, G;
    double mplanet, mplanetearth, radplanet;
    double long_apparent, dlat, dlong, rdotn;
    double semimaj, ecc, inc, longascend, argper, time;
    double PsToPo, Teff, Porbit, Pspin, dt;

    FILE * output;

    printf("  \n");
    printf("*********************************************** \n");
    printf("    Flux in Spin-Orbit Resonant Planets \n");
    printf("    Date Created : 7th August 2013 \n");
    printf("*********************************************** \n");
    printf("  \n");

    // Read in parameters file

    fileString = "spinorbit_coupled_planet.params";

    printf("Reading in parameter file");

    strcpy(inputFile, fileString.c_str());

    ifstream myfile(inputFile);
    // check that the file exists as cpp does not check
    //and it just returns 0 if it doesnt exist
    if (myfile == 0)
	{
	cout << "No Input file found" << endl;
	return 0;
	}

    // Then loop through each line using getline and then
    //assign to vectors
    while (getline(myfile, line))
	{
	istringstream iss(line);
	iss >> par;

	if(par =="OutputPrefix"){

	    iss >> prefixString;
	    prefixString = prefixString + '.';
	}

	if (par == "NTime")
	    {
	    iss >> nTime;

	    }

	if (par == "NLatitude")
	    {
	    iss >> nLatitude;

	    }

	if (par == "NLongitude")
	    {
	    iss >> nLongitude;

	    }

	if (par == "PsToPo")
	    {
	    iss >> PsToPo;

	    }

	if (par == "StarMass")
	    {
	    iss >> mstar;

	    }

	if (par == "StarLum")
	    {
	    iss >> lstar;

	    }

	if (par == "StarTeff")
	    {
	    iss >> Teff;
	    }

	if (par == "PlanetMass")
	    {
	    iss >> mplanetearth; // Read in in earth masses
	    mplanet = mplanetearth*mEarthToMSol; // convert to Solar masses

	    }

	if (par == "SemiMajorAxis")
	    {
	    iss >> semimaj;
	    }
	if (par == "Eccentricity")
	    {
	    iss >> ecc;

	    }

	if (par == "Inclination")
	    {
	    iss >> inc;

	    }

	if (par == "AscendingNode")
	    {
	    iss >> longascend;

	    }

	if (par == "Periapsis")
	    {
	    iss >> argper;

	    }

	}
    myfile.close();

    cout << "Read in " << fileString << endl;

    // Number of zeros for output files

    double nzeros = int(log10(nTime)+1);

    // set up Star and Planet Objects


    totalMass = mstar + mplanet;
    G = 4.0*pi*pi;

    Vector3D starpos(0.0, 0.0, 0.0);
    Vector3D starvel(0.0, 0.0, 0.0);

    time = 0.0;

    starname = "Star";
    planetname = "Planet";
    nlambda = 1000;

    Star star(starname, mstar, radstar, starpos, starvel, Teff, nlambda);
    Planet planet(planetname, mplanet, radplanet, semimaj, ecc, inc, time,
	    longascend, argper, G, totalMass);

    cout << "Planet Position, Velocity " << endl;
    cout << planet.getPosition().elements[0] << "  "
	    << planet.getPosition().elements[1] << "  "
	    << planet.getPosition().elements[2] << "  " << endl;
    cout << planet.getVelocity().elements[0] << "  "
	    << planet.getVelocity().elements[1] << "  "
	    << planet.getVelocity().elements[2] << "  " << endl;

    Porbit = sqrt(4.0 * pi * pi * semimaj * semimaj * semimaj / (G * totalMass));
    Pspin = Porbit * PsToPo;

    dt = Porbit / float(nTime);

    // Set up arrays to store fluxes
    double flux[nLongitude][nLatitude];

    // Arrays to store latitude and longitude

    double latitude[nLatitude];
    double longitude[nLongitude];

    dlat = pi/float(nLatitude);
    dlong = 2.0*pi/float(nLongitude);

    for (int j = 0; j < nLongitude; j++)
	{
	longitude[j] = j * dlong;
	}
    for (int k = 0; k < nLatitude; k++)
	{
	latitude[k] = k * dlat;
	}

    // Begin Loop over time

    for (int i = 0; i < nTime; i++)
	{

	// Update time
	time = i * dt;

	// Update planet position
	planet.calcTrueAnomaly(G, totalMass, time);
	planet.calcVectorFromOrbit(G, totalMass);

	Vector3D pos = planet.getPosition().unitVector();

	double magpos = pos.magVector();

	Vector3D unitpos = pos.unitVector();

	// Now begin calculation over latitude and longitude points

	for (int j = 0; j < nLongitude; j++)
	    {
	    // Account for planet rotation
	    long_apparent = longitude[j] + 2.0 * pi * time / Pspin;
	    for (int k = 0; k < nLatitude; k++)

		{

		// construct surface normal vector

		Vector3D surface(sin(latitude[k]) * cos(long_apparent), sin(
			latitude[k]) * sin(long_apparent), cos(latitude[k]));

		surface = surface.unitVector();

		// If necessary, rotate surface vector by obliquity (TODO)
		//surface = surface.rotateX(-obliquity);

		// take the dot product with the unit position vector

		rdotn = unitpos.dotProduct(surface);

		// Calculate fluxes

		if (rdotn > 0.0)
		    {

		    flux[j][k] = lstar * rdotn / (4.0 * pi * magpos * magpos);
		    }
		else
		    {
		    flux[j][k] = 0.0;
		    }

		}
	    }

	//Send the data to file

	// Create filename

	ostringstream convert;
	convert << i+1;

	numString = convert.str();

	while (numString.length() < nzeros)
		{
		numString = "0" + numString;
		}

	fileString = prefixString+numString;

	strcpy(outputFile, fileString.c_str());

	output = fopen(outputFile, "w");
	fprintf(output, "%+.4E %i %i %+.4E  %+.4E  %+.4E %+.4E  %+.4E\n", time, nLatitude,nLongitude, magpos, planet.getTrueAnomaly(), planet.getPosition().elements[0],planet.getPosition().elements[1],planet.getPosition().elements[2]);

	for (int j = 0; j < nLongitude; j++)
	    {
	    for (int k = 0; k < nLatitude; k++)
		{
		fprintf(output, "%+.4E  %+.4E  %+.4E\n", longitude[j],
			latitude[k], flux[j][k]);
		}
	    }
	fflush(output);
	fclose(output);

	}

    // End of program
    return 0;
    }
