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


// Reads in the parameters of the input file
int readParameters(ifstream &inputStream, string &prefixString, int &nTime, int &nLatitude, int &nLongitude, double &PsToPo, 
		   double &mstar, double &lstar, double &Teff, double &mplanet, double &albedo, double &obliquity, double &semi_maj, double &ecc,
		   double &inc, double &longascend, double &argper);

// Calculates the flux between star and planet
void calcFlux(Star &star, Planet &planet, double &longitude, double &latitude, double &hourAngle, 
	      double &flux, double &altitude, double &azimuth,
	      double &time, double &Pspin, double &obliquity);

double safeAcos(double x); // Simple function to stop acos becoming infinite if abs(x) > 1.0

int main()
{

  int success;    
  string fileString, numString, prefixString;
  char inputFile[100], outputFile[100];
      
  int nTime, nLongitude, nLatitude, nlambda;
  string starname, planetname;
  double mstar, radstar, lstar, totalMass, G;
  double mplanet, radplanet, obliquity;
  double dlat, dlong, albedo;
  double semimaj, ecc, inc, longascend, argper, time;
  double PsToPo, Teff, Porbit, Pspin, dt;

  double hplanck = 6.626e-34;
  double c = 2.9e8;
  double stefan = 5.67e-8;

  FILE * output, *outputlog;

  printf("  \n");
  printf("*********************************************** \n");
  printf("    Flux in Spin-Orbit Resonant Planets \n");
  printf("    Date Created : 7th August 2013 \n");
  printf("*********************************************** \n");
  printf("  \n");

  // Read in parameters file

  fileString = "spinorbit_coupled_planet.params";

  printf("Reading in parameter file \n");

  strcpy(inputFile, fileString.c_str());

  // Call Read in function

  ifstream myfile(inputFile);
  success = readParameters(myfile, prefixString, nTime, nLatitude, nLongitude, PsToPo, 
			   mstar, lstar, Teff, mplanet, albedo, obliquity, semimaj, ecc,
			   inc, longascend, argper);
  if(success==-1) {
    return 0;
  }

  cout << "Read in " << fileString << endl;

  // Number of zeros for output files

  double nzeros = int(log10(nTime) + 1);

  // set up Star and Planet Objects

  totalMass = mstar + mplanet;
  G = 4.0 * pi * pi;

  Vector3D starpos(0.0, 0.0, 0.0);
  Vector3D starvel(0.0, 0.0, 0.0);

  time = 0.0;

  starname = "Star";
  planetname = "Planet";
  nlambda = 1000;

  Star star(starname, mstar, radstar, starpos, starvel, Teff, nlambda);
  Planet planet(planetname, mplanet, radplanet, semimaj, ecc, inc, time,
		longascend, argper, G, totalMass);

  star.setLuminosity(lstar);

  // Set up time parameters

  Porbit
    = sqrt(4.0 * pi * pi * semimaj * semimaj * semimaj
	   / (G * totalMass));
  Pspin = Porbit * PsToPo;

  dt = Porbit / float(nTime);

  // Set up array to store fluxes, and altitude and azimuth of star

  double flux[nLongitude][nLatitude];
  double darkness[nLongitude][nLatitude];
  double planetTeff[nLongitude][nLatitude];
  double Ngamma[nLongitude][nLatitude];
  double altitude[nLongitude][nLatitude];
  double azimuth[nLongitude][nLatitude];

  // Set up arrays to store latitude and longitude

  double latitude[nLatitude];
  double longitude[nLongitude];
  double hourAngle[nLongitude];

  dlat = pi / float(nLatitude);
  dlong = 2.0 * pi / float(nLongitude);

  for (int j = 0; j < nLongitude; j++)
    {
      longitude[j] = j * dlong;
      for(int k=0;k<nLatitude; k++)
	  {
	  darkness[j][k]=0.0;
	  }
    }
  for (int k = 0; k < nLatitude; k++)
    {
      latitude[k] = k * dlat;
    }

  // Set up variables for calculating photon fluxes

  double lambda_peak = 2.9e-3/star.getTeff(); // Peak wavelength in SI (Wien's Law)
  double micromol = 1.0e6/6.02214e23; // multiply with to find per micromol
  double Ephot = lambda_peak/(hplanck*c); // multiply with to divide by photon energy

  // Set up output log file

  fileString = prefixString + "log";

  strcpy(outputFile, fileString.c_str());

  outputlog = fopen(outputFile, "w");

  // Begin Loop over time

  for (int i = 0; i < nTime; i++)
    {

      // Update time
      time = i * dt;

      // Update planet position
      planet.calcTrueAnomaly(G, totalMass, time);
      planet.calcVectorFromOrbit(G, totalMass);

      // Loop over longitude and latitude

      for (int j=0; j<nLongitude; j++)
	{
	  for (int k=0; k< nLatitude; k++)
	    {
	      flux[j][k] = 0.0; 

	      // Call flux calculation function
	
	      calcFlux(star,planet,longitude[j], latitude[k],hourAngle[j], 
		       flux[j][k], altitude[j][k], azimuth[j][k],
		       time, Pspin, obliquity);

	      if(flux[j][k]==0.0)
		  {
		  darkness[j][k] = darkness[j][k]+dt;
		  }

	      flux[j][k] = flux[j][k]*lsol/(AU*AU);

	      // Calculate equilibrium temperature (assuming a fixed albedo)

	      planetTeff[j][k] = pow(flux[j][k]*(1.0-albedo)/stefan, 0.25);

	      // Calculate Maximum Photon Flux Density (in micro mol m^-2 s^-1)

	      Ngamma[j][k] = flux[j][k]*Ephot*micromol;

	    }
	}

      //Send the data to file

      // First, write data to log file

      fprintf(outputlog,
	      "%+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E \n",
	      time, Porbit, Pspin, (planet.getPosition()).magVector(), planet.getTrueAnomaly(),
	      planet.getPosition().elements[0],
	      planet.getPosition().elements[1],
	      planet.getPosition().elements[2]);
      fflush(outputlog);

      // Now, write snapshot of latitude/longitude to file

      // Create filename

      ostringstream convert;
      convert << i + 1;

      numString = convert.str();

      while (numString.length() < nzeros)
	{
	  numString = "0" + numString;
	}

      fileString = prefixString + numString;

      strcpy(outputFile, fileString.c_str());

      output = fopen(outputFile, "w");
      fprintf(output, "%+.4E %i %i \n", time, nLatitude, nLongitude);

      for (int j = 0; j < nLongitude; j++)
	{
	  for (int k = 0; k < nLatitude; k++)
	    {
	      fprintf(output, "%+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E \n",
		      longitude[j], latitude[k], flux[j][k], planetTeff[j][k], Ngamma[j][k], darkness[j][k],
		      altitude[j][k], azimuth[j][k], hourAngle[j]);
	    }
	}
      fflush(output);
      fclose(output);

    }

  fclose(outputlog);
  // End of program
  return 0;
}


int readParameters(ifstream &inputStream, string &prefixString, int &nTime, int &nLatitude, int &nLongitude, double &PsToPo, 
		   double &mstar, double &lstar, double &Teff, double &mplanet, double &albedo, double &obliquity, double &semimaj, double &ecc,
		   double &inc, double &longascend, double &argper)
{

  string par, line;
  double mplanetearth;
  double mEarthToMSol = 3.0034e-6;

  // check that the file exists as cpp does not check
  //and it just returns 0 if it doesnt exist
  if (inputStream == 0)
    {
      cout << "No Input file found" << endl;
      return -1;
    }

  // Then loop through each line using getline and then
  //assign to vectors
  while (getline(inputStream, line))
    {
      istringstream iss(line);
      iss >> par;

      if (par == "OutputPrefix")
	{

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
	  mplanet = mplanetearth * mEarthToMSol; // convert to Solar masses

	}

      if (par == "PlanetAlbedo")
	{
	  iss >> albedo;

	}

      if (par == "Obliquity")
	{
	  iss >> obliquity;
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
  inputStream.close();
  return 0;
}


// Function takes in Star and Planet Objects as input, calculates longitude/latitude flux map from star

void calcFlux(Star &star, Planet &planet,double &longitude, double &latitude, double &hourAngle, 
	      double &flux, double &altitude, double &azimuth,
	      double &time, double &Pspin, double &obliquity)
{
  double rdotn,declination, long_apparent;
  double magpos,lstar;
  Vector3D pos, unitpos;
  Vector3D zvector(0.0,0.0,1.0);

  // Get position of planet relative to star
  pos = (star.getPosition()).relativeVector(planet.getPosition());
  magpos = pos.magVector();
  unitpos = pos.unitVector();
  lstar = star.getLuminosity();

  // Longitude corresponding to noon - assumes noon at t=0 is 0
  //noon = fmod(2.0 * pi * time / Pspin, 2.0 * pi);
  
  // Declination of the Sun - angle between planet's position vector and equator (at noon)
  
  Vector3D decVector(unitpos.elements[0], unitpos.elements[1],
		     unitpos.elements[2]);

  // Rotate this vector if planet has non-zero obliquity
  if (obliquity != 0.0)
    {
      decVector.rotateX(obliquity);
    }

  rdotn = unitpos.dotProduct(decVector);
  declination = safeAcos(rdotn);

  // Now begin calculation over latitude and longitude points

  // Rotate planet according to its spin period

  long_apparent = fmod(longitude + 2.0 * pi * time / Pspin, 2.0 * pi);

  // Calculate hour angle - distance between current longitude and noon

  // Distance between longitude and noon = hourAngle
  // hour angle between -180 and +180

  Vector3D longSurface(cos(long_apparent), sin(long_apparent), 0.0);

  rdotn = unitpos.dotProduct(longSurface);
  hourAngle = safeAcos(rdotn);

  if((unitpos.crossProduct(longSurface)).dotProduct(zvector) > 0.0)
    {
      hourAngle = -hourAngle;
    }

  //cout << "        longitude  " << longitude << "  hour angle " << hourAngle << " Declination "<< declination << endl;

  // construct surface normal vector at this longitude and latitude

  Vector3D surface(sin(latitude) * cos(long_apparent), sin(latitude) * sin(long_apparent), cos(latitude));

  surface = surface.unitVector();

  // If necessary, rotate surface vector by obliquity
  
  if (obliquity != 0.0)
    {
      surface.rotateX(obliquity);
    }
  
  // take the dot product with the unit position vector
  
  rdotn = unitpos.dotProduct(surface);
 
  // Calculate fluxes
  // if position.surface is less than zero, long/lat location is not illuminated
  
  if (rdotn > 0.0)
    {
      
      flux = flux + lstar * rdotn / (4.0 * pi * magpos * magpos);
    } 

  // Calculate altitude and azimuthal position on sky, and angular size
  // Formulae not exactly as used normally
  // As we measure latitudes from 0 to 180, not -90 to 90, some sines have become cosines, and vice versa
  // (some extra plus and minus signs as a result)
  
  //altitude = cos(declination) * cos(latitude) * cos(hourAngle) + sin(declination) * sin(latitude);  // These calculations assume lat->-90,90 deg

  altitude = - cos(declination) * cos(hourAngle) * sin(latitude) + sin(declination) * cos(latitude);  // These calculations assume lat->0,180 deg
  altitude = asin(altitude);

  // Azimuth angle measured from north, clockwise

  if(cos(altitude)*sin(latitude)!=0.0){

	//azimuth = (sin(declination) - sin(altitude)*cos(latitude))/(cos(altitude)*cos(latitude)); // These calculations assume lat->-90, 90 deg
    azimuth = (sin(altitude)*sin(latitude) - sin(declination))/(cos(altitude)*sin(latitude)); // These calculations assume lat->0, 180 deg
    azimuth = safeAcos(azimuth);
  }
  else
    {
      azimuth = 0.0;
    }

  // If hour angle positive (afternoon) azimuth is > 180
  
  if(hourAngle>0.0)
    {
      azimuth = 2.0*pi - azimuth;
    }
  
}




// Simple function to stop acos becoming infinite if abs(x) > 1.0
double safeAcos(double x)
{
  if (x < -1.0)
    x = -1.0;
  else if (x > 1.0)
    x = 1.0;
  return acos(x);
}
