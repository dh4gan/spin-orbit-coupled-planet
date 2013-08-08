/*
 * Star.h
 *
 *  Created on: Nov 8, 2012
 *      Author: dh4gan
 */

#include "Body.h"

using namespace std;

#ifndef STAR_H_
#define STAR_H_

class Star: public Body {
public:

	Star();
	Star(string &namestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel);
	Star(string &namestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel, double &T, int &n);
	virtual ~Star();

	void setLuminosity(double lum){luminosity = lum;}
	double getLuminosity() {return luminosity;}

	void setTeff(double T){Teff=T;}
	double getTeff() {return Teff;}

	void setLambdaMin(double l){lambda_min=l;}
	double getLambdaMin() {return lambda_min;}

	void setLambdaMax(double l){lambda_max=l;}
	double getLambdaMax() {return lambda_max;}

	void setNLambda(int n){nlambda = n;}
	int getNLambda(){return nlambda;}

	void setInnerHZ(double r){innerHZ = r;}
	double getInnerHZ(){return innerHZ;}

	void setOuterHZ(double r){outerHZ = r;}
	double getOuterHZ(){return innerHZ;}

	vector<double> getILambda(){return I_lambda;}

	// Standard cloning method
	virtual Star* Clone() { return new Star(*this); }

	// Calculation Methods

	void calculateBlackbodySpectrum();
	void calculateSingleHZ();

protected:

	double luminosity; // Luminosity of Star IN SOLAR UNITS
	double Teff;
	double lambda_min;
	double lambda_max;
	double innerHZ;
	double outerHZ;

	int nlambda;
	vector<double> I_lambda;

};

#endif /* STAR_H */

