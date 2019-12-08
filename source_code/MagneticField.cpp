#include "MagneticField.h"
#include <cmath>


MagneticField::MagneticField(double& amplitude, long& norient) : value(amplitude), N(norient)
{
	projections.reserve( 3*norient );
	RandGens* prandgen = new RandGens;
	for (long f = 0; f < norient; ++f) { randomOrient(*prandgen); }
	delete prandgen;
}

void MagneticField::randomOrient(RandGens& randgen)
{
	double fphi = randgen.uniformRand(0.0, 2.0*PI);
	//double ftheta = randgen.sineweighted(); // sine weighted distribution on the (0, pi/2) interval
	double ftheta;
	if (randgen.uniformRand() <= 0.5) ftheta = randgen.sineweighted();
	else                              ftheta = PI - randgen.sineweighted();
	

	projections.push_back( sin(ftheta)*cos(fphi) );
	projections.push_back( sin(ftheta)*sin(fphi) );
	projections.push_back( cos(ftheta) );
}

double MagneticField::zprojection(long const& num) const
{
	double zprojection = projections[3*num+2];
	return zprojection;
}

double* MagneticField::projection(long const& num) const
{
	double* projection = new double[3];
	projection[0] = projections[3*num];
	projection[1] = projections[3*num+1];
	projection[2] = projections[3*num+2];
	return projection;
}

double* MagneticField::frameB(EulerRotation const& RotationMatrix, double const initialVector[]) const
{
	double* newVector = new double[3];
	for(size_t i = 0; i < 3; i++) { 
		//newVector[i] = RotationMatrix.R[i][0] * initialVector[0] + RotationMatrix.R[i][1] * initialVector[1] + RotationMatrix.R[i][2] * initialVector[2]; // changed recently
		newVector[i] = RotationMatrix.RT[i][0] * initialVector[0] + RotationMatrix.RT[i][1] * initialVector[1] + RotationMatrix.RT[i][2] * initialVector[2];
	}
	return newVector;
}