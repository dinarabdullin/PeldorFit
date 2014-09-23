#ifndef MAGNETICFIELD_H
#define MAGNETICFIELD_H

#include <vector>
using std::vector;
#include "Definitions.h"
#include "RandGens.h"
#include "Rotations.h"


class MagneticField
{
private:
	vector<double> projections;

public:
	MagneticField(double& amplitude, long& norient);

	// Generate random orientations of the magnetic field unit vector in a half sphere
	void randomOrient(RandGens& randgen);

	// Get the z-component of magnetic field vector
	double zprojection(long const& num) const;

	// Get the magnetic field vector
	double* projection(long const& num) const;
	
	// Calculate the projection of the magnetic field vector in the spin B frame
	double* frameB(EulerRotation const& RotationMatrix, double const initialVector[]) const;

    double value;
	long N;
};
 
#endif