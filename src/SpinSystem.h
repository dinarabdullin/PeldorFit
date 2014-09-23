#ifndef SPINSYSTEM_H
#define SPINSYSTEM_H

#include <vector>
using std::vector;
#include "RandGens.h"
#include "Definitions.h"


class SpinSystem
{
private:
	double F;
	double g[3];
    vector<int> n;
    vector<double> I;
	vector<double> A;
	bool gStrain_on;
	bool AStrain_on;
	double gStrain[3];
	double AStrain[3];
	double lwpp;

public:
	size_t nstates;        // Number of microstates
    size_t ncomp;          // Number of spectral components
    vector<int> intensity; // Intensities of the individual spectral components

	SpinSystem(double g_in[3], vector<int>& n_in, vector<double>& I_in, vector<double>& A_in,
	           double gStrain_in[3], double AStrain_in[3], double& lwpp_in);

	// Compute effective g-factor
	double effectiveGfactor(double const fprojection[]) const;

	// Compute resonance frequencies
	double* resonanceFrequencies(double const fprojection[], double const& fvalue, RandGens* const prandgen) const;
};

#endif
