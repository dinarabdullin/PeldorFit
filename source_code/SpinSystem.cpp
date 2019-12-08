#include "SpinSystem.h"
#include <random>
#include <cmath>


SpinSystem::SpinSystem(double g_in[3], vector<int>& n_in, vector<double>& I_in, vector<double>& A_in,
                       double gStrain_in[3], double AStrain_in[3], double& lwpp_in)
	: n(n_in), I(I_in), A(A_in), lwpp(lwpp_in)
{
	F = 13.99624; // F = betta/h in GHz/T

	for (size_t i = 0; i < 3; i++) {
		g[i] = g_in[i];
		gStrain[i] = gStrain_in[i];
		AStrain[i] = AStrain_in[i];
	}
	if ((gStrain[0] != 0) || (gStrain[1] != 0) || (gStrain[2] != 0)) gStrain_on = true;
	else gStrain_on = false;
	if ((AStrain[0] != 0) || (AStrain[1] != 0) || (AStrain[2] != 0)) AStrain_on = true;
	else AStrain_on = false;

	// Consider the cases of 0, 1, 2 different nuclei separately
	if ( n.empty() ) { 
		nstates = 1;
		ncomp = 1;
		intensity.reserve(ncomp);
		intensity.push_back(1); 
	}
	else if ( n.size() == 1 ) {
		nstates = (size_t) ( pow( (2*I[0] + 1), n[0]) );
		ncomp = (size_t) ( 2*I[0]*n[0] + 1 );
		intensity.reserve(ncomp);
		int nucleusIndex = static_cast<int>( 2*I[0] - 1 );
		int numberIndex = n[0] - 1;
		for (size_t i = 0; i < ncomp; i++) intensity.push_back( relativeIntensity[nucleusIndex][numberIndex][i] );
	}
	else if ( n.size() == 2 ) {
		nstates = (size_t) ( pow( (2*I[0] + 1), n[0]) * pow( (2*I[1] + 1), n[1]) );
		size_t ncomp0 = (size_t) ( 2*I[0]*n[0] + 1 );
		size_t ncomp1 = (size_t) ( 2*I[1]*n[1] + 1 );
		ncomp = ncomp0 * ncomp1;
		intensity.reserve(ncomp);
		int nucleusIndex0 = static_cast<int>( 2*I[0] - 1 );
		int numberIndex0 = n[0] - 1;
		int nucleusIndex1 = static_cast<int>( 2*I[1] - 1 );
		int numberIndex1 = n[1] - 1;
		for (size_t i = 0; i < ncomp0; i++) {
			for (size_t j = 0; j < ncomp1; j++) {
				intensity.push_back( relativeIntensity[nucleusIndex0][numberIndex0][i] * relativeIntensity[nucleusIndex1][numberIndex1][j] );
			}
		}
	}
}

double SpinSystem::effectiveGfactor(double const fprojection[]) const
{
	double gEff = sqrt( pow((g[0]*fprojection[0]), 2) + pow((g[1]*fprojection[1]), 2) + pow((g[2]*fprojection[2]), 2) );
	return gEff;
}

double* SpinSystem::resonanceFrequencies(double const fprojection[], double const& fvalue, RandGens* const prandgen) const
{
	double* resFreq = new double[ncomp];

	// Initialize a normal distribution function
	std::normal_distribution<double> unormaldistr(0,1);
	
	// g-factor
	double gEff(0), dgEff(0);
	gEff = sqrt( pow((g[0]*fprojection[0]), 2) + pow((g[1]*fprojection[1]), 2) + pow((g[2]*fprojection[2]), 2) );
	if (gStrain_on) {
		dgEff = (g[0]*gStrain[0]*fprojection[0]*fprojection[0] + g[1]*gStrain[1]*fprojection[1]*fprojection[1] + g[2]*gStrain[2]*fprojection[2]*fprojection[2]) / gEff;
		gEff += 0.424661*dgEff * unormaldistr(prandgen->prng); // FWHM = 2.35482 * sigma
	}

	// Consider the cases of 0, 1, 2 different nuclei separately
	// The case of 0 nuclei
	if ( n.empty() ) {
		// Inhomogenious broadering 
		double df_inhom = 1e-3 * 0.5*lwpp * unormaldistr(prandgen->prng); 
		resFreq[0] =  F * gEff * fvalue + df_inhom;
	}
	// The case of 1 sort of nuclei
	if ( n.size() == 1 ) {
		// Hyperfine coupling constant (GHz)
		double AEff(0), dAEff(0);
		AEff = 1e-3 * sqrt( pow((A[0]*fprojection[0]), 2) + pow((A[1]*fprojection[1]), 2) + pow((A[2]*fprojection[2]), 2) );
		if (AStrain_on) { 
			dAEff = 1e-6 * (A[0]*AStrain[0]*fprojection[0]*fprojection[0] + A[1]*AStrain[1]*fprojection[1]*fprojection[1] + A[2]*AStrain[2]*fprojection[2]*fprojection[2]) / AEff;
			AEff += 0.424661*dAEff * unormaldistr(prandgen->prng);
		}
		// Inhomogenious broadering 
		double df_inhom(0);
		// Nuclear quantum number
		double m = -I[0]*n[0]; 
		for (size_t i = 0; i < ncomp; ++i) {
			// Inhomogenious broadering
			df_inhom = 1e-3 * 0.5*lwpp * unormaldistr(prandgen->prng);
			resFreq[i] =  F * gEff * fvalue + AEff * m + df_inhom;
			++m;
		}
	}
	// The case of 2 sorts of nuclei
	if ( n.size() == 2 ) {
		// Hyperfine coupling constants (GHz)
		double AEff0(0), dAEff0(0), AEff1(0);
		AEff0 = 1e-3 * sqrt( pow((A[0]*fprojection[0]), 2) + pow((A[1]*fprojection[1]), 2) + pow((A[2]*fprojection[2]), 2) );
		AEff1 = 1e-3 * sqrt( pow((A[3]*fprojection[0]), 2) + pow((A[4]*fprojection[1]), 2) + pow((A[5]*fprojection[2]), 2) );
		if (AStrain_on) { 
			dAEff0 = 1e-6 * (A[0]*AStrain[0]*fprojection[0]*fprojection[0] + A[1]*AStrain[1]*fprojection[1]*fprojection[1] + A[2]*AStrain[2]*fprojection[2]*fprojection[2]) / AEff0;
			AEff0 += 0.424661*dAEff0 * unormaldistr(prandgen->prng);
		}
		// Inhomogenious broadering 
		double df_inhom(0);
		// Number of spectral components for each nucleus
		size_t ncomp0 = (size_t) ( 2*I[0]*n[0] + 1 );
		size_t ncomp1 = (size_t) ( 2*I[1]*n[1] + 1 );
		// Nuclear quantum numbers
		double m0 = -I[0] * n[0];
		double m1 = -I[1] * n[1];
		for (size_t i = 0; i < ncomp0; ++i) {
			m1 = -I[1] * n[1];
			for (size_t j = 0; j < ncomp1; ++j) {
				df_inhom = 1e-3 * 0.5*lwpp * unormaldistr(prandgen->prng);
				resFreq[ncomp1*i+j] =  F * gEff * fvalue + AEff0 * m0  + AEff1 * m1 + df_inhom;
				++m1;
			}
			++m0;
		}
	}

	return resFreq;
}
