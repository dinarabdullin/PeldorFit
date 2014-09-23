#ifndef ENVIROMENT_H
#define ENVIROMENT_H

#include <vector>
using std::vector;
#include <string>
using std::string;
#include "Population.h"
#include "SpinSystem.h"
#include "MagneticField.h"
#include "Pulses.h"
#include "Rotations.h"
#include "RandGens.h"


// Enviroment class
class Enviroment 
{        
private: 
    size_t nof;
	vector<bool> pswitches;
	double limit;
	SpinSystem const* pspinA;
	SpinSystem const* pspinB;
	vector<MagneticField> magneticFields;
	vector<Pulses> pulseSequences;
    vector<vector<double>> expTime;
    vector<vector<double>> expSignal;
    vector<double> expMean;
	vector<double> expVariance;

public:
	Enviroment(size_t const& nOffset, vector<bool> const& paramswitches);

	void createSpins(double g1[3],vector<int>& n1,vector<double>& I1,vector<double>& A1,double gStrain1[3],double AStrain1[3],double& lwpp1,
	                 double g2[3],vector<int>& n2,vector<double>& I2,vector<double>& A2,double gStrain2[3],double AStrain2[3],double& lwpp2);

	void createField(vector<double>& valueField, long& nAvg);

	void createPulses(vector<double>& detPiLength, vector<double>& detPiHalfLength, vector<double>& pumpPiLength,
	                  vector<double>& detFreq, vector<double>& pumpFreq);
    
    void setData(vector<vector<double>> const& times, vector<vector<double>> const& signals, int& scoreIndex);
    
    vector<double> computeTimetrace(const vector<double> &parameters, 
									vector<double> const& timepoints,
									SpinSystem const& spinA, 
									SpinSystem const& spinB, 
									MagneticField const& magneticField, 
									Pulses const& pulses,
									long const& nAvg,
									size_t const& offset) const;
    
	// Calculate a RMSD
    double rmsd(vector<double> const& X, vector<double> const& expY, vector<double> const& calcY) const;
    
	// Calculate a RMSD / Pearson's coefficient
    double rmsd_PearsonCC(vector<double> const& X, vector<double> const& expY, vector<double> const& calcY,
	                      double const& expYMean, double const& expYVar) const;
    
	// Calculate a Pearson's coefficient
    double PearsonCC(vector<double> const& X, vector<double> const& expY, vector<double> const& calcY,
					 double const& expYMean, double const& expYVar) const;
    
	// Score a population
	void runScoring(Population& population, int const& scoreIndex, long const& nAvg) const;

	// Record the score
	void recordScore(double const& score, string const& file_root, size_t const& numOfPop) const;
	
	// Record optimized parameters
	void recordParameters(vector<double> const& parameters, string const& file_root) const;
	
	// Record the fit
	void recordFit(vector<double> const& parameters, string const& file_root, size_t const& numOfPop) const;

	// Record the EPR spectrum of the spin pair
    void recordSpectrum(vector<double> const& parameters, string const& file_root) const;	

	// Record the symmetry-related solutions
	void recordSymmetricSolutions(vector<double> const& parameters, vector<double> const& ranges, int const& scoreIndex, string const& file_root) const;
	
	// Calculate the modulation depth parameter for individual conformers
    vector<double> computeModDepth(vector<double> const& parameters, SpinSystem const& spinA, SpinSystem const& spinB,
                                   MagneticField const& magneticField, Pulses const& pulses) const;
    
    // Record the modulation depth parameter vs dipolar angle
    void recordModDepth(vector<double> const& parameters, string const& file_root) const;

    // Record error profile
    void recordErrorProfile(vector<double> const& parameters, vector<int> const& variables, vector<double> const& ranges, string const& file_root) const;
};
 
#endif