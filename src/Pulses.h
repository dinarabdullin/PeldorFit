#ifndef PULSES_H
#define PULSES_H

#include "Definitions.h"

class Pulses
{
private:
	double detFreq;
	double pumpFreq;
	double detFreqB1;
	double detFreqB1s;
	double pumpFreqB1;

public:
	Pulses(double& m_detPiLength, double& m_detPiHalfLength, double& m_pumpPiLength, double& m_detFreq, double& m_pumpFreq);

	double detExcitation(double const& freq) const;

	double pumpExcitation(double const& freq) const;

	double decay(double const& dipFreq) const;
    
    double detPiLength;
	double detPiHalfLength;
	double pumpPiLength;
};

#endif

