#include "Pulses.h"
#include <cmath>


Pulses::Pulses(double& m_detPiLength, double& m_detPiHalfLength, double& m_pumpPiLength, double& m_detFreq, double& m_pumpFreq) 
	: detPiLength(m_detPiLength), detPiHalfLength(m_detPiHalfLength), pumpPiLength(m_pumpPiLength), detFreq(m_detFreq), pumpFreq(m_pumpFreq)
{
	detFreqB1 = 0.25 / m_detPiHalfLength;
	detFreqB1s = 0.5 / m_detPiLength;
	pumpFreqB1 = 0.5 / m_pumpPiLength;
}

double Pulses::detExcitation(double const& freq) const
{
	double freqEff(0), freqEffs(0), excitation(0);
	if (detFreqB1 == detFreqB1s) {
		freqEff = sqrt( pow(detFreq-freq, 2) + pow(detFreqB1, 2) );
		excitation = pow( (detFreqB1/freqEff) * sin(2*PI * freqEff * detPiHalfLength), 5); 
	}
	else {
		freqEff = sqrt( pow(detFreq-freq, 2) + pow(detFreqB1, 2) );
		freqEffs = sqrt( pow(detFreq-freq, 2) + pow(detFreqB1s, 2) );
		excitation = (detFreqB1/freqEff) * sin(2*PI * freqEff * detPiHalfLength) * pow( (detFreqB1s/freqEffs) * sin(0.5 * 2*PI * freqEffs * detPiLength), 4);
	}
	return fabs(excitation);
}

double Pulses::pumpExcitation(double const& freq) const
{
	double freqEff = sqrt( pow(pumpFreq-freq, 2) + pow(pumpFreqB1, 2) );
	double excitation = pow( (pumpFreqB1/freqEff) * sin(0.5 * 2*PI * freqEff * pumpPiLength), 2);
	return fabs(excitation);
}

double Pulses::decay(double const& dipFreq) const
{
	return exp( -pow(dipFreq/(2*PI*17), 2) ); // For the 32 ns detection pulse and 12 ns pump pulse
}