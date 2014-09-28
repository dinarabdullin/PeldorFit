#include "Enviroment.h"
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::setw;
using std::setprecision;
#include <fstream>
using std::fstream;
using std::ios;
using std::ofstream;
#include <sstream>
using std::ostringstream;
#include <numeric>
using std::accumulate;
#include <cmath>
#ifdef __gnu_linux__
	#include <sys/time.h>
#endif
#include "tbb/tbb.h"
using tbb::concurrent_vector;
using tbb::parallel_for;
using tbb::mutex;
static mutex mtx;


Enviroment::Enviroment(size_t const& nOffset, vector<bool> const& paramswitches) 
	: nof(nOffset), pswitches(paramswitches)
{
	limit = 1.e-6;
}

void Enviroment::createSpins(double g1[3],vector<int>& n1,vector<double>& I1,vector<double>& A1,double gStrain1[3],double AStrain1[3],double& lwpp1,
	                         double g2[3],vector<int>& n2,vector<double>& I2,vector<double>& A2,double gStrain2[3],double AStrain2[3],double& lwpp2)
{
	pspinA = new SpinSystem(g1,n1,I1,A1,gStrain1,AStrain1,lwpp1);
	pspinB = new SpinSystem(g2,n2,I2,A2,gStrain2,AStrain2,lwpp2);
}

void Enviroment::createField(vector<double>& valueField, long& nAvg)
{
	magneticFields.reserve(nof);
	for (size_t i = 0; i < nof ; ++i) {
		MagneticField const* tempField = new MagneticField(valueField[i],nAvg);
		magneticFields.push_back(*tempField);
		delete tempField;
	}
}

void Enviroment::createPulses(vector<double>& detPiLength,vector<double>& detPiHalfLength,vector<double>& pumpPiLength,
							  vector<double>& detFreq,vector<double>& pumpFreq)
{
	pulseSequences.reserve(nof);
	for (size_t i = 0; i < nof; ++i) {
		Pulses const* tempPulses = new Pulses(detPiLength[i],detPiHalfLength[i],pumpPiLength[i],detFreq[i],pumpFreq[i]);
		pulseSequences.push_back(*tempPulses);
		delete tempPulses;
	}
}

void Enviroment::setData(vector<vector<double>> const& times,vector<vector<double>> const& signals,int& scoreIndex)
{
	// Set the experimental data
	expTime.reserve(nof);
	expSignal.reserve(nof);
	for (size_t i = 0; i < times.size(); ++i) {
		expTime.push_back(times[i]);
		expSignal.push_back(signals[i]);
	}
	
	// Some calculations for the Pearson's product-moment coefficient
	double t_expMean, t_expVariance;
	expMean.reserve(nof);
	expVariance.reserve(nof);
	for (size_t i = 0; i < nof; ++i) {
		t_expMean = accumulate(expSignal[i].begin(), expSignal[i].end(), 0.0) / (double) expSignal[i].size();
		t_expVariance = 0.0;
		for (size_t j = 0; j < expSignal[i].size(); ++j) t_expVariance += pow((expSignal[i][j] - t_expMean), 2);
		expMean.push_back(t_expMean);
		expVariance.push_back(t_expVariance);
	}
}

vector<double> Enviroment::computeTimetrace(const vector<double> &parameters, 
											vector<double> const& timepoints,
											SpinSystem const& spinA, 
											SpinSystem const& spinB, 
											MagneticField const& magneticField, 
											Pulses const& pulses,
											long const& nAvg,
											size_t const& offset) const
{		
	// Initialize a calculated signal with 0 values
	size_t const nTime = timepoints.size(); // Set the number of time points
	double* calcSignal = new double[nTime];
	for (size_t t = 0; t < nTime; ++t) { calcSignal[t] = 0.0; }
	
	// Set the values of geometric parameters
	double dist1_mean,  dist1_width,  dist2_mean,  dist2_width;
	double theta1_mean, theta1_width, theta2_mean, theta2_width;
	double phi1_mean,   phi1_width,   phi2_mean,   phi2_width;
	double alpha1_mean, alpha1_width, alpha2_mean, alpha2_width;
	double betta1_mean, betta1_width, betta2_mean, betta2_width;
	double gamma1_mean, gamma1_width, gamma2_mean, gamma2_width;
	double ratio;
	double jmean, jwidth;
	double damping;
	int np(-1);
	if ( pswitches[0] ) {
		dist1_mean = parameters[++np];
		if ( pswitches[1] ) dist1_width = parameters[++np];
		else          dist1_width = 0;
	}
	else { dist1_mean = 0; dist1_width = 0; }
	if ( pswitches[2] ) {
		theta1_mean = parameters[++np];
		if ( pswitches[3] ) theta1_width = parameters[++np];
		else          theta1_width = 0;
	}
	else { theta1_mean = 0; theta1_width = 0; }
	if ( pswitches[4] ) {
		phi1_mean = parameters[++np];
		if ( pswitches[5] ) phi1_width = parameters[++np];
		else          phi1_width = 0;
	}
	else { phi1_mean = 0; phi1_width = 0; }
	if ( pswitches[6] ) {
		alpha1_mean = parameters[++np];
		if ( pswitches[7] ) alpha1_width = parameters[++np];
		else          alpha1_width = 0;
	}
	else { alpha1_mean = 0; alpha1_width = 0; }
	if ( pswitches[8] ) {
		betta1_mean = parameters[++np];
		if ( pswitches[9] ) betta1_width = parameters[++np];
		else          betta1_width = 0;
	}
	else { betta1_mean = 0; betta1_width = 0; }
	if ( pswitches[10] ) {
		gamma1_mean = parameters[++np];
		if ( pswitches[11] ) gamma1_width = parameters[++np];
		else           gamma1_width = 0;
	}
	else { gamma1_mean = 0; gamma1_width = 0; }
	if ( pswitches[12] ) {
		dist2_mean = parameters[++np];
		if ( pswitches[13] ) dist2_width = parameters[++np];
		else           dist2_width = dist1_width;
	}
	else { dist2_mean = dist1_mean; dist2_width = dist1_width; }
	if ( pswitches[14] ) {
		theta2_mean = parameters[++np];
		if ( pswitches[15] ) theta2_width = parameters[++np];
		else           theta2_width = theta1_width;
	}
	else { theta2_mean = theta1_mean; theta2_width = theta1_width; }
	if ( pswitches[16] ) {
		phi2_mean = parameters[++np];
		if ( pswitches[17] ) phi2_width = parameters[++np];
		else           phi2_width = phi1_width;
	}
	else { phi2_mean = phi1_mean; phi2_width = phi1_width; }
	if ( pswitches[18] ) {
		alpha2_mean = parameters[++np];
		if ( pswitches[19] ) alpha2_width = parameters[++np];
		else           alpha2_width = alpha1_width;
	}
	else { alpha2_mean = alpha1_mean; alpha2_width = alpha1_width; }
	if ( pswitches[20] ) {
		betta2_mean = parameters[++np];
		if ( pswitches[21] ) betta2_width = parameters[++np];
		else           betta2_width = betta1_width;
	}
	else { betta2_mean = betta1_mean; betta2_width = betta1_width; }
	if ( pswitches[22] ) {
		gamma2_mean = parameters[++np];
		if ( pswitches[23] ) gamma2_width = parameters[++np];
		else           gamma2_width = gamma1_width;
	}
	else { gamma2_mean = gamma1_mean; gamma2_width = gamma1_width;}
	if ( pswitches[24] ) ratio = parameters[++np];
	else           ratio = 1;
	// J coupling constant
	if ( pswitches[25] ) jmean = parameters[++np];
	else           jmean = 0;
	if ( pswitches[26] ) jwidth = parameters[++np];
	else           jwidth = 1e-12;
	// Excitation damping
	if ( pswitches[27] == false ) damping = 1.0;
	else {
		if (pswitches.size() == 28) damping = parameters[++np];
		else                  damping = parameters[np+offset+1];
	}

	////////////////////////////////////////////////////////////////////////////
	//                      COMPUTE THE TIME TRACE                            //
	////////////////////////////////////////////////////////////////////////////

	RandGens* prandgen = new RandGens;            // Initialize random generator
	std::normal_distribution<double> unormal_distribution(0,1); // Initialize a normal distribution
	size_t const compA = spinA.ncomp;             // Set the number of spectral components of the spin A
	double const statesA = (double) spinA.nstates;// Set the number of microstates of the spin A
	size_t const compB = spinB.ncomp;             // Set the number of spectral components of the spin B
	double const statesB = (double) spinB.nstates;// Set the number of microstates of the spin B
	bool first_ensemble;
	double distance(0),theta(0), phi(0), alpha(0), betta(0), gamma(0);
	double gA(0), gB(0);
	double* projectionA(nullptr);
	double* projectionB(nullptr);
	double* freqA(nullptr);
	double* freqB(nullptr);
	double detExcitA(0), detExcitB(0), pumpExcitA(0), pumpExcitB(0);
	double dipProj(0);
	double D(0);
	double JC(0);
	double dipFreq(0);
	double ampl(0);
		
	// Itterate over the conformers
	for (long n = 0; n < nAvg; n++) {

		// Choose between 2 ensembles of parameters
		if (prandgen->uniformRand() <= ratio) first_ensemble = true;
		else                                  first_ensemble = false;

		// Generate the set of geometrical parameters of the system
		if (first_ensemble) {       // if the 1st ensemble of parameters was chosen
			if (dist1_width == 0)  distance = dist1_mean;
			else                   distance = dist1_mean + dist1_width * unormal_distribution(prandgen->prng);
			if (theta1_width == 0) theta = theta1_mean;
			else                   theta = prandgen->uniformRand(theta1_mean-0.5*theta1_width, theta1_mean+0.5*theta1_width);
			if (phi1_width == 0)   phi = phi1_mean;
			else                   phi = prandgen->uniformRand(phi1_mean-0.5*phi1_width, phi1_mean+0.5*phi1_width);
			if (alpha1_width == 0) alpha = alpha1_mean;
			else                   alpha = prandgen->uniformRand(alpha1_mean-0.5*alpha1_width, alpha1_mean+0.5*alpha1_width);
			if (betta1_width == 0) betta = betta1_mean;
			else                   betta = prandgen->uniformRand(betta1_mean-0.5*betta1_width, betta1_mean+0.5*betta1_width);
			if (gamma1_width == 0) gamma = gamma1_mean;
			else                   gamma = prandgen->uniformRand(gamma1_mean-0.5*gamma1_width, gamma1_mean+0.5*gamma1_width);
		}
		else {                      // if the 2nd ensemble of parameters was chosen
			if (dist2_width == 0)  distance = dist2_mean;
			else                   distance = dist2_mean + dist2_width * unormal_distribution(prandgen->prng);
			if (theta2_width == 0) theta = theta2_mean;
			else                   theta = prandgen->uniformRand(theta2_mean- 0.5*theta2_width, theta2_mean+0.5*theta2_width);
			if (phi2_width == 0)   phi = phi2_mean;
			else                   phi = prandgen->uniformRand(phi2_mean- 0.5*phi2_width, phi2_mean+0.5*phi2_width);
			if (alpha2_width == 0) alpha = alpha2_mean;
			else                   alpha = prandgen->uniformRand(alpha2_mean-0.5*alpha2_width, alpha2_mean+0.5*alpha2_width);
			if (betta2_width == 0) betta = betta2_mean;
			else                   betta = prandgen->uniformRand(betta2_mean-0.5*betta2_width, betta2_mean+0.5*betta2_width);
			if (gamma2_width == 0) gamma = gamma2_mean;
			else                   gamma = prandgen->uniformRand(gamma2_mean-0.5*gamma2_width, gamma2_mean+0.5*gamma2_width);
		}
		/*
		if (first_ensemble) {       // if the 1st ensemble of parameters was chosen
			if (dist1_width == 0)  distance = dist1_mean;
			else                   distance = dist1_mean + dist1_width * unormal_distribution(prandgen->prng);
			if (theta1_width == 0) theta = theta1_mean;
			else                   theta = theta1_mean + theta1_width * unormal_distribution(prandgen->prng);
			if (phi1_width == 0)   phi = phi1_mean;
			else                   phi = phi1_mean + phi1_width * unormal_distribution(prandgen->prng);
			if (alpha1_width == 0) alpha = alpha1_mean;
			else                   alpha = alpha1_mean + alpha1_width * unormal_distribution(prandgen->prng);
			if (betta1_width == 0) betta = betta1_mean;
			else                   betta = betta1_mean + betta1_width * unormal_distribution(prandgen->prng);
			if (gamma1_width == 0) gamma = gamma1_mean;
			else                   gamma = gamma1_mean + gamma1_width * unormal_distribution(prandgen->prng);
		}
		else {                      // if the 2nd ensemble of parameters was chosen
			if (dist2_width == 0)  distance = dist2_mean;
			else                   distance = dist2_mean + dist2_width * unormal_distribution(prandgen->prng);
			if (theta2_width == 0) theta = theta2_mean;
			else                   theta = theta2_mean + theta2_width * unormal_distribution(prandgen->prng);
			if (phi2_width == 0)   phi = phi2_mean;
			else                   phi = phi2_mean + phi2_width * unormal_distribution(prandgen->prng);
			if (alpha2_width == 0) alpha = alpha2_mean;
			else                   alpha = alpha2_mean + alpha2_width * unormal_distribution(prandgen->prng);
			if (betta2_width == 0) betta = betta2_mean;
			else                   betta = betta2_mean + betta2_width * unormal_distribution(prandgen->prng);
			if (gamma2_width == 0) gamma = gamma2_mean;
			else                   gamma = gamma2_mean + gamma2_width * unormal_distribution(prandgen->prng);
		}
		*/

		// Rotation matrix between the spin A and spin B coordinate systems
		EulerRotation const* pFrameRotation = new EulerRotation(alpha,betta,gamma);
		
		// Compute the orientation of the magnetic field in the frame of the spin A
		projectionA = magneticField.projection(n);
		// Compute the orientation of the magnetic field in the frame of the spin B
		projectionB = magneticField.frameB(*pFrameRotation, projectionA);
		// Compute the cosine of the dipolar angle
		dipProj = pow( (projectionA[0] * cos(phi)*sin(theta) + projectionA[1] * sin(phi)*sin(theta) + projectionA[2] * cos(theta)), 2);

		// Compute an effective g-factor of the spin A
		gA = spinA.effectiveGfactor(projectionA);
		// Compute an effective g-factor of the spin B
		gB = spinB.effectiveGfactor(projectionB);

		// Compute the resonance frequences of the spin A
		freqA = spinA.resonanceFrequencies(projectionA, magneticField.value, prandgen);
		// Compute the resonance frequence of the spin B
		freqB = spinB.resonanceFrequencies(projectionB, magneticField.value, prandgen);

		// Compute the excitation efficiency of the spin A by detection pulses
		detExcitA = 0;
		for (size_t k = 0; k < compA; k++) detExcitA += static_cast<double>(spinA.intensity[k]) * pulses.detExcitation( freqA[k] );
		detExcitA /= statesA;
		// Compute the excitation efficiency of the spin B by detection pulses
		detExcitB = 0;
		for (size_t k = 0; k < compB; k++) detExcitB += static_cast<double>(spinB.intensity[k]) * pulses.detExcitation( freqB[k] );
		detExcitB /= statesB;
			
		// Include J coupling
		//JC = jmean + jstd * unormal_distribution(prandgen->prng);
		JC = prandgen->uniformRand(jmean-0.5*jwidth, jmean+0.5*jwidth);

		// Check if the spins are excited by detection pulses for the current field orientation
		if ( (detExcitA >= limit) && (detExcitB >= limit) ) {
			ampl += (detExcitA + detExcitB);
			// Compute the excitation efficiency of the spin A by pump pulse
			pumpExcitA = 0;
			for (size_t k = 0; k < compA; k++) pumpExcitA += static_cast<double>(spinA.intensity[k]) * pulses.pumpExcitation( freqA[k] );
			pumpExcitA /= statesA;
			pumpExcitA *= damping;
			// Compute the excitation efficiency of the spin B by pump pulse
			pumpExcitB = 0;
			for (size_t k = 0; k < compB; k++) pumpExcitB += static_cast<double>(spinB.intensity[k]) * pulses.pumpExcitation( freqB[k] );
			pumpExcitB /= statesB;
			pumpExcitB *= damping;
			// Check if the spins are excited by the pump pulse for the current field orientation
			if ( (pumpExcitA >= limit) || (pumpExcitB >= limit) ) {
				// Compute the cosine of the dipolar angle
				//dipProj = pow( magneticField.zprojection(n), 2);
				// Calculate the D constant in the dipolar frequency expression:
				D = 52.04 * gA * gB / pow(2.0023, 2);
				dipFreq =  2*PI * D * (1.0 - 3.0*dipProj) / pow(distance,3) + 2*PI * JC;
				for (size_t t = 0; t < nTime; ++t) {
					calcSignal[t] += (detExcitA*pumpExcitB + detExcitB*pumpExcitA) * (1.0 - cos(dipFreq*timepoints[t]));
					//calcSignal[t] += (detExcitA*(1-pumpExcitA) * pumpExcitB*(1-detExcitB) + detExcitB*(1-pumpExcitB) * pumpExcitA*(1-detExcitA)) * (1.0 - cos(dipFreq*timepoints[t]));
				}
			}
		}	
		else if ( (detExcitA >= limit) && (detExcitB < limit) ) {
			ampl += detExcitA;
			// Compute the excitation efficiency of the spin A by pump pulse
			pumpExcitA = 0;
			for (size_t k = 0; k < compA; k++) pumpExcitA += static_cast<double>(spinA.intensity[k]) * pulses.pumpExcitation( freqA[k] );
			pumpExcitA /= statesA;
			pumpExcitA *= damping;
			// Compute the excitation efficiency of the spin B by pump pulse
			pumpExcitB = 0;
			for (size_t k = 0; k < compB; k++) pumpExcitB += static_cast<double>(spinB.intensity[k]) * pulses.pumpExcitation( freqB[k] );
			pumpExcitB /= statesB;
			pumpExcitB *= damping;
			// Check if the spins are excited by the pump pulse for the current field orientation
			if (pumpExcitB >= limit) {
				// Compute the cosine of the dipolar angle
				//dipProj = pow( magneticField.zprojection(n), 2);
				// Calculate the D constant in the dipolar frequency expression:
				D = 52.04 * gA * gB / pow(2.0023, 2);
				dipFreq =  2*PI * D * (1.0 - 3.0*dipProj) / pow(distance,3) + 2*PI * JC;
				for (size_t t = 0; t < nTime; ++t) {
					calcSignal[t] += detExcitA*pumpExcitB * (1.0 - cos(dipFreq*timepoints[t])); 
					//calcSignal[t] += detExcitA*(1-pumpExcitA)*pumpExcitB * (1.0 - cos(dipFreq*timepoints[t]));
				}
			}
		}
		else if ( (detExcitA < limit) && (detExcitB >= limit) ) {
			ampl += detExcitB;
			// Compute the excitation efficiency of the spin A by pump pulse
			pumpExcitA = 0;
			for (size_t k = 0; k < compA; k++) pumpExcitA += static_cast<double>(spinA.intensity[k]) * pulses.pumpExcitation( freqA[k] );
			pumpExcitA /= statesA;
			pumpExcitA *= damping;
			// Compute the excitation efficiency of the spin B by pump pulse
			pumpExcitB = 0;
			for (size_t k = 0; k < compB; k++) pumpExcitB += static_cast<double>(spinB.intensity[k]) * pulses.pumpExcitation( freqB[k] );
			pumpExcitB /= statesB;
			pumpExcitB *= damping;
			// Check if the spins are excited by the pump pulse for the current field orientation
			if ( pumpExcitA >= limit) {
				// Compute the cosine of the dipolar angle
				//dipProj = pow( magneticField.zprojection(n), 2);
				// Calculate the D constant in the dipolar frequency expression:
				D = 52.04 * gA * gB / pow(2.0023, 2);
				dipFreq = 2*PI * D * (1.0 - 3.0*dipProj) / pow(distance,3) + 2*PI * JC;
				for (size_t t = 0; t < nTime; ++t) {
					calcSignal[t] += detExcitB*pumpExcitA * (1.0 - cos(dipFreq*timepoints[t]));
					//calcSignal[t] += detExcitB*(1-pumpExcitB)*pumpExcitA * (1.0 - cos(dipFreq*timepoints[t])); 
				}
			}
		}

		// Refresh some variables
		delete [] projectionA;
		delete [] projectionB;
		delete [] freqA;
		delete [] freqB;
		projectionA = nullptr;
		projectionB = nullptr;
		freqA = nullptr;
		freqB = nullptr;
		//delete pFieldRotation;
		delete pFrameRotation;
	} // end of the conformer's loop
			
	// Calculate a normalized timetrace and save it to the vector
	vector<double> calcSignal_formated; calcSignal_formated.reserve(nTime);
	double norm = 1 / ampl;
	for (size_t t = 0; t < nTime; ++t) { 
		calcSignal[t] = (ampl - calcSignal[t]) * norm;
		calcSignal_formated.push_back( calcSignal[t] );
	}
	
	// Delete dynamically allocated objects
	delete prandgen;
	delete [] calcSignal; calcSignal = nullptr;
	
	return calcSignal_formated;
}

double Enviroment::rmsd(vector<double> const& X, vector<double> const& expY, vector<double> const& calcY) const
{
	size_t const nTime = X.size();
	double rmsd(0);
	for (size_t t = 0; t < nTime; ++t) { rmsd += pow( (expY[t] - calcY[t]), 2 ); }
	rmsd = sqrt( rmsd / nTime );
	return rmsd;
}

double Enviroment::rmsd_PearsonCC(vector<double> const& X, vector<double> const& expY, 
	                              vector<double> const& calcY, double const& expYMean, 
	                              double const& expYVar) const
{
	size_t const nTime = X.size();
	double calcYMean(0), calcYVar(0), rmsd(0), lcc(0), rmsd_lcc(0);
	
	calcYMean = accumulate(calcY.begin(), calcY.end(), 0.0) / (double) nTime;
	
	for (size_t t = 0; t < nTime; ++t) { 
		rmsd += pow( (expY[t] - calcY[t]), 2 );
		lcc += (calcY[t] - calcYMean) * (expY[t] - expYMean);
		calcYVar += pow( (calcY[t] - calcYMean), 2 );
	}
	lcc /= sqrt(calcYVar * expYVar);
	rmsd = sqrt( rmsd / nTime );
	rmsd_lcc = rmsd / lcc;
	
	return rmsd_lcc;
}

double Enviroment::PearsonCC(vector<double> const& X, 
							 vector<double> const& expY,
							 vector<double> const& calcY,
							 double const& expYMean, 
							 double const& expYVar) const
{
	size_t const nTime = X.size();
	double calcYMean(0), calcYVar(0), lcc(0);
	
	calcYMean = accumulate(calcY.begin(), calcY.end(), 0.0) / (double) nTime;
	
	for (size_t t = 0; t < nTime; ++t) { 
		lcc += (calcY[t] - calcYMean) * (expY[t] - expYMean);
		calcYVar += pow( (calcY[t] - calcYMean), 2 );
	}
	lcc /= sqrt(calcYVar * expYVar);
	
	return (1.0-lcc);
}

void Enviroment::runScoring(Population& population, int const& scoreIndex, long const& nAvg) const
{
	parallel_for( size_t(0), population.size, [&] (size_t i) {
		
		vector<double> calcSignal; calcSignal.reserve(1000);
		double score(0);
		//double scoreWeight(0);

		for (size_t offset = 0; offset < nof; ++offset) {
			calcSignal = computeTimetrace(population.chromosomes[i].genes,expTime[offset],*pspinA,*pspinB,magneticFields[offset],pulseSequences[offset],nAvg,offset);
			/*
			scoreWeight = expSignal[0].size() / ( (1.0 - expMean[offset]) * expSignal[offset].size() );
			if (scoreIndex == 0) score += scoreWeight * rmsd(expTime[offset],expSignal[offset],calcSignal);
			if (scoreIndex == 1) score += scoreWeight * rmsd_PearsonCC(expTime[offset],expSignal[offset],calcSignal,expMean[offset],expVariance[offset]);
			if (scoreIndex == 2) score += scoreWeight * PearsonCC(expTime[offset],expSignal[offset],calcSignal,expMean[offset],expVariance[offset]);
			*/
			if (scoreIndex == 0) score += rmsd(expTime[offset],expSignal[offset],calcSignal);
			if (scoreIndex == 1) score += rmsd_PearsonCC(expTime[offset],expSignal[offset],calcSignal,expMean[offset],expVariance[offset]);
			if (scoreIndex == 2) score += PearsonCC(expTime[offset],expSignal[offset],calcSignal,expMean[offset],expVariance[offset]);
			
			calcSignal.clear();
		}
		population.chromosomes[i].fitness = score;
	});
}

void Enviroment::recordScore(double const& score, 
							 string const& file_root, 
							 size_t const& numOfPop) const
{
	ostringstream filename;
	filename << file_root << "score.dat";
	fstream file;
	file.open(filename.str().c_str(), ios::out|ios::app);
	file << setprecision(6) << setw(10) << numOfPop << setw(20) << score << endl;
	file.close();
}

void Enviroment::recordParameters(vector<double> const& parameters, 
								  string const& file_root) const
{
	ostringstream filename;
	filename << file_root << "param.dat";
	fstream file;
	file.open(filename.str().c_str(), ios::out|ios::app);
	int np(-1);
	for (size_t i = 0; i < pswitches.size(); ++i) {
		if ( pswitches[i] ) {
			if ((i == 0) || (i == 1) || (i == 12) || (i == 13) || (i >= 24)) 
				file << setprecision(4) << setw(12) << parameters[++np]; 
			else
				file << setprecision(4) << setw(12) << parameters[++np] * 180/PI;
		}
	}
	file << endl;
	file.close();
}

void Enviroment::recordFit(vector<double> const& parameters,
						   string const& file_root, 
						   size_t const& numOfPop) const
{
	parallel_for( size_t(0), nof, [&] (size_t offset) {
		
		// Calculate the timetrace
		size_t nTime = expTime[offset].size();
		vector<double> calcSignal; calcSignal.reserve(nTime);
		calcSignal = computeTimetrace(parameters,expTime[offset],*pspinA,*pspinB,magneticFields[offset],pulseSequences[offset],magneticFields[offset].N,offset);
		
		// Save the the timetrace to the file
		ostringstream filename;
		filename << file_root << "fit" << numOfPop << "_" << offset << ".dat";
		ofstream file(filename.str().c_str());
		if (file.is_open()) {
			for (size_t t = 0; t < nTime; ++t) { 
				file << setprecision(6) << setw(20) << expTime[offset][t] << setw(20) << expSignal[offset][t] << setw(20) << calcSignal[t] << endl; 
			}
			file.close();
		}

		calcSignal.clear();
		filename.str(string()); filename.clear();
	});
}

void Enviroment::recordSpectrum(vector<double> const& parameters, string const& file_root) const
{
	// Set the file in which the the spectrum will be recorded
	ostringstream filename;
	filename << file_root << "spectrum.dat";
	fstream file;
	file.open(filename.str().c_str(), ios::out|ios::app);

	// Initialize some variables
	RandGens* prandgen = new RandGens;	
	size_t const compA = pspinA->ncomp;	
	size_t const compB = pspinB->ncomp;
	double* projectionA(nullptr);
	double* freqA(nullptr);
	double* freqB(nullptr);
	double dipProj(0);
	double detExcitA(0), detExcitB(0), pumpExcitA(0), pumpExcitB(0);
	
	// Determine the lowest and the highest frequencies
	double freqMax, freqMin;
	for (long n = 0; n < magneticFields[0].N; n++) {
		// Select one of projections of magnetic field
		projectionA = magneticFields[0].projection(n);
		// Compute the resonance frequences of the spin A
		freqA = pspinA->resonanceFrequencies(projectionA, magneticFields[0].value, prandgen);
		// Compute the resonance frequence of the spin B
		freqB = pspinB->resonanceFrequencies(projectionA, magneticFields[0].value, prandgen);
		// Search for the lowest and the highest frequencies
		if (n == 0) {
			freqMin = freqA[0];
			freqMax = freqA[0];
		}
		for (size_t k = 0; k < compA; k++) {
			if (freqA[k] < freqMin) freqMin = freqA[k];
			if (freqA[k] > freqMax) freqMax = freqA[k];
		}
		for (size_t k = 0; k < compB; k++) {
			if (freqB[k] < freqMin) freqMin = freqB[k];
			if (freqB[k] > freqMax) freqMax = freqB[k];
		}
	}

	// Slightly increase the frequency range
	freqMin -= 0.100;
	freqMax += 0.100; 
	
	// Create the axes of the spectrum
	double df = 0.001;                                 // 1 MHz step
	freqMin = ceil(freqMin/df) * df;                  // round to the 0.001 precision
	freqMax = ceil(freqMax/df) * df;                  // round to the 0.001 precision
	int nf = static_cast<int>( (freqMax-freqMin)/df ); // number of intervals 
	vector<double> fAxis; fAxis.reserve(nf);           // x-axis
	vector<double> specAxis; specAxis.reserve(nf);	   // y-axis
	vector<double> detAxis; detAxis.reserve(nf);	   // y1-axis
	vector<double> pumpAxis; pumpAxis.reserve(nf);	   // y2-axis
	for (int i = 0; i < nf; i++) {
		fAxis.push_back( freqMin + 0.5*df*(2*i+1) );   // initialize x-axis
		specAxis.push_back( 0.0 );                     // initialize y-axis
		detAxis.push_back( 0.0 );                      // initialize y1-axis
		pumpAxis.push_back( 0.0 );                     // initialize y2-axis
	}

	// Correlate frequencies with their indeces
	int indexMin = static_cast<int>( freqMin/df ) + 1; 
	int index;

	// Calculate the spectrum
	for (long n = 0; n < magneticFields[0].N; n++) {
		// Select one of projections of magnetic field
		projectionA = magneticFields[0].projection(n);
		// Compute the resonance frequences of the spin A
		freqA = pspinA->resonanceFrequencies(projectionA, magneticFields[0].value, prandgen);
		// Compute the resonance frequence of the spin B
		freqB = pspinB->resonanceFrequencies(projectionA, magneticFields[0].value, prandgen);
		// Append the calculated frequencies to the spectrum
		for (size_t k = 0; k < compA; k++) {
			index = static_cast<int>( freqA[k]/df ) - indexMin;
			specAxis[index] += static_cast<double>(pspinA->intensity[k]) / static_cast<double>(pspinA->nstates);
			// Compute the excitation efficiencies
			detExcitA = static_cast<double>(pspinA->intensity[k]) * pulseSequences[0].detExcitation( freqA[k] );
			pumpExcitA = static_cast<double>(pspinA->intensity[k]) * pulseSequences[0].pumpExcitation( freqA[k] );
			detAxis[index] += detExcitA;
			pumpAxis[index] += pumpExcitA;
		}
		for (size_t k = 0; k < compB; k++) {
			index = static_cast<int>( freqB[k]/df ) - indexMin;
			specAxis[index] += static_cast<double>(pspinB->intensity[k]) / static_cast<double>(pspinB->nstates);
			// Compute the excitation efficiencies
			detExcitB = static_cast<double>(pspinB->intensity[k]) * pulseSequences[0].detExcitation( freqB[k] );
			pumpExcitB = static_cast<double>(pspinB->intensity[k]) * pulseSequences[0].pumpExcitation( freqB[k] );
			detAxis[index] += detExcitB;
			pumpAxis[index] += pumpExcitB;
		}
		// Refresh some variables
		delete [] projectionA;
		delete [] freqA;
		delete [] freqB;
		projectionA = nullptr;
		freqA = nullptr;
		freqB = nullptr;
	}
	
	// Normilize the spectrum
	double specAxisMax, detAxisMax, pumpAxisMax;
	specAxisMax = *max_element(specAxis.begin(), specAxis.end());
	detAxisMax  = *max_element(detAxis.begin(), detAxis.end());
	pumpAxisMax = *max_element(pumpAxis.begin(), pumpAxis.end());

	// Write the spectrum into the file
	for (int i = 0; i < nf; i++) {
		file << setprecision(8) << setw(20) << fAxis[i] << setw(20) << specAxis[i]/specAxisMax << setw(20) << detAxis[i]/detAxisMax << setw(20) << pumpAxis[i]/pumpAxisMax << endl;
	}
	file.close();

	// Delete dynamically allocated objects
	delete prandgen;
}

void Enviroment::recordSymmetricSolutions(vector<double> const& parameters, vector<double> const& ranges, int const& scoreIndex, string const& file_root) const
{
	// Set the optimized angles
	int count(0);
	double theta(0), phi(0), alpha(0), betta(0), gamma(0);
	for (size_t i = 0; i < pswitches.size(); ++i) {
		if ( pswitches[i] ) {
			if (i == 2)  theta = parameters[count]; 
			if (i == 4)  phi   = parameters[count];
			if (i == 6)  alpha = parameters[count];
			if (i == 8)  betta = parameters[count];
			if (i == 10) gamma = parameters[count];
			++count;
		}
	}

	// Calculate the orientation of the radius vector
	double radiusvector[3] = {cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)};

	// Calculate the rotation matrix from frame A to frame B
	EulerRotation const* pRframe = new EulerRotation(alpha,betta,gamma);

	// Create the 180° rotation matrices
	EulerRotation const* pRx = new EulerRotation(0.0, PI, 0.0);
	EulerRotation const* pRy = new EulerRotation(PI, PI, 0.0);
	EulerRotation const* pRz = new EulerRotation(PI, 0.0, 0.0);

	// Create the array of symmetry-related solutions
	size_t const local_size = 16;
	double all_angles[local_size][5];
	for (size_t i = 0; i < local_size; ++i) {
		all_angles[i][0] = theta;
		all_angles[i][1] = phi;
		all_angles[i][2] = alpha;
		all_angles[i][3] = betta;
		all_angles[i][4] = gamma;
	}
	
	// Initialize some supplementary variables
	double theta_new(0), phi_new(0), alpha_new(0), betta_new(0), gamma_new(0);
	double radiusvector_new[3];
	double Rframe_new[3][3];
	int solutionNumber(0);

	// Rotation of spin A about gx
	//----------------------------
	pRx->multiplyByVector(radiusvector, true, radiusvector_new);
	// Convert the radius vector into the angular parameters
	theta_new = acos( radiusvector_new[2] );
    if ((theta_new == 0) || (theta_new == PI)) phi_new = 0;
    else {
        if ( radiusvector_new[0] == 0 ) {
            if ( radiusvector_new[1] == sin(theta_new) ) phi_new = 0.5*PI;
            else                                         phi_new = -0.5*PI;
		}
        else phi_new = atan2( radiusvector_new[1], radiusvector_new[0] );
	}
	if (phi_new < 0) phi_new += 2*PI;
	// Save new parameters
	++solutionNumber;
	all_angles[solutionNumber][0] = theta_new;
	all_angles[solutionNumber][1] = phi_new;

	// Rotation of spin A about gy
	//----------------------------
	pRy->multiplyByVector(radiusvector, true, radiusvector_new);
	// Convert the radius vector into the angular parameters
	theta_new = acos( radiusvector_new[2] );
    if ((theta_new == 0) || (theta_new == PI)) phi_new = 0;
    else {
        if ( radiusvector_new[0] == 0 ) {
            if ( radiusvector_new[1] == sin(theta_new) ) phi_new = 0.5*PI;
            else                                         phi_new = -0.5*PI;
		}
        else phi_new = atan2( radiusvector_new[1], radiusvector_new[0] );
	}
	if (phi_new < 0) phi_new += 2*PI;
	// Save new parameters
	++solutionNumber;
	all_angles[solutionNumber][0] = theta_new;
	all_angles[solutionNumber][1] = phi_new;

	// Rotation of spin A about gz
	//----------------------------
	pRz->multiplyByVector(radiusvector, true, radiusvector_new);
	// Convert the radius vector into the angular parameters
	theta_new = acos( radiusvector_new[2] );
    if ((theta_new == 0) || (theta_new == PI)) phi_new = 0;
    else {
        if ( radiusvector_new[0] == 0 ) {
            if ( radiusvector_new[1] == sin(theta_new) ) phi_new = 0.5*PI;
            else                                         phi_new = -0.5*PI;
		}
        else phi_new = atan2( radiusvector_new[1], radiusvector_new[0] );
	}
	if (phi_new < 0) phi_new += + 2*PI;
	// Save new parameters
	++solutionNumber;
	all_angles[solutionNumber][0] = theta_new;
	all_angles[solutionNumber][1] = phi_new;

	// Rotation of spin B about gx
	//----------------------------
	pRframe->multiplyByMatrix(pRx->R, false, Rframe_new);
	// Convert the radius vector into the angular parameters
	betta_new = acos( Rframe_new[2][2] );
    if (betta_new == 0) {
        gamma_new = 0;
        alpha_new = atan2(Rframe_new[1][0], Rframe_new[1][1]);
	}
    else {
        gamma_new = atan2(Rframe_new[2][0]/sin(betta_new), Rframe_new[2][1]/sin(betta_new));
        alpha_new = atan2(Rframe_new[0][2]/sin(betta_new), -Rframe_new[1][2]/sin(betta_new));
	}
    if (alpha_new < 0) alpha_new += 2*PI;
    if (gamma_new < 0) gamma_new += 2*PI;
	// Save new parameters
	++solutionNumber;
	all_angles[solutionNumber][2] = alpha_new;
	all_angles[solutionNumber][3] = betta_new;
	all_angles[solutionNumber][4] = gamma_new;

	// Rotation of spin B about gy
	//----------------------------
	pRframe->multiplyByMatrix(pRy->R, false, Rframe_new);
	// Convert the radius vector into the angular parameters
	betta_new = acos( Rframe_new[2][2] );
    if (betta_new == 0) {
        gamma_new = 0;
        alpha_new = atan2(Rframe_new[1][0], Rframe_new[1][1]);
	}
    else {
        gamma_new = atan2(Rframe_new[2][0]/sin(betta_new), Rframe_new[2][1]/sin(betta_new));
        alpha_new = atan2(Rframe_new[0][2]/sin(betta_new), -Rframe_new[1][2]/sin(betta_new));
	}
    if (alpha_new < 0) alpha_new += 2*PI;
    if (gamma_new < 0) gamma_new += 2*PI;
	// Save new parameters
	++solutionNumber;
	all_angles[solutionNumber][2] = alpha_new;
	all_angles[solutionNumber][3] = betta_new;
	all_angles[solutionNumber][4] = gamma_new;

	// Rotation of spin B about gz
	//----------------------------
	pRframe->multiplyByMatrix(pRz->R, false, Rframe_new);
	// Convert the radius vector into the angular parameters
	betta_new = acos( Rframe_new[2][2] );
    if (betta_new == 0) {
        gamma_new = 0;
        alpha_new = atan2(Rframe_new[1][0], Rframe_new[1][1]);
	}
    else {
        gamma_new = atan2(Rframe_new[2][0]/sin(betta_new), Rframe_new[2][1]/sin(betta_new));
        alpha_new = atan2(Rframe_new[0][2]/sin(betta_new), -Rframe_new[1][2]/sin(betta_new));
	}
    if (alpha_new < 0) alpha_new += 2*PI;
    if (gamma_new < 0) gamma_new += 2*PI;
	// Save new parameters
	++solutionNumber;
	all_angles[solutionNumber][2] = alpha_new;
	all_angles[solutionNumber][3] = betta_new;
	all_angles[solutionNumber][4] = gamma_new;

	// Combined rotations
	for (size_t i = 1; i <= 3; ++i) {
		for (size_t j = 4; j <= 6; ++j) {
			++solutionNumber;
			all_angles[solutionNumber][0] = all_angles[i][0];
			all_angles[solutionNumber][1] = all_angles[i][1];
			all_angles[solutionNumber][2] = all_angles[j][2];
			all_angles[solutionNumber][3] = all_angles[j][3];
			all_angles[solutionNumber][4] = all_angles[j][4];
		}
	}

	// Create a local population
	Population* local_population = new Population(local_size,parameters.size(),ranges,pswitches);
	// Substitute the genes of the created chromosomes by the optimized genes
	for (size_t i = 0; i < local_size; ++i) local_population->chromosomes[i].genes = parameters;
	// Introduce changes according to the symmetry-related solutions
	int count1(0);
	for (size_t i = 0; i < pswitches.size(); ++i) {
		if ( pswitches[i] ) {
			if (i == 2)  { for (size_t j = 0; j < local_size; ++j) local_population->chromosomes[j].genes[count1] = all_angles[j][0]; } 
			if (i == 4)  { for (size_t j = 0; j < local_size; ++j) local_population->chromosomes[j].genes[count1] = all_angles[j][1]; } 
			if (i == 6)  { for (size_t j = 0; j < local_size; ++j) local_population->chromosomes[j].genes[count1] = all_angles[j][2]; } 
			if (i == 8)  { for (size_t j = 0; j < local_size; ++j) local_population->chromosomes[j].genes[count1] = all_angles[j][3]; } 
			if (i == 10) { for (size_t j = 0; j < local_size; ++j) local_population->chromosomes[j].genes[count1] = all_angles[j][4]; } 
			++count1;
		}
	}

	// Score the generated population
	runScoring(*local_population, scoreIndex, magneticFields[0].N);

	// Save to the file
	ostringstream filename;
	filename << file_root << "symmetricSolutions.dat";
	ofstream file(filename.str().c_str());
	if (file.is_open()) {
		int np(-1);
		for (size_t j = 0; j < local_size; ++j) {
			np = -1;
			for (size_t i = 0; i < pswitches.size(); ++i) {
				if ( pswitches[i] ) {
					if ((i == 0) || (i == 1) || (i == 12) || (i == 13) || (i >= 24)) 
						file << setprecision(4) << setw(12) << local_population->chromosomes[j].genes[++np]; 
					else
						file << setprecision(4) << setw(12) << local_population->chromosomes[j].genes[++np] * 180/PI;
				}
			}
			file << setprecision(6) << setw(12) << local_population->chromosomes[j].fitness << endl; 			
		}
		file.close();
	}

	// Delete dynamically allocated objects
	delete pRframe;
	delete pRx;
	delete pRy;
	delete pRz;
	delete local_population;
}

vector<double> Enviroment::computeModDepth(vector<double> const& parameters, SpinSystem const& spinA, SpinSystem const& spinB,
										   MagneticField const& magneticField, Pulses const& pulses) const
{		
	// Initialize the array for the modulation depth parameter
	double modDepths_abs[90];
	for (int i = 0; i < 90; ++i) modDepths_abs[i] = 0;
	vector<double> modDepths_norm; modDepths_norm.reserve(90);

	// Set the values of geometric parameters
	double dist1_mean,  dist1_width,  dist2_mean,  dist2_width;
	double theta1_mean, theta1_width, theta2_mean, theta2_width;
	double phi1_mean,   phi1_width,   phi2_mean,   phi2_width;
	double alpha1_mean, alpha1_width, alpha2_mean, alpha2_width;
	double betta1_mean, betta1_width, betta2_mean, betta2_width;
	double gamma1_mean, gamma1_width, gamma2_mean, gamma2_width;
	double ratio;
	int np(-1);
	if ( pswitches[0] ) {
		dist1_mean = parameters[++np];
		if ( pswitches[1] ) dist1_width = parameters[++np];
		else          dist1_width = 0;
	}
	else { dist1_mean = 0; dist1_width = 0; }
	if ( pswitches[2] ) {
		theta1_mean = parameters[++np];
		if ( pswitches[3] ) theta1_width = parameters[++np];
		else          theta1_width = 0;
	}
	else { theta1_mean = 0; theta1_width = 0; }
	if ( pswitches[4] ) {
		phi1_mean = parameters[++np];
		if ( pswitches[5] ) phi1_width = parameters[++np];
		else          phi1_width = 0;
	}
	else { phi1_mean = 0; phi1_width = 0; }
	if ( pswitches[6] ) {
		alpha1_mean = parameters[++np];
		if ( pswitches[7] ) alpha1_width = parameters[++np];
		else          alpha1_width = 0;
	}
	else { alpha1_mean = 0; alpha1_width = 0; }
	if ( pswitches[8] ) {
		betta1_mean = parameters[++np];
		if ( pswitches[9] ) betta1_width = parameters[++np];
		else          betta1_width = 0;
	}
	else { betta1_mean = 0; betta1_width = 0; }
	if ( pswitches[10] ) {
		gamma1_mean = parameters[++np];
		if ( pswitches[11] ) gamma1_width = parameters[++np];
		else           gamma1_width = 0;
	}
	else { gamma1_mean = 0; gamma1_width = 0; }
	if ( pswitches[12] ) {
		dist2_mean = parameters[++np];
		if ( pswitches[13] ) dist2_width = parameters[++np];
		else           dist2_width = dist1_width;
	}
	else { dist2_mean = dist1_mean; dist2_width = dist1_width; }
	if ( pswitches[14] ) {
		theta2_mean = parameters[++np];
		if ( pswitches[15] ) theta2_width = parameters[++np];
		else           theta2_width = theta1_width;
	}
	else { theta2_mean = theta1_mean; theta2_width = theta1_width; }
	if ( pswitches[16] ) {
		phi2_mean = parameters[++np];
		if ( pswitches[17] ) phi2_width = parameters[++np];
		else           phi2_width = phi1_width;
	}
	else { phi2_mean = phi1_mean; phi2_width = phi1_width; }
	if ( pswitches[18] ) {
		alpha2_mean = parameters[++np];
		if ( pswitches[19] ) alpha2_width = parameters[++np];
		else           alpha2_width = alpha1_width;
	}
	else { alpha2_mean = alpha1_mean; alpha2_width = alpha1_width; }
	if ( pswitches[20] ) {
		betta2_mean = parameters[++np];
		if ( pswitches[21] ) betta2_width = parameters[++np];
		else           betta2_width = betta1_width;
	}
	else { betta2_mean = betta1_mean; betta2_width = betta1_width; }
	if ( pswitches[22] ) {
		gamma2_mean = parameters[++np];
		if ( pswitches[23] ) gamma2_width = parameters[++np];
		else           gamma2_width = gamma1_width;
	}
	else { gamma2_mean = gamma1_mean; gamma2_width = gamma1_width;}
	if ( pswitches[24] ) ratio = parameters[++np];
	else           ratio = 1;

	////////////////////////////////////////////////////////////////////////////
	//            COMPUTE THE MODULATION DEPTH PARAMETER                      //
	////////////////////////////////////////////////////////////////////////////

	RandGens* prandgen = new RandGens;            // Initialize random generator
	std::normal_distribution<double> unormal_distribution(0,1); // Initialize a normal distribution
	long const nAvg = magneticField.N;            // Set the number of averages
	size_t const compA = spinA.ncomp;             // Set the number of spectral components of the spin A
	double const statesA = (double) spinA.nstates;// Set the number of microstates of the spin A
	size_t const compB = spinB.ncomp;             // Set the number of spectral components of the spin B
	double const statesB = (double) spinB.nstates;// Set the number of microstates of the spin B
	bool first_ensemble;
	double theta(0), phi(0), alpha(0), betta(0), gamma(0);
	double* projectionA(nullptr);
	double* projectionB(nullptr);
	double* freqA(nullptr);
	double* freqB(nullptr);
	double detExcitA(0), detExcitB(0), pumpExcitA(0), pumpExcitB(0);
	double modDepth(0);
	double dipAngle(0);
	size_t ndipAngle(0);

	// Itterate over the conformers
	for (long n = 0; n < nAvg; n++) {

		// Choose between 2 ensembles of parameters
		if (prandgen->uniformRand() <= ratio) first_ensemble = true;
		else first_ensemble = false;

		// Generate the set of geometrical parameters of the system
		/**/
		if (first_ensemble) {       // if the 1st ensemble of parameters was chosen
			if (theta1_width == 0) theta = theta1_mean;
			else                   theta = prandgen->uniformRand(theta1_mean-0.5*theta1_width, theta1_mean+0.5*theta1_width);
			if (phi1_width == 0)   phi = phi1_mean;
			else                   phi = prandgen->uniformRand(phi1_mean-0.5*phi1_width, phi1_mean+0.5*phi1_width);
			if (alpha1_width == 0) alpha = alpha1_mean;
			else                   alpha = prandgen->uniformRand(alpha1_mean-0.5*alpha1_width, alpha1_mean+0.5*alpha1_width);
			if (betta1_width == 0) betta = betta1_mean;
			else                   betta = prandgen->uniformRand(betta1_mean-0.5*betta1_width, betta1_mean+0.5*betta1_width);
			if (gamma1_width == 0) gamma = gamma1_mean;
			else                   gamma = prandgen->uniformRand(gamma1_mean-0.5*gamma1_width, gamma1_mean+0.5*gamma1_width);
		}
		else {                      // if the 2nd ensemble of parameters was chosen
			if (theta2_width == 0) theta = theta2_mean;
			else                   theta = prandgen->uniformRand(theta2_mean- 0.5*theta2_width, theta2_mean+0.5*theta2_width);
			if (phi2_width == 0)   phi = phi2_mean;
			else                   phi = prandgen->uniformRand(phi2_mean- 0.5*phi2_width, phi2_mean+0.5*phi2_width);
			if (alpha2_width == 0) alpha = alpha2_mean;
			else                   alpha = prandgen->uniformRand(alpha2_mean-0.5*alpha2_width, alpha2_mean+0.5*alpha2_width);
			if (betta2_width == 0) betta = betta2_mean;
			else                   betta = prandgen->uniformRand(betta2_mean-0.5*betta2_width, betta2_mean+0.5*betta2_width);
			if (gamma2_width == 0) gamma = gamma2_mean;
			else                   gamma = prandgen->uniformRand(gamma2_mean-0.5*gamma2_width, gamma2_mean+0.5*gamma2_width);
		}
        /**/
		/*
		if (first_ensemble) {       // if the 1st ensemble of parameters was chosen
			if (theta1_width == 0) theta = theta1_mean;
			else                   theta = theta1_mean + theta1_width * unormal_distribution(prandgen->prng);
			if (phi1_width == 0)   phi = phi1_mean;
			else                   phi = phi1_mean + phi1_width * unormal_distribution(prandgen->prng);
			if (alpha1_width == 0) alpha = alpha1_mean;
			else                   alpha = alpha1_mean + alpha1_width * unormal_distribution(prandgen->prng);
			if (betta1_width == 0) betta = betta1_mean;
			else                   betta = betta1_mean + betta1_width * unormal_distribution(prandgen->prng);
			if (gamma1_width == 0) gamma = gamma1_mean;
			else                   gamma = gamma1_mean + gamma1_width * unormal_distribution(prandgen->prng);
		}
		else {                      // if the 2nd ensemble of parameters was chosen
			if (theta2_width == 0) theta = theta2_mean;
			else                   theta = theta2_mean + theta2_width * unormal_distribution(prandgen->prng);
			if (phi2_width == 0)   phi = phi2_mean;
			else                   phi = phi2_mean + phi2_width * unormal_distribution(prandgen->prng);
			if (alpha2_width == 0) alpha = alpha2_mean;
			else                   alpha = alpha2_mean + alpha2_width * unormal_distribution(prandgen->prng);
			if (betta2_width == 0) betta = betta2_mean;
			else                   betta = betta2_mean + betta2_width * unormal_distribution(prandgen->prng);
			if (gamma2_width == 0) gamma = gamma2_mean;
			else                   gamma = gamma2_mean + gamma2_width * unormal_distribution(prandgen->prng);
		}
		*/
                        
		// Rotation matrix between the spin A and spin B coordinate systems
		EulerRotation const* pFrameRotation = new EulerRotation(alpha,betta,gamma);

		// Compute the orientation of the magnetic field in the frame of the spin A
		projectionA = magneticField.projection(n);
		// Compute the orientation of the magnetic field in the frame of the spin B
		projectionB = magneticField.frameB(*pFrameRotation, projectionA);

		// Compute the resonance frequences of the spin A
		freqA = spinA.resonanceFrequencies(projectionA, magneticField.value, prandgen);
		// Compute the resonance frequence of the spin B
		freqB = spinB.resonanceFrequencies(projectionB, magneticField.value, prandgen);

		// Compute the excitation efficiency of the spin A by detection pulses
		detExcitA = 0;
		for (size_t k = 0; k < compA; k++) detExcitA += static_cast<double>(spinA.intensity[k]) * pulses.detExcitation( freqA[k] );
		detExcitA /= statesA;
		// Compute the excitation efficiency of the spin B by detection pulses
		detExcitB = 0;
		for (size_t k = 0; k < compB; k++) detExcitB += static_cast<double>(spinB.intensity[k]) * pulses.detExcitation( freqB[k] );
		detExcitB /= statesB;
		// Compute the excitation efficiency of the spin A by pump pulse
		pumpExcitA = 0;
		if (detExcitB >= limit) {
			for (size_t k = 0; k < compA; k++) pumpExcitA += static_cast<double>(spinA.intensity[k]) * pulses.pumpExcitation( freqA[k] );
			pumpExcitA /= statesA;
		}
		// Compute the excitation efficiency of the spin B by pump pulse
		pumpExcitB = 0;
		if (detExcitA >= limit) {
			for (size_t k = 0; k < compB; k++) pumpExcitB += static_cast<double>(spinB.intensity[k]) * pulses.pumpExcitation( freqB[k] );
			pumpExcitB /= statesB;
		}
		
		// Compute the modulation depth parameter
		modDepth = (detExcitA * pumpExcitB + detExcitB * pumpExcitA);
		//modDepth = detExcitA*(1-pumpExcitA) * pumpExcitB*(1-detExcitB) + detExcitB*(1-pumpExcitB) * pumpExcitA*(1-detExcitA);

		// Compute the dipolar angle
		dipAngle = 180/PI * acos(projectionA[0] * cos(phi)*sin(theta) + projectionA[1] * sin(phi)*sin(theta) + projectionA[2] * cos(theta));
		if (dipAngle < 0)  dipAngle = -dipAngle;
		if (dipAngle > 90) dipAngle = 180 - dipAngle;
		
		// Compute the number of the angle
		if (dipAngle == 90) ndipAngle = 89;
		else ndipAngle = (size_t) floor( dipAngle );

		// Append the calculated modulation depth parameter to the corresponding dipolar angle
		modDepths_abs[ndipAngle] += modDepth;

		// Refresh some variables
		delete [] projectionA;
		delete [] projectionB;
		delete [] freqA;
		delete [] freqB;
		projectionA = nullptr;
		projectionB = nullptr;
		freqA = nullptr;
		freqB = nullptr;
		delete pFrameRotation;
	} // end of the conformer's loop

	// Delete dynamically allocated objects
	delete prandgen;

	// Normilize the angular dependence of the modulation depth parameter to 1
	double modDepths_int(0);
	for (int i = 0; i < 90; ++i) modDepths_int += modDepths_abs[i];
	modDepths_int *= PI/180;
	for (int i = 0; i < 90; ++i) modDepths_norm.push_back( modDepths_abs[i] / modDepths_int );

	return modDepths_norm;
}

void Enviroment::recordModDepth(vector<double> const& parameters, string const& file_root) const
{
	parallel_for( size_t(0), nof, [&] (size_t offset) {

		vector<int> angles; angles.reserve(90);
		for (int i = 0; i < 90; ++i) angles.push_back( i+1 );

		vector<double> lambdas; lambdas.reserve(90);
		lambdas = computeModDepth(parameters,*pspinA,*pspinB,magneticFields[offset],pulseSequences[offset]);
		
		// Save the the timetrace to the file
		ostringstream filename;
		filename << file_root << "lambda" << offset << ".dat";
		ofstream file(filename.str().c_str());
		if (file.is_open()) {
			for (size_t i = 0; i < 90; ++i) file << setprecision(6) << setw(15) << angles[i] << setw(15) << lambdas[i] << endl; 
			file.close();
		}
	});
}

void Enviroment::recordErrorProfile(vector<double> const& parameters, vector<int> const& variables, 
									vector<double> const& ranges, string const& file_root) const
{
	size_t const local_size = 2048;
	size_t const local_length = parameters.size();
	
	// Define which parameters are optimized 
	vector<int> paramIndices; paramIndices.reserve(variables.size());
	int count(-1);
	for (size_t i = 0; i < pswitches.size(); ++i) {
		if ( pswitches[i] ) ++count;
		for (size_t j = 0; j < variables.size(); ++j) { 
			if (i == variables[j]) { paramIndices.push_back(count); break; }
		}
	}
	
	// Create a new bounds
	vector<double> local_bounds; local_bounds.reserve(2*local_length);
	for (size_t i = 0; i < local_length; ++i) {
		local_bounds.push_back( parameters[i] );
		local_bounds.push_back( parameters[i] + 1.e-12 );
	}
	for (size_t i = 0; i < paramIndices.size(); ++i) {
		local_bounds[ 2*paramIndices[i] ]   = ranges[ 2*i ];
		local_bounds[ 2*paramIndices[i]+1 ] = ranges[ 2*i+1 ];
	}
	//cout << "Bounds for the error profile:" << endl;
	//for (size_t i = 0; i < local_length; ++i) cout << setw(15) << local_bounds[2*i] << setw(15) << local_bounds[2*i+1] << endl;
	
	// Create a population
	Population* local_population = new Population(local_size,local_length,local_bounds,pswitches);
	
	// Run scoring
	int errorScore = 0; // Use RMSD to calculate error profile
	runScoring(*local_population, errorScore, magneticFields[0].N);

	// Range chromosomes
	local_population->range();

	// Save to the file
	ostringstream filename;
	filename << file_root << "errorprofile";
	for (size_t i = 0; i < variables.size(); ++i) filename << "_" << variables[i];
	filename << ".dat";
	ofstream file(filename.str().c_str());
	if (file.is_open()) {
		for (size_t i = 0; i < local_size; ++i) {
			for (size_t j = 0; j < variables.size(); ++j) {
				if ((variables[j] == 0) || (variables[j] == 1) || (variables[j] == 12) || (variables[j] == 13) || (variables[j] >= 24))
					file << setprecision(6) << setw(15) << local_population->chromosomes[i].genes[paramIndices[j]];
				else
					file << setprecision(6) << setw(15) << local_population->chromosomes[i].genes[paramIndices[j]] * 180/PI;
			}
			file << setprecision(6) << setw(15) << local_population->chromosomes[i].fitness << endl; 
		}
		file.close();
	}
	
	delete local_population;
}

