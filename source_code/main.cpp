#ifdef WIN32
	#pragma comment(lib, "libconfig")
#endif
#include "libconfig.h"
#include "SpinSystem.h"
#include "MagneticField.h"
#include "Pulses.h"
#include "Enviroment.h"

#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::setw;
#include <vector>
using std::vector;
#include <iterator>
#include <string>
using std::string;
using std::getline;
#include <fstream>
using std::ifstream;
#include <sstream>
using std::stringstream;
#include <numeric>
using std::accumulate;
#ifdef __gnu_linux__
	#include <sys/time.h>
#endif


int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	//                           READ INI FILE                                //
	////////////////////////////////////////////////////////////////////////////
  
	// Check that the number of arguments equals 2	
	if ( argc != 2 ) {
		std::cerr << "Incorrect number of arguments!" << endl;
        return EXIT_FAILURE;
	}	
	
	// Read input data from the config file
	config_t cfg, *cf;
	cf = &cfg;
    config_init(cf);
	if (!config_read_file(cf, argv[1])) {
		std::cerr <<  config_error_file(cf) << ": " << config_error_line(cf) << " - " << config_error_text(cf) << endl;
		config_destroy(cf);
		return EXIT_FAILURE;
	}
	// Enable the auto convertion of data types
	config_set_auto_convert(cf, true);
	
	// Experimental parameters--------------------------------------------------
	
	const config_setting_t *experimentals = config_lookup(cf, "experimentals");
    size_t nOffset = config_setting_length(experimentals);
    cout << endl << "Number of PELDOR timetraces: " << nOffset << endl;
	
	const config_setting_t *set;
	const char *path = NULL;
	string dataline;
	double timeval, signalval;
	vector<double> t_expTime;
	vector<double> t_expSignal;
	vector<vector<double>> expTime;   expTime.reserve(nOffset);
	vector<vector<double>> expSignal; expSignal.reserve(nOffset);
	double t_detPiLength, t_detPiHalfLength, t_pumpPiLength, t_detFreq, t_pumpFreq, t_valueField;
	vector<double> detPiLength;     detPiLength.reserve(nOffset);
	vector<double> detPiHalfLength; detPiHalfLength.reserve(nOffset);
	vector<double> pumpPiLength;    pumpPiLength.reserve(nOffset);
	vector<double> detFreq;         detFreq.reserve(nOffset);
	vector<double> pumpFreq;        pumpFreq.reserve(nOffset);
	vector<double> valueField;      valueField.reserve(nOffset);
	
	for (unsigned int i = 0; i < nOffset; ++i) {
		set = config_setting_get_elem(experimentals, i);
		// Read experimental data
		path = config_setting_get_string( config_setting_get_member(set, "filename") );
		ifstream datafile(path, std::ios::in);
		if (! datafile.is_open()) {
			std::cerr << "Could not open the data file!" << endl;
			return EXIT_FAILURE;
		}
		else {
			while (getline(datafile, dataline)) {
				stringstream datastream(dataline);
				datastream >> timeval >> signalval;
				t_expTime.push_back(timeval);
				t_expSignal.push_back(signalval);
			}
			expTime.push_back(t_expTime);
			expSignal.push_back(t_expSignal);
			t_expTime.clear();
			t_expSignal.clear();
		}
		// Read experimental parameters
		t_detPiLength = config_setting_get_float( config_setting_get_member(set, "detPiLength") );
		t_detPiHalfLength = config_setting_get_float( config_setting_get_member(set, "detPiHalfLength") );
		t_pumpPiLength = config_setting_get_float( config_setting_get_member(set, "pumpPiLength") );
		t_detFreq = config_setting_get_float( config_setting_get_member(set, "detFreq") );
		t_pumpFreq = config_setting_get_float( config_setting_get_member(set, "pumpFreq") );
		t_valueField = config_setting_get_float( config_setting_get_member(set, "magnField") );
		detPiLength.push_back(t_detPiLength);
		detPiHalfLength.push_back(t_detPiHalfLength);
		pumpPiLength.push_back(t_pumpPiLength);
		detFreq.push_back(t_detFreq);
		pumpFreq.push_back(t_pumpFreq);
		valueField.push_back(t_valueField);
	}
	
	// Spectroscopic parameters-------------------------------------------------

	double spinA_g[3];
	vector<int> spinA_n;
	vector<double> spinA_I;
	vector<double> spinA_A;
	const config_setting_t *spinA = config_lookup(cf, "spinA");
	for (size_t i = 0; i < 3; ++i) {
		spinA_g[i] = config_setting_get_float_elem( config_setting_get_member(spinA, "g"), i );
	}
	size_t sizeA = config_setting_length( config_setting_get_member(spinA, "n") );
	for (size_t i = 0; i < sizeA; ++i) {
		spinA_n.push_back( config_setting_get_int_elem( config_setting_get_member(spinA, "n"), i ) );
		spinA_I.push_back( config_setting_get_float_elem( config_setting_get_member(spinA, "I"), i ) );
		spinA_A.push_back( config_setting_get_float_elem( config_setting_get_member(spinA, "A"), 3*i+0 ) );
		spinA_A.push_back( config_setting_get_float_elem( config_setting_get_member(spinA, "A"), 3*i+1 ) );
		spinA_A.push_back( config_setting_get_float_elem( config_setting_get_member(spinA, "A"), 3*i+2 ) );
	}
	double spinA_gStrain[3] = {0,0,0};
	double spinA_AStrain[3] = {0,0,0};
	double spinA_lwpp(0);
	if ( config_setting_length(config_setting_get_member(spinA, "gStrain")) != 0 ) {
		for (size_t i = 0; i < 3; ++i) spinA_gStrain[i] = config_setting_get_float_elem( config_setting_get_member(spinA, "gStrain"), i );
	}
	if ( config_setting_length(config_setting_get_member(spinA, "AStrain")) != 0 ) {
		for (size_t i = 0; i < 3; ++i) spinA_AStrain[i] = config_setting_get_float_elem( config_setting_get_member(spinA, "AStrain"), i );
	}
	spinA_lwpp = config_setting_get_float( config_setting_get_member(spinA, "lwpp") );
	
	double spinB_g[3];
	vector<int> spinB_n;
	vector<double> spinB_I;
	vector<double> spinB_A;
	const config_setting_t *spinB = config_lookup(cf, "spinB");
	for (size_t i = 0; i < 3; ++i) {
		spinB_g[i] = config_setting_get_float_elem( config_setting_get_member(spinB, "g"), i );
	}
	size_t sizeB = config_setting_length( config_setting_get_member(spinB, "n") );
	for (size_t i = 0; i < sizeB; ++i) {
		spinB_n.push_back( config_setting_get_int_elem( config_setting_get_member(spinB, "n"), i ) );
		spinB_I.push_back( config_setting_get_float_elem( config_setting_get_member(spinB, "I"), i ) );
		spinB_A.push_back( config_setting_get_float_elem( config_setting_get_member(spinB, "A"), 3*i+0 ) );
		spinB_A.push_back( config_setting_get_float_elem( config_setting_get_member(spinB, "A"), 3*i+1 ) );
		spinB_A.push_back( config_setting_get_float_elem( config_setting_get_member(spinB, "A"), 3*i+2 ) );
	}
	double spinB_gStrain[3] = {0,0,0};
	double spinB_AStrain[3] = {0,0,0};
	double spinB_lwpp(0);
	if ( config_setting_length(config_setting_get_member(spinB, "gStrain")) != 0 ) {
		for (size_t i = 0; i < 3; ++i) spinB_gStrain[i] = config_setting_get_float_elem( config_setting_get_member(spinB, "gStrain"), i );
	}
	if ( config_setting_length(config_setting_get_member(spinB, "AStrain")) != 0 ) {
		for (size_t i = 0; i < 3; ++i) spinB_AStrain[i] = config_setting_get_float_elem( config_setting_get_member(spinB, "AStrain"), i );
	}
	spinB_lwpp = config_setting_get_float( config_setting_get_member(spinB, "lwpp") );

	
	// Fitting parameters-------------------------------------------------
	
	// Set the path for the data to be saved
	const char *pfileroot = NULL;
	config_lookup_string(cf, "OutputDirectory", &pfileroot);
	string fileroot(pfileroot);

	// Set the fitting parameters
	const config_setting_t *fitParam = config_lookup(cf, "parameters");
	int nFitParam = config_setting_length(fitParam);
	if ((nFitParam != 28) && (nFitParam != 27+nOffset)) std::cerr << "Wrong number of fitting parameters!" << endl;
	// Set the fitting parameters
	vector<bool> paramSwitches; paramSwitches.reserve(nFitParam);
	const config_setting_t *oneParam;
	int chrLength(0);
	for (int i = 0; i < nFitParam; ++i) {
		// Read one parameter
		oneParam = config_setting_get_elem(fitParam, i);
		// Check if this parameter is enabled
		if ( config_setting_get_int(config_setting_get_member(oneParam,"opt")) ) {
			paramSwitches.push_back(true);
			// Increase the number of fitting parameters by 1
			++chrLength;
		}
		else {
			paramSwitches.push_back(false);
		} 
	}
	// Set the ranges for the fitting parameters
	double lbound, hbound;
	vector<double> allbounds; allbounds.reserve(2*chrLength);
	for (int i = 0; i < nFitParam; ++i) {
		// Check if this parameter is enabled
		if ( paramSwitches[i] ) {
			// Read one parameter
			oneParam = config_setting_get_elem(fitParam, i);
			// Read the range for this parameter
			lbound = config_setting_get_float_elem( config_setting_get_member(oneParam,"range"), 0 );
			hbound = config_setting_get_float_elem( config_setting_get_member(oneParam,"range"), 1 );
			// Save the ranges to the vector allbounds
			if ((i == 0) || (i == 1) || (i == 12) || (i == 13) || (i >= 24)) {
				allbounds.push_back(lbound);
				allbounds.push_back(hbound);
			}
			else {
				allbounds.push_back(lbound * PI/180);
				allbounds.push_back(hbound * PI/180);
			}
		}	
	}

	// Scoring function
	int scoreIndex;	
	config_lookup_int(cf, "scoreIndex", &scoreIndex);

	// Additional capabilities--------------------------------------------

	// Record the EPR spectrum of the spin pair
	int Spectrum(0);
	config_lookup_bool(cf, "Spec", &Spectrum);

	// Create symmetry-related solutions
	int SymSol(0);
	config_lookup_bool(cf, "SymSol", &SymSol);

	// Recording the angular dependence of the modulation depth parameter
	int modDepth(0);
	config_lookup_bool(cf, "ModDepth", &modDepth);

	// Recording the error profile for the fitting parameters
	int errorProfile(0);
	config_lookup_bool(cf, "ErrorProfile", &errorProfile);
	
	vector<int> errorVar;    // Variables of the error profile 
	vector<double> errorBounds;	// Ranges of variables of the error profile 
	if (errorProfile) {
		const config_setting_t *errorParam = config_lookup(cf, "ErrorProfileVar");
		// Number of variables
		int nErrorParam = config_setting_length(errorParam);
		// List of variables and their bounds
		if (nErrorParam == 0) {
			nErrorParam = 2;
			errorVar.reserve(2);
			errorVar.push_back( 0 ); 
			errorVar.push_back( 1 ); 
			errorBounds.reserve(4);
			errorBounds.push_back( allbounds[0] );
			errorBounds.push_back( allbounds[1] );
			errorBounds.push_back( allbounds[2] ); 
			errorBounds.push_back( allbounds[3] );
		}
		else {
			errorVar.reserve(nErrorParam);
			errorBounds.reserve(2*nErrorParam);
			const config_setting_t *oneErrorParam;
			for (int i = 0; i < nErrorParam; ++i) {
				oneErrorParam = config_setting_get_elem(errorParam, i);
				// Read the index of a variable
				errorVar.push_back( config_setting_get_int(config_setting_get_member(oneErrorParam,"var")) );
				// Read the bounds
				lbound = config_setting_get_float_elem( config_setting_get_member(oneErrorParam,"range"), 0 );
				hbound = config_setting_get_float_elem( config_setting_get_member(oneErrorParam,"range"), 1 );
				if ((errorVar[i] == 0) || (errorVar[i] == 1) || (errorVar[i] == 12) || (errorVar[i] == 13) || (errorVar[i] >= 24)) {
					errorBounds.push_back( lbound );
					errorBounds.push_back( hbound );
				}
				else {
					errorBounds.push_back( lbound * PI/180 );
					errorBounds.push_back( hbound * PI/180 );
				}
			}
		}
	}
	
	// Genetic algorithm-------------------------------------------------

	config_setting_t *genetic = config_lookup(cf, "genetic");
	size_t nGenetic = config_setting_length( config_setting_get_member(genetic, "genMax") );
	vector<int> popNum; popNum.reserve(nGenetic);
	vector<long> nAvg; nAvg.reserve(nGenetic);
	for (size_t i = 0; i < nGenetic; ++i) {
		popNum.push_back( config_setting_get_int_elem( config_setting_get_member(genetic, "genMax"), i ) );
		nAvg.push_back(  config_setting_get_int64_elem( config_setting_get_member(genetic, "nAvg"), i ) );
	}
	long nAvgMax = *std::max_element(nAvg.begin(), nAvg.end());
	size_t popSize = config_setting_get_int( config_setting_get_member(genetic, "genSize") );
	double pCros = config_setting_get_float( config_setting_get_member(genetic, "pCros") );
	double pExchange = config_setting_get_float( config_setting_get_member(genetic, "pExch") );
	double pMut = config_setting_get_float( config_setting_get_member(genetic, "pMut") );
	double pUniformVsSinglePoint = config_setting_get_float( config_setting_get_member(genetic, "ratioCros") );
	double pSmallCreepVsRandom = config_setting_get_float( config_setting_get_member(genetic, "ratioMut") );
	bool save_best_chromosome = config_setting_get_bool( config_setting_get_member(genetic, "elitism") );

	config_destroy(cf);
	
	////////////////////////////////////////////////////////////////////////////
	//                         GENETIC ALGORITHM                              //
	////////////////////////////////////////////////////////////////////////////
	
	cout << endl << "*** Optimization via genetic algorithm ***" << endl;
		// Start timer
	#ifdef __gnu_linux__
		struct timeval start, end;
		gettimeofday(&start, NULL);
	#endif
	
	// Generate an enviroment object
	Enviroment* penviroment = new Enviroment(nOffset,paramSwitches);
	// Create spin systems
	penviroment->createSpins(spinA_g,spinA_n,spinA_I,spinA_A,spinA_gStrain,spinA_AStrain,spinA_lwpp,
	                         spinB_g,spinB_n,spinB_I,spinB_A,spinB_gStrain,spinB_AStrain,spinB_lwpp);
	// Create a magnetic field
	penviroment->createField(valueField,nAvgMax);
	// Create a microwave pulses
	penviroment->createPulses(detPiLength,detPiHalfLength,pumpPiLength,detFreq,pumpFreq);
	// Set an experimental dataset
	penviroment->setData(expTime,expSignal,scoreIndex);

	// Generate a primary population
	size_t p(1); // Number of the population
	Population* ppopulation = new Population(popSize,chrLength,allbounds,paramSwitches);
	cout << "Generation no. " << p;
	// Apply selection pressure for the initial population
	penviroment->runScoring(*ppopulation,scoreIndex,nAvg[0]);
	// Range chromosomes in the population
	ppopulation->range();
	cout << "\tScore = " << ppopulation->chromosomes[0].fitness << endl;
	// Record statistics
	penviroment->recordScore(ppopulation->chromosomes[0].fitness,fileroot,p);
	penviroment->recordParameters(ppopulation->chromosomes[0].genes,fileroot);
	penviroment->recordFit(ppopulation->chromosomes[0].genes,fileroot,p);
	
	// Record the EPR spectrum of the spin pair
	if (Spectrum) {
		cout << endl << "*** Recording the EPR spectrum of the spin pair ***" << endl;
		penviroment->recordSpectrum(ppopulation->chromosomes[0].genes,fileroot);
		cout << " Done! " << endl << endl;
	}

	int pm = p;
	for (size_t m = 0; m < nGenetic; ++m) {
		while (pm < popNum[m]) {
			// Generate next population
			++pm;
			++p;
			ppopulation->produceOffspring(pCros,pMut,pUniformVsSinglePoint,pSmallCreepVsRandom,pExchange,save_best_chromosome);
			cout << "Generation no. " << p;
			// Apply selection pressure for the new population
			penviroment->runScoring(*ppopulation,scoreIndex,nAvg[m]);
			// Range chromosomes in the population
			ppopulation->range();
			cout << "\tScore = " << ppopulation->chromosomes[0].fitness << endl;
			// Record statistics
			penviroment->recordScore(ppopulation->chromosomes[0].fitness,fileroot,p);
			penviroment->recordParameters(ppopulation->chromosomes[0].genes,fileroot);
			if ((p % 50) == 0) penviroment->recordFit(ppopulation->chromosomes[0].genes,fileroot,p);
		}
		pm = 0;
	}
	cout << " Done! " << endl;

	////////////////////////////////////////////////////////////////////////////
	//                         RECORD STATISTICS                              //
	////////////////////////////////////////////////////////////////////////////

	// Record symmetry-related solutions
	if (SymSol) {
		cout << endl << "*** Creating the symmetry-related solutions ***" << endl;
		penviroment->recordSymmetricSolutions(ppopulation->chromosomes[0].genes,allbounds,scoreIndex,fileroot);
		cout << " Done! " << endl;
	}

	// Record the modulation depth parameter vs dipolar angle
	if (modDepth) {
		cout << endl << "*** Recording the angular dependence of the modulation depth parameter ***" << endl;
		penviroment->recordModDepth(ppopulation->chromosomes[0].genes,fileroot);
		cout << " Done! " << endl;
	}
	
	// Record the error profile for the distance parameters
	if (errorProfile) {
		cout << endl << "*** Recording the RMSD surface ***" << endl;
		penviroment->recordErrorProfile(ppopulation->chromosomes[0].genes,errorVar,errorBounds,fileroot);
		cout << " Done! " << endl;
	}

	// Delete objects
	delete ppopulation;
	delete penviroment;
	
	// Stop timer
	#ifdef __gnu_linux__
		gettimeofday(&end, NULL);
		cout << "*** Calculation took " << ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / (3600.0 * 1.e6) << " hour(s)! ***" << endl << endl;
	#endif

	return 0;
}