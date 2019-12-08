#ifndef POPULATION_H
#define POPULATION_H

#include "Chromosome.h"
#include "RandGens.h"
#include "tbb/tbb.h"
using tbb::concurrent_vector;
using tbb::parallel_for;


class Population 
{
private:
	int length;
	vector<int> indices;

public:        
	// Initialization
	Population(const size_t& popSize, const int& chrLength, const vector<double>& genesRanges, vector<bool> const& paramswitches);

	// Generate an initial population
	void generateInitialChromosomes();

	// Generate an initial population
	void produceOffspring(const double& pCros, const double& pMut, const double& pUniformVsSinglePoint, 
	                      const double& pSmallCreepVsRandom, const double& pExchange, bool save_best_chromosome);

	// Range the chromosomes in the population according to their fitness
	void range();

	// Check in the angular distributions
	void acceptDistributions();
	
	// Tournament selection
	int tournamentSelection(RandGens* prandgen);

	// Single point crossover
	void singlePointCrossover(Chromosome* pChild1, Chromosome* pChild2, RandGens* prandgen);

	// Exchange crossover
	void uniformCrossover(Chromosome* pChild1, Chromosome* pChild2, const double& pExchange, RandGens* prandgen);

	// Random mutation
	void randomMutation(Chromosome* pChild, const double& pMut, RandGens* prandgen);

	// Small creep mutation
	void smallCreepMutation(Chromosome* pChild, const double& pMut, RandGens* prandgen);
    
    size_t size;
	//vector<Chromosome> chromosomes;
	//vector<Chromosome> childs;
	concurrent_vector<Chromosome> chromosomes;
	concurrent_vector<Chromosome> childs;
    vector<double> ranges;
};
 
#endif