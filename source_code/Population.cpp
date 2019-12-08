#include "Population.h"
#include <cmath>
#include <algorithm>
using std::sort;


Population::Population(const size_t& popSize, const int& chrLength, const vector<double>& genesRanges, vector<bool> const& paramswitches):
	size(popSize), length(chrLength), ranges(genesRanges)
{
    chromosomes.reserve(size);
    childs.reserve(size);
	// This vector specifies which genes correspond to which geometric parameters
	indices.reserve(paramswitches.size());
	int count(0);
	for (size_t i = 0; i < paramswitches.size(); ++i) {
		if (paramswitches[i]) { indices.push_back(count); ++count; }
		else                  { indices.push_back(0);              }
	}

    generateInitialChromosomes();
	//acceptDistributions();
}

// Generate an initial population
void Population::generateInitialChromosomes() 
{
	RandGens* prandgen = new RandGens;
	for (size_t i = 0; i < size; ++i) {
		Chromosome* pchromosome = new Chromosome(length);
		pchromosome->generateInitialChromosome(ranges, *prandgen);
		chromosomes.push_back(*pchromosome);
		delete pchromosome;
	}
	delete prandgen;
}

// Generate an initial population
void Population::produceOffspring(const double& pCros, 
								  const double& pMut, 
								  const double& pUniformVsSinglePoint, 
								  const double& pSmallCreepVsRandom, 
								  const double& pExchange, 
								  bool save_best_chromosome) 
{	
	// Initialize random generator
	RandGens* prandgen = new RandGens;

	// Safeguard the strongest chromosomes (Elitism)
	size_t numToSafeguard(0);
	if (save_best_chromosome) {
		if ( size % 2 == 0 ) {
			childs.push_back(chromosomes[0]);
			childs.push_back(chromosomes[1]);
			numToSafeguard = 2;
		}
		else {
			childs.push_back(chromosomes[0]);
			numToSafeguard = 1;
		}
	}

	// Select parents and produce new chromosomes
	for (size_t i = numToSafeguard; i < size; i += 2) 
	{
		// Tournament selection
		Chromosome* pChild1 = new Chromosome(length);
		Chromosome* pChild2 = new Chromosome(length);
		*pChild1 = chromosomes[tournamentSelection(prandgen)];
		*pChild2 = chromosomes[tournamentSelection(prandgen)];
		// cout << "Here" << setw(10) << pChild1->fitness << setw(10) << pChild2->fitness << endl;
			
		// Crossover
		if (prandgen->uniformRand() <= pCros) {
			// Choose the crossover mode
			if (prandgen->uniformRand() >= pUniformVsSinglePoint) { singlePointCrossover(pChild1, pChild2, prandgen); }
			else { uniformCrossover(pChild1, pChild2, pExchange, prandgen); }
		}

		// Mutation
		// Choose the mutation mode
		if (prandgen->uniformRand() >= pSmallCreepVsRandom) { randomMutation(pChild1, pMut, prandgen); }
		else { smallCreepMutation(pChild1, pMut, prandgen); }
		if (prandgen->uniformRand() >= pSmallCreepVsRandom) { randomMutation(pChild2, pMut, prandgen); }
		else { smallCreepMutation(pChild2, pMut, prandgen); }
			
		// Save new chromosomes to the array
		childs.push_back(*pChild1);
		childs.push_back(*pChild2);
		delete pChild1;
		delete pChild2;
	}
	chromosomes = childs;
	acceptDistributions();

	childs.clear();
	delete prandgen;
}

// Range the chromosomes in the population according to their fitness
void Population::range()
{
	sort(chromosomes.begin(), chromosomes.end());
}

// Check in the angular distributions
void Population::acceptDistributions() 
{
	int widthsIndices[10] = {3,5,7,9,11,15,17,19,21,23}; // Indices of the widths
	int widthIndex(0), meanIndex(0);
	double paramMin(0), paramMax(0);

	for (size_t i = 0; i < size; ++i) { // Run over all chromosomes
		for (int j = 0; j < 10; ++j) { // Run over all widths
			widthIndex = indices[widthsIndices[j]]; // Current index of the width parameter in the parameters vector
			meanIndex = widthIndex - 1;
			if (widthIndex) {
				paramMin = chromosomes[i].genes[meanIndex] - 0.5 * chromosomes[i].genes[widthIndex];
				paramMax = chromosomes[i].genes[meanIndex] + 0.5 * chromosomes[i].genes[widthIndex];
				if (paramMin < ranges[2*meanIndex]) {
					chromosomes[i].genes[meanIndex] = 0.5 * (paramMax + ranges[2*meanIndex]);
					chromosomes[i].genes[widthIndex] = paramMax - ranges[2*meanIndex];
				}
				if (paramMax > ranges[2*meanIndex+1]) {
					chromosomes[i].genes[meanIndex] = 0.5 * (paramMin + ranges[2*meanIndex+1]);
					chromosomes[i].genes[widthIndex] = ranges[2*meanIndex+1] - paramMin;
				}
			}
		}
	}
}

// Tournament selection
int Population::tournamentSelection(RandGens* prandgen) 
{
	int numOfCandidate1 = prandgen->uniformRand(size)-1;
	int numOfCandidate2 = prandgen->uniformRand(size)-1;
	if (chromosomes[numOfCandidate1].fitness <= chromosomes[numOfCandidate2].fitness) { return numOfCandidate1; }
	else { return numOfCandidate2; }
}

// Single point crossover
void Population::singlePointCrossover(Chromosome* pChild1, 
									  Chromosome* pChild2, 
									  RandGens* prandgen) 
{
	int point = prandgen->uniformRand(length)-1;
	vector<double> temp = pChild1->genes;
	for (int i = 0; i < length; ++i) { 
		if (i <= point) { pChild1->genes[i] = pChild2->genes[i]; }
		else { pChild2->genes[i] = temp[i]; }
	}
		
	return;
}

// Exchange crossover
void Population::uniformCrossover(Chromosome* pChild1, 
								  Chromosome* pChild2, 
								  const double& pExchange, 
								  RandGens* prandgen) 
{
	for (int i = 0; i < length; ++i) {
		if (prandgen->uniformRand() >= pExchange) { 
			double temp = pChild1->genes[i];
			pChild1->genes[i] = pChild2->genes[i]; 
			pChild2->genes[i] = temp; 
		}
	}
	return;
}

// Random mutation
void Population::randomMutation(Chromosome* pChild, 
								const double& pMut, 
								RandGens* prandgen) 
{
	for (int i = 0; i < length; ++i) {
		if (prandgen->uniformRand() <= pMut) { pChild->genes[i] = pChild->generateRandomGene(ranges[2*i], ranges[2*i+1], *prandgen); }
	}
}

// Small creep mutation
void Population::smallCreepMutation(Chromosome* pChild, 
									const double& pMut, 
									RandGens* prandgen) 
{
	double increment(0);
	double temp(0);
	for (int i = 0; i < length; ++i) {
		if (prandgen->uniformRand() <= pMut) {
			increment = 0.1 * (ranges[2*i+1] - ranges[2*i]); // 10% of the whole range
			for (;;) {
				temp = pChild->genes[i] + prandgen->uniformRand(-increment, increment);
				if ( (temp > ranges[2*i]) && (temp < ranges[2*i+1]) ) { break; }
			}
			pChild->genes[i] = temp;
		}
	}
}
