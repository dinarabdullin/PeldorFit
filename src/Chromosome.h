#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "Definitions.h"
#include "RandGens.h"
#include <vector>
using std::vector;


// Chromosome class
class Chromosome 
{
public:
	size_t length;
	vector<double> genes;
	double fitness;
	//atomic<double> fitness;

	// Initialization
	Chromosome(int& chrLength);
	
	// Generate an initial chromosome
	void generateInitialChromosome(const vector<double>& genesRanges, RandGens& randgen);
	
	// Generates random gene
	double generateRandomGene(const double& min, const double& max, RandGens& randgen) const;

	// Overload less operator
	bool operator < (const Chromosome& A) const { return this->fitness < A.fitness; }
};

#endif