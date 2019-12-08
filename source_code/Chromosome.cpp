#include "Chromosome.h"


Chromosome::Chromosome(int& chrLength): length(chrLength) 
{
    genes.reserve(length);
}
	
// Generate an initial chromosome
void Chromosome::generateInitialChromosome(const vector<double>& genesRanges, RandGens& randgen) 
{
	double gene(0);
	for (size_t i = 0; i < length; ++i) {
		gene = generateRandomGene(genesRanges[2*i], genesRanges[2*i+1], randgen);
		genes.push_back(gene); 
	}
}
	
// Generate a random gene
double Chromosome::generateRandomGene(const double& min, const double& max, RandGens& randgen) const
{
	return randgen.uniformRand(min, max);
}