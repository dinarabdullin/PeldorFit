#ifndef RANDGENS_H
#define RANDGENS_H

#include <random> // for std::normal_distribution
#include <cmath> // for acos
#include "MersenneTwister.h"


class RandGens
{
private:

public:
	MersenneTwister prng;

	RandGens(void) 
	{
		prng.seed_manual();
	}

	// Generate a random number in (0,1].
	double uniformRand() 
	{ 
		return (double)(prng()) / (double)(prng.max()); 
	}

	// Generate a random number in a real interval (a,b].
	double uniformRand(double a, double b) 
	{ 
            return ( (b - a) * uniformRand() + a ); 
	}

	// Generate a random integer in the range [1,n]
	long uniformRand(long n) 
	{
            if (n < 0) n = -(n);
            if (n == 0) return 0;
            long guard = (long)(uniformRand() * n) + 1;
            return (guard > n)? n : guard;
	}

	// Generate a normally distributed random numbers 
	double normal(double mean, double stdev)
	{
		std::normal_distribution<double> f(mean, stdev);
		return f(prng);
	}

	double sineweighted(void) 
	{
		return acos( uniformRand() ); 
	}
};


#endif