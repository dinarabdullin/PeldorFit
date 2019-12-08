#ifndef MERSENNE_TWISTER_H
#define MERSENNE_TWISTER_H

#include <random>
#include <cstdint> // for data types
#include <algorithm> // for std::generate_n
#include <functional> // std::ref

class MersenneTwister
{
private:
	std::mt19937_64 prng;
public:
	typedef std::mt19937_64::result_type result_type;
	
	MersenneTwister(void) 
	{
	}

	void seed_manual(void)
	{
		std::uint_least64_t seed_data[std::mt19937_64::state_size];
		std::random_device r;
		std::generate_n(seed_data, std::mt19937_64::state_size, std::ref(r));
		std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
		prng.seed(seq);
	}

	static result_type min() { return 0; }
 
    static result_type max() { return std::numeric_limits<uint_fast64_t>::max(); }

	result_type operator()() { return prng(); }
};

#endif