// RandomSampler
// - Use to generate a set of random numbers quickly
//
// TODO 
// - Default constructor
//   - Use a random sequence from Random
// - Convert Engine to template parameter

#ifndef BOOTSTRAP_SUBSAMPLER_H
#define BOOTSTRAP_SUBSAMPLER_H

#include <algorithm>
#include <functional>
#include <memory>
#include <random>

#include "Random.h"

template<class Distribution>
class RandomSampler
{
 public:
	using Engine       = std::mt19937;
	using SeedSequence = std::seed_seq;

	using value_type = typename Distribution::result_type;
	//using Distribution = std::uniform_int_distribution<int>;

	RandomSampler(
		const Distribution& distribution,
		const std::vector<int>& seeds
	);

	// Produces a random sample
	// - TODO: return by value?
	void generate(
		const int num_samples,            // number of samples to generate
		std::vector<value_type>& samples  // samples generated
	);

 private:
	// Distribution from which to sample
	Distribution distribution_;

	// Since std::seed_seq is not copyable, need to store in a unique_ptr to enable reseeding.
	// Store seeds separately.
	std::vector<int> seeds_;
	std::unique_ptr<SeedSequence> seed_sequence_ptr_ = nullptr;

	// TODO possible to allocate statically?
	std::unique_ptr<Engine> mt_engine_ptr_ = nullptr;

};


template<class Distribution>
RandomSampler<Distribution>::RandomSampler(
	const Distribution& distribution, const std::vector<int>& seeds
):
	distribution_(distribution),
	seeds_(seeds),
	seed_sequence_ptr_( new SeedSequence(seeds_.begin(), seeds_.end()) ),
	mt_engine_ptr_( new Engine(*seed_sequence_ptr_) )
{}


template<class Distribution>
void RandomSampler<Distribution>::generate(const int num_samples, std::vector<value_type>& samples)
{
	// Quick function lambda for generating random numbers
	// - TODO: Private member function?
	auto number_generator = [this]() { return distribution_(*mt_engine_ptr_); };
	//auto number_generator = [this]() { return (*distribution_ptr_)(*mt_engine_ptr_); };

	samples.resize(num_samples);
	std::generate( samples.begin(), samples.end(), number_generator );
}

#endif // ifndef BOOTSTRAP_SUBSAMPLER_H
