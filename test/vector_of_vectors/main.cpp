// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

//#include "main.h"

#include <array>
#include <random>
#include <tuple>
#include <vector>

#include "numeric/Assert.h"
#include "numeric/Random.h"
#include "numeric/RandomSampler.h"
#include "numeric/VectorOfVectors.h"

template<typename T>
using TestElement = std::tuple<unsigned,unsigned,T>;  // (outer, inner, value)

int main(int argc, char* argv[])
{
	unsigned outer_size = 10;
	unsigned min_inner_size = 5;

	VectorOfVectors<int> vec(outer_size, min_inner_size);
	vec.checkInternalConsistency();

	//std::cout << "vec = \n" << vec << "\n";

	//----- Basic allocation -----//

	std::cout << "Testing basic allocation..." << std::endl;

	// Outer vector
	FANCY_ASSERT( vec.size() == outer_size,
	             "wrong size (expected " << outer_size << ", got " << vec.size() << ")" );
	unsigned min_total_capacity = outer_size*min_inner_size;
	FANCY_ASSERT( vec.capacity() >= min_total_capacity,
	             "bad total capacity (expected >=" << min_total_capacity << ", got " << vec.capacity() << ")" );

	// Inner vectors
	for ( unsigned i=0; i<outer_size; ++i ) {
		FANCY_ASSERT( vec.size(i) == 0,
		              "wrong size (expected " << 0 << ", got " << vec.size(i) << ")" );
		FANCY_ASSERT( vec.capacity(i) >= min_inner_size,
		              "bad subvector capacity (expected >= " << min_inner_size << ", got " << vec.capacity(i) << ")" );
	}
	std::cout << "  Done" << std::endl;


	//----- resize a subvector -----//

	// Grow
	unsigned index = 2;
	unsigned test_size = outer_size + 3;
	vec.resize(index, test_size);
	vec.checkInternalConsistency();
	FANCY_ASSERT( vec.size(index) == test_size,
	              "wrong size (expected " << test_size << ", got " << vec.size(index) << ")" );
	FANCY_ASSERT( vec.capacity(index) >= test_size,
	              "bad capacity (expected >= " << test_size << ", got " << vec.capacity(index) << ")" );


	// Shrink
	index = 2;
	test_size = 5;
	vec.resize(index, test_size);
	vec.checkInternalConsistency();
	FANCY_ASSERT( vec.size(index) == test_size,
	              "wrong size (expected " << test_size << ", got " << vec.size(index) << ")" );
	FANCY_ASSERT( vec.capacity(index) >= test_size,
	              "bad capacity (expected >= " << test_size << ", got " << vec.capacity(index) << ")" );

	// Stay the same
	vec.resize(index, test_size);
	vec.checkInternalConsistency();
	FANCY_ASSERT( vec.size(index) == test_size,
	              "wrong size (expected " << test_size << ", got " << vec.size(index) << ")" );
	FANCY_ASSERT( vec.capacity(index) >= test_size,
	              "bad capacity (expected >= " << test_size << ", got " << vec.capacity(index) << ")" );


	//----- Test iteration -----//

	// Setting and getting values with operator()
	index = 1;
	test_size = 4;
	vec.resize(index, 4);
	for ( unsigned j=0; j<test_size; ++j ) {
		int val = 2*j;
		vec(index, j) = val;
		FANCY_ASSERT( vec(index, j) == val, "bad value (expected " << val << ", got " << vec(index,j) << ")" );
	}
	vec.checkInternalConsistency();

	// Check the number of elements (according to an iterator)
	FANCY_ASSERT( (vec.end(index) - vec.begin(index)) == test_size, "bad size" );

	// Iterators
	int j = 0;
	for ( auto it = vec.begin(index); it != vec.end(index); ++it, ++j ) {
		int val = 2*j;	
		FANCY_ASSERT( *it == 2*j,"bad value (expected " << val << ", got " << *it << ")" );
	}
	vec.checkInternalConsistency();


	//----- push_back() inner -----//

	// Set some values *after* where insertion will happen
	std::vector<TestElement<int>> known_elements;
	unsigned later_index = 7;
	vec.reset(later_index);
	test_size = 3;
	for ( unsigned j=0; j<test_size; ++j ) {
		int val = 4*j - 3;
		vec.push_back(later_index, val);
		known_elements.push_back( std::make_tuple(later_index, j, val) );
	}
	vec.checkInternalConsistency();	

	std::cout << "dump(vec) = \n";
	vec.dump(std::cout);
	std::cout << "\n";

	// Force the reallocation of a subvector
	index = later_index - 3;
	unsigned num = vec.capacity(index) + 3;
	for ( unsigned j=0; j<num; ++j ) {
		vec.push_back(index, j);
	}
	vec.checkInternalConsistency();
	for ( unsigned j=0; j<num; ++j ) {
		FANCY_ASSERT( vec(index, j) == static_cast<int>(j),
		              "bad value (expected " << j << ", got " << vec(index, j) << ")" );
	}

	// Check for the values added earlier
	for ( unsigned k=0; k<known_elements.size(); ++k ) {
		int i   = std::get<0>(known_elements[k]);
		int j   = std::get<1>(known_elements[k]);
		int val = std::get<2>(known_elements[k]);
		FANCY_ASSERT( vec(i, j) == val,
		              "bad value (expected " << val << ", got " << vec(i, j) << ")" );
	}

	std::cout << "dump(vec) = \n";
	vec.dump(std::cout);
	std::cout << "\n";


	//----- push_back outer -----//

	// TODO



	//----- Reset -----//

	vec.reset();
	for ( unsigned i=0; i<vec.size(); ++i ) {
		FANCY_ASSERT( vec.size(i) == 0,
		              "bad value (expected " << 0 << ", got " << vec.size(i) << ")" );
	}


	//----- Neighbor search -----//

	std::cout << "Test using neighbor search\n";

	static constexpr int N_DIM = 3;
	using Real3     = std::array<double,N_DIM>;
	using Positions = std::vector<Real3>;

	// Sample positions
	int num_particles = 1000;
	double box_length = 1.0;
	double r_cut = 0.25;  // cutoff

	// Generate random numbers within the confines of a cubic box
	using Distribution = std::uniform_real_distribution<double>;
	Distribution distribution(0.0, box_length);
	RandomSampler<Distribution> sampler( distribution, Random::getDebugSequence() );
	std::vector<double> sample;
	sampler.generate(N_DIM*num_particles, sample);

	// Set up random positions
	Positions positions(num_particles);
	for ( int i=0; i<num_particles; ++i ) {
		for ( int d=0; d<N_DIM; ++d ) {
			positions[i][d] = sample[i*N_DIM + d];
		}
	}

	// Build a neighbor list using different data structures
	VectorOfVectors<int> neighbor_list(num_particles, 2);
	std::vector<std::vector<int>> ref_list(num_particles);
	for ( int i=0; i<num_particles; ++i ) {
		for ( int j=i+1; j<num_particles; ++j ) {
			// Distance between particles
			double r = 0.0;
			for ( int d=0; d<N_DIM; ++d ) {
				double dx = positions[j][d] - positions[i][d];
				r += dx*dx;
			}
			r = sqrt(r);

			if ( r <= r_cut ) {
				// Save in each structure
				neighbor_list.push_back(i, j);
				neighbor_list.push_back(j, i);

				ref_list[i].push_back(j);
				ref_list[j].push_back(i);
			}
		}
	}

	// Compare results
	for ( int i=0; i<num_particles; ++i ) {
		FANCY_ASSERT( neighbor_list.size(i) == ref_list[i].size(), "size mismatch" );
		if ( neighbor_list.size(i) == 0 ) {
			std::cerr << "particle " << i << " has no neighbors: this is HIGHLY unlikely!\n";
		}

		int num = neighbor_list.size(i);
		for ( int j=0; j<num; ++j ) {
			FANCY_ASSERT( neighbor_list(i,j) == ref_list[i][j], "element mismatch" );
		}
	}


	std::cout << "Done\n";

	return 0;
};
