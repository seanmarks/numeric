#include "numeric/Assert.h"
#include "numeric/VectorOfVectors.h"

#include <tuple>


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

	unsigned index = 2;
	unsigned test_size = 13;
	vec.resize(index, test_size);
	vec.checkInternalConsistency();
	FANCY_ASSERT( vec.size(index) == test_size,
	              "wrong size (expected " << test_size << ", got " << vec.size(index) << ")" );
	FANCY_ASSERT( vec.capacity(index) >= test_size,
	              "bad capacity (expected >= " << test_size << ", got " << vec.capacity(index) << ")" );

	std::cout << "dump(vec) = \n";
	vec.dump(std::cout);
	std::cout << "\n";


	//----- push_back() inner -----//

	// Set some values to add and test
	std::vector<TestElement<int>> known_elements;

	/*
	test_tuples.push_back( {{ index-1, 0, 100 }} );
	test_tuples.push_back( {{ index-1, 0, 100 }} );
	test_tuples.push_back( {{ 0, 2, 100 }} );
	*/

	// Force a reallocation
	index = 3;
	for ( unsigned j=0; j<min_inner_size+1; ++j ) {
		vec.push_back(index, j);
	}
	vec.checkInternalConsistency();
	for ( unsigned j=0; j<min_inner_size+1; ++j ) {
		FANCY_ASSERT( vec(index, j) == j,
		              "bad value (expected " << j << ", got " << vec(index, j) << ")" );
	}

	std::cout << "dump(vec) = \n";
	vec.dump(std::cout);
	std::cout << "\n";

	// TODO
	// - Force a reallocation and check for valid state afterwards


	//----- push_back outer -----//

	// TODO





	return 0;
};
