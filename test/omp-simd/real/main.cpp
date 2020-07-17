// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "main.h"

int main(int argc, char* argv[])
{
	//-----------------//
	//----- Setup -----//
	//-----------------//

	using Real = double;

	using namespace numeric;
	using AlignedVector = std::vector<Real, AlignedAlloc<Real>>;

	// Convert 'argv'
	std::vector<std::string> args(argc);
	for ( int c=0; c<argc; ++c ) {
		args[c] = std::string( argv[c] );
	}

	// Number of values
	int num_neigh = 5;
	int l = 6;
	int len = num_neigh*l;  // default
	if ( argc > 1 ) {
		len = std::stoi( args[1] );
	}
	len = 5000;

	// Number of iterations to perform
	int num_iterations = 100000;
	//int num_iterations = 200000;
	//int num_iterations = 1000000;

	aligned::Allocator<Real> loose_allocator;
	std::cout << "OpenMP SIMD Testing\n"
	          << "      alignment:  " << loose_allocator.get_alignment() << "\n"
	          << "  vector length:  " << len << "\n"
	          << "     iterations:  " << num_iterations << "\n"
	          << std::endl;

	Timer timer_std, timer_new;
	std::string header;


	//-----  Addition -----//

	std::vector<Real> a_std(len, 4), b_std(len, 3), output_std(len, 1.0);
	AlignedVector     a_new(len, 4), b_new(len, 3), output_new(len, 1.0);

	header = "Addition (+)";
	
	timer_std.start();
	for ( int k=0; k<num_iterations; ++k ) {
		for ( int i=0; i<len; ++i ) {
			output_std[i] = a_std[i] + b_std[i];
		}
	}
	timer_std.stop();

	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		aligned::simd::real::add(a_new.data(), b_new.data(), len, output_new.data());
	}
	timer_new.stop();

	check(output_new, output_std);
	comparePerformance(header, rmsd(output_std, output_new), timer_std, timer_new);

	return 0;
}
