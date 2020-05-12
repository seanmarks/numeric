#include "main.h"

template<typename T>
using AlignedAlloc = numeric::aligned::CacheAlignedAllocator<T>;

int main(int argc, char* argv[])
{
	//-----------------//
	//----- Setup -----//
	//-----------------//

	using Real = double;

	// numeric
	using namespace numeric;
	using RealVector    = std::vector<Real, AlignedAlloc<Real>>;
	using ComplexVector = aligned::ComplexVector<Real>;

	// STL
	using Complex          = std::complex<Real>;
	using StdRealVector    = std::vector<Real>;
	using StdComplexVector = std::vector<Complex>;

	// Convert 'argv'
	std::vector<std::string> args(argc);
	for ( int c=0; c<argc; ++c ) {
		args[c] = std::string( argv[c] );
	}


	//-----------------------------//
	//----- Interface Testing -----//
	//-----------------------------//

	{
		Real value;
		unsigned size;

		// Construction
		size  = 10;
		value = 3.0;
		ComplexVector z(size, value);
		FANCY_ASSERT( z.size() == size, "incorrect size: " << z.size() );
		for ( unsigned i=0; i<size; ++i ) {
			FANCY_ASSERT( z.real(i) == value, "incorrect value: " << z.real(i) );
			FANCY_ASSERT( z.imag(i) == 0.0,   "incorrect value: " << z.imag(i) );
		}

		// Fill
		value = 5.0;
		z.fill(value);
		FANCY_ASSERT( z.size() == size, "incorrect size: " << z.size() );
		for ( unsigned i=0; i<size; ++i ) {
			FANCY_ASSERT( z.real(i) == value, "incorrect value: " << z.real(i) );
			FANCY_ASSERT( z.imag(i) == 0.0,   "incorrect value: " << z.imag(i) );
		}

		/*
		// TODO
		a_new.fill(5.0);
		std::cout << a_new << std::endl;

		a_new.fill(Complex{{1.0, -1.0}});
		std::cout << a_new << std::endl;

		a_new.resize(15);

		std::cout << a_new;
		std::cout << "size = " << a_new.size() << std::endl;

		a_new.real(2) = 1.5;
		a_new.imag(4) = 7.0;
		std::cout << a_new << std::endl;

		a_new(5) = Complex{{20.0, 30.0}};
		std::cout << a_new << std::endl;

		Complex a_5 = a_new[5];
		std::cout << a_5 << std::endl;

		const ComplexVector& b_new = a_new;
		std::cout << "b_new[5] = " << b_new[5] << std::endl;
		Complex b_5 = b_new[5];
		std::cout << "b_5 = " << b_5 << std::endl;

		std::cout << a_5 + b_new[0] << std::endl;
		*/
	}


	//------------------------//
	//----- SIMD Testing -----//
	//------------------------//

	std::cout << "//-------------------------//" << std::endl;
	std::cout << "//----- ComplexVector -----//" << std::endl;
	std::cout << "//-------------------------//" << std::endl;
	std::cout  << std::endl;

	// Number of values
	int num_neigh = 5;
	int l = 6;
	int len = num_neigh*l;  // default
	len = 5000;

	// Number of iterations to perform
	int num_iterations = 100000;
	//int num_iterations = 200000;
	//int num_iterations = 1000000;

	// Input handling
	if ( argc > 1 ) {
		len = std::stoi( args[1] );
	}
	FANCY_ASSERT(len >= 1, "invalid length");
	if ( argc > 2 ) {
		num_iterations = std::stoi( args[2] );
	}
	FANCY_ASSERT(len >= 1, "invalid length");


	aligned::Allocator<Real> loose_allocator;
	std::cout << "OpenMP SIMD Testing\n"
	          << "      alignment:  " << loose_allocator.get_alignment() << "\n"
	          << "  vector length:  " << len << "\n"
	          << "     iterations:  " << num_iterations << "\n"
	          << std::endl;

	Timer timer_old, timer_new;
	std::string header;

	// Sample values
	Complex alpha = {{ 1.0,  2.0 }};
	Complex beta  = {{ 3.0, -4.0 }};
	Complex one   = {{ 1.0,  0.0 }};
	//Complex unit  = {{ 1.0/sqrt(2.0), -1.0/sqrt(2.0) }};
	Complex unit  = {{ 0.0, -1.0 }};
	StdComplexVector u_old(len, alpha);  ComplexVector u_new      = u_old;
	StdComplexVector v_old(len, beta);   ComplexVector v_new      = v_old;
	StdComplexVector output_old(len);    ComplexVector output_new = output_old;


	std::cout << "//----- Vector <op> Vector -----//" << std::endl;

	header = "Add:  output = u + v";

	timer_old.start();
	for ( int k=0; k<num_iterations; ++k ) {
		for ( int i=0; i<len; ++i ) {
			output_old[i] = u_old[i] + v_old[i];
		}
	}
	timer_old.stop();
	std::cout << "  u_old[0] = " << u_old[0]          << std::endl;

	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		aligned::simd::complex::add(u_new.data(), v_new.data(), len, output_new.data());
	}
	timer_new.stop();
	std::cout << "  u_new[0] = " << Complex(u_new[0]) << std::endl;

	check(output_new, output_old);
	comparePerformance(header, rmsd(output_new, output_old), timer_new, timer_old);


	header = "Mul:  output = u * v  (element-wise)";

	timer_old.start();
	for ( int k=0; k<num_iterations; ++k ) {
		for ( int i=0; i<len; ++i ) {
			output_old[i] = u_old[i]*v_old[i];
		}
	}
	timer_old.stop();
	std::cout << "  output_old[0] = " << output_old[0] << std::endl;

	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		aligned::simd::complex::multiply(u_new.data(), v_new.data(), len, output_new.data());
	}
	timer_new.stop();
	std::cout << "  output_new[0] = " << Complex(output_new[0]) << std::endl;

	check(output_new, output_old);
	comparePerformance(header, rmsd(output_new, output_old), timer_new, timer_old);


	header = "Div:  output = u * v  (element-wise)";

	timer_old.start();
	for ( int k=0; k<num_iterations; ++k ) {
		for ( int i=0; i<len; ++i ) {
			output_old[i] = u_old[i]/v_old[i];
		}
	}
	timer_old.stop();
	std::cout << "  output_old[0] = " << output_old[0] << std::endl;

	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		aligned::simd::complex::divide(u_new.data(), v_new.data(), len, output_new.data());
	}
	timer_new.stop();
	std::cout << "  output_new[0] = " << Complex(output_new[0]) << std::endl;

	check(output_new, output_old);
	comparePerformance(header, rmsd(output_new, output_old), timer_new, timer_old);


	std::cout << "//----- Vector <op> Scalar (in-place) -----//" << std::endl;


	header = "Mul:  output = alpha*output  (element-wise)";

	output_old.assign(len, alpha);
	timer_old.start();
	for ( int k=0; k<num_iterations; ++k ) {
		for ( int i=0; i<len; ++i ) {
			output_old[i] = unit*output_old[i];
		}
	}
	timer_old.stop();
	std::cout << "  output_old[0] = " << output_old[0] << std::endl;

	output_new.assign(len, alpha);
	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		aligned::simd::complex::left_multiply_in_place(unit.real(), unit.imag(), len, output_new.data());
	}
	timer_new.stop();
	std::cout << "  output_new[0] = " << Complex(output_new[0]) << std::endl;

	check(output_new, output_old);
	comparePerformance(header, rmsd(output_new, output_old), timer_new, timer_old);

	return 0;
}
