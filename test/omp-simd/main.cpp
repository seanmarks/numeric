
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

#include "omp-simd/Aligned.h"
#include "omp-simd/CommonTypes.h"
#include "omp-simd/ComplexVector.h"
#include "omp-simd/Timer.h"
#include "omp-simd/VectorComplexN.h"


void comparePerformance(std::string header, const double rmsd, const Timer& timer_old, const Timer& timer_new)
{
	double t_old_ms = timer_old.get_duration_ms();
	double t_new_ms = timer_new.get_duration_ms();

	double t_old_us = timer_old.get_duration_us();
	double t_new_us = timer_new.get_duration_us();
	double speedup = 100.0*(t_old_us/t_new_us - 1.0);

	if ( header.empty() ) {
		header = "Compare";
	}
	else if ( header.back() == '\n' ) {
		header.pop_back();
	}

	std::cout << header << "\n"
	          << "  (rmsd = " << rmsd << ")\n"
	          << "  old: " << t_old_ms << " ms (" << t_old_us << " us)\n"
	          << "  new: " << t_new_ms << " ms (" << t_new_us << " us)\n"
	          << "  speedup: " << speedup << "%\n"
	          << std::endl;
}

// Cache-aligned vector
template<typename T>
using Alloc  = Aligned::CacheAlignedAllocator<T>;
template<typename T>
using Vector = std::vector<T, Alloc<T>>;

// Custom complex vector types
template<typename T> 
using ComplexVector = Aligned::ComplexVector<T,std::vector,Alloc>;
template<typename T, std::size_t N> 
using VectorComplexN = Aligned::VectorComplexN<T,N,std::vector,Alloc>;

// Standard types
template<typename T>
using Complex = std::complex<T>;
template <typename T>
using StdComplexVector = std::vector<Complex<T>>;
//using StdComplexVector = CommonTypes::AlignedVector<std::complex<T>>;


template<typename T>
T rmsd(const StdComplexVector<T>& std_vec, const ComplexVector<T>& new_vec)
{
	T sum = 0.0;
	unsigned length = std_vec.size();
	assert( length == new_vec.size() );

	for ( unsigned i=0; i<length; ++i ) {
		sum += std::norm( std_vec[i] - new_vec(i) );  // std::norm returns |z|^2

		// FIXME DEBUG
		//std::cout << "i=" << i << ": " << std_vec[i] << "   " << new_vec(i) << "\n";
	}
	return sqrt( sum/static_cast<double>(length) );
}


template <typename T, std::size_t N>
using StdVectorComplexN = std::vector< std::array<Complex<T>,N> >;

template<typename T, std::size_t N>
T rmsd(const StdVectorComplexN<T,N>& std_vec, const VectorComplexN<T,N>& new_vec)
{
	T sum = 0.0;
	unsigned length = std_vec.size();
	assert( length == new_vec.size() );

	for ( unsigned i=0; i<length; ++i ) {
		for ( unsigned n=0; n<N; ++n ) {
			sum += std::norm( std_vec[i][n] - new_vec(i,n) );  // std::norm returns |z|^2
		}
	}
	return sqrt( sum/static_cast<double>(length*N) );
}


int main(int argc, char* argv[])
{
	//-----------------//
	//----- Setup -----//
	//-----------------//

	using Real = double;
	using CommonTypes::Complex;

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

	// Number of iterations to perform
	//int num_iterations = 100000;
	int num_iterations = 200000;
	//int num_iterations = 1000000;

	Aligned::Allocator<Complex> loose_allocator;
	std::cout << "OpenMP SIMD Testing\n"
	          << "      alignment:  " << loose_allocator.get_alignment() << "\n"
	          << "  vector length:  " << len << "\n"
	          << "     iterations:  " << num_iterations << "\n"
	          << std::endl;

	Timer timer_std, timer_new;


	//-------------------------//
	//----- ComplexVector -----//
	//-------------------------//

	// Vectors of new type
	ComplexVector<Real> a(len), b(len), c(len), d(len);
	ComplexVector<Real> output(len), zeros(len), ones(len);
	for ( int i=0; i<len; ++i ) {
		a.setValue(i, Complex(     i,            i + 1.0));
		b.setValue(i, Complex(     i - 1.0,      i - 5.0));
		c.setValue(i, Complex(     i + 2.0,  2.0*i + 3.0));
		d.setValue(i, Complex(-3.0*i + 4.0,  7.0*i - 3.0));

		output.setValue(i, Complex(0.0, 0.0));
		zeros.setValue(i, Complex(0.0, 0.0));
		ones.setValue(i, Complex(1.0, 1.0));
	}

	// Create std::vectors with the same values for comparison
	auto a_std = a.convert<StdComplexVector<Real>>();
	auto b_std = b.convert<StdComplexVector<Real>>();
	auto c_std = c.convert<StdComplexVector<Real>>();
	auto d_std = d.convert<StdComplexVector<Real>>();
	auto output_std = output.convert<StdComplexVector<Real>>();
	auto zeros_std  = zeros.convert<StdComplexVector<Real>>();
	auto ones_std   = ones.convert<StdComplexVector<Real>>();

	std::string header;


	//-----  Addition -----//

	header = "ComplexVector addition (in place)";
	output_std = ones_std;
	timer_std.start();
	for ( int k=0; k<num_iterations; ++k ) {
		#pragma omp simd
		for ( int i=0; i<len; ++i ) {
			//output_std[i] = a_std[i] + b_std[i];
			output_std[i] += a_std[i];
			//output_std[i] = a_std[i] + b_std[i] + d_std[i];
		}
	}
	timer_std.stop();

	output = ones;
	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		//output = a;  output += b;
		//output += a;  //  output += b;  output += d;
		//Aligned::add(a, b, output);
		Aligned::add(output, a, output);  // adds in place
	}
	timer_new.stop();
	comparePerformance(header, rmsd(output_std, output), timer_std, timer_new);

	header = "Add complex value to complex vectors";
	Complex z0(1.0, 0.0);
	output_std = ones_std;
	timer_std.start();
	for ( int k=0; k<num_iterations; ++k ) {
		#pragma omp simd
		for ( int i=0; i<len; ++i ) {
			output_std[i] += z0;
		}
	}
	timer_std.stop();
	output = ones;
	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		output += z0;
	}
	timer_new.stop();
	comparePerformance(header, rmsd(output_std, output), timer_std, timer_new);


	//-----  Multiplication -----//

	header = "Multiply complex vectors";
	output_std = ones_std;
	timer_std.start();
	for ( int k=0; k<num_iterations; ++k ) {
		#pragma omp simd
		for ( int i=0; i<len; ++i ) {
			output_std[i] = a_std[i]*b_std[i];
		}
	}
	timer_std.stop();
	output = ones;
	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		Aligned::multiply(a, b, output);
	}
	timer_new.stop();
	comparePerformance(header, rmsd(output_std, output), timer_std, timer_new);

	header = "Multiply by complex value";
	output_std = ones_std;
	timer_std.start();
	for ( int k=0; k<num_iterations; ++k ) {
		#pragma omp simd
		for ( int i=0; i<len; ++i ) {
			output_std[i] *= z0;
		}
	}
	timer_std.stop();
	output = ones;
	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		output *= z0;
	}
	timer_new.stop();
	comparePerformance(header, rmsd(output_std, output), timer_std, timer_new);


	//----- Division -----//

	header = "ComplexVector division";
	output_std = ones_std;
	timer_std.start();
	for ( int k=0; k<num_iterations; ++k ) {
		#pragma omp simd
		for ( int i=0; i<len; ++i ) {
			output_std[i] = a_std[i]/b_std[i];
		}
	}
	timer_std.stop();

	output = ones;
	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		Aligned::divide(a, b, output);  // output = a/b
	}
	timer_new.stop();
	comparePerformance(header, rmsd(output_std, output), timer_std, timer_new);


	//--------------------------//
	//----- VectorComplexN -----//
	//--------------------------//

	std::cout << std::endl;

	using VectorComplex3 = VectorComplexN<Real,3>;
	VectorComplex3 x(len), y(len), z(len);  // 'len' Complex3 values each (equivalent)

	using StdVectorComplex3 = StdVectorComplexN<Real,3>;
	StdVectorComplex3 x_std(len), y_std(len), z_std(len);

	header = "Multiply by complex vector";
	timer_std.start();
	for ( int k=0; k<num_iterations; ++k ) {
		//#pragma omp simd
		for ( int i=0; i<len; ++i ) {
			#pragma omp simd
			for ( int d=0; d<3; ++d ) {
				x_std[i][d] *= a_std[i];
			}
		}
	}
	timer_std.stop();
	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		x.multiplyChunk(a);
	}
	timer_new.stop();
	comparePerformance(header, rmsd(x_std, x), timer_std, timer_new);

	return 0;
}




/*

using namespace CommonTypes;

#pragma omp declare simd
void fma_Complex3(Complex3& result, const Complex3& a, const Complex3& b, const Complex3& c) {
	for ( int d=0; d<DIM_; ++d ) {
		result[d] = a[d] + b[d]*c[d];
	}
}

int test_OpenMP()
{
	Timer timer;
	int num_iterations = 100000;

	int l = 6;
	int num_neighbors = 5;
	//int num_neighbors = 500;

	int num_complex3 = l*num_neighbors;
	int num_complex  = DIM_*num_complex3;
	int num_doubles  = 2*num_complex;

	std::cout << "Each iteration: " << num_complex << " Complex" << std::endl;
	std::cout << "Data" << std::endl;
	std::cout << "  sizeof(double)   = " << sizeof(double)   << std::endl;
	std::cout << "  sizeof(Complex)  = " << sizeof(Complex)  << std::endl;
	std::cout << "  sizeof(Complex3) = " << sizeof(Complex3) << std::endl;
	std::cout << std::endl;

	Complex a_dim( 1.0, -1.0);
	Complex b_dim( 1.0,  1.0);
	Complex c_dim(-2.0,  2.0);
	Complex zero(0.0, 0.0);


	Complex3 zero_complex3, a_i, b_i, c_i;
	a_i.fill(a_dim);
	b_i.fill(b_dim);
	c_i.fill(c_dim);
	zero_complex3.fill(zero);

	AlignedVector<Complex3> a, b, c, result;
	a.assign(num_complex3, a_i);
	b.assign(num_complex3, b_i);
	c.assign(num_complex3, c_i);
	result.assign(num_complex3, zero_complex3);

	// Normal
	timer.start();
	for ( int i=0; i<num_iterations; ++i ) {
		for ( int k=0; k<num_complex3; ++k ) {
			for ( int d=0; d<DIM_; ++d ) {
				result[k][d] = a[k][d] + b[k][d]*c[k][d];
			}
		}
	}
	timer.stop();
	std::cout << "Normal loop: " << timer.get_duration_ms() << " ms"
	          << " (" << timer.get_duration_us() << " us)" << std::endl;

	AlignedVector<Complex3> a2 = a, b2 = b, c2 = c, result2 = result;

	// With OpenMP SIMD
	timer.start();
	for ( int i=0; i<num_iterations; ++i ) {
		#pragma omp simd
		for ( int k=0; k<num_complex3; ++k ) {
			for ( int d=0; d<DIM_; ++d ) {
				result2[k][d] = a2[k][d] + b2[k][d]*c2[k][d];
			}
		}
	}
	timer.stop();
	std::cout << "openmp-simd: " << timer.get_duration_ms() << " ms"
	          << " (" << timer.get_duration_us() << " us)" << std::endl;

	// With OpenMP SIMD collapse(n)
	timer.start();
	//#pragma omp simd collapse(3)
	for ( int i=0; i<num_iterations; ++i ) {
		#pragma omp simd collapse(2)
		for ( int k=0; k<num_complex3; ++k ) {
			for ( int d=0; d<DIM_; ++d ) {
				result2[k][d] = a2[k][d] + b2[k][d]*c2[k][d];
			}
		}
	}
	timer.stop();
	std::cout << "openmp-simd-collapse(3): " << timer.get_duration_ms() << " ms"
	          << " (" << timer.get_duration_us() << " us)" << std::endl;

	// With OpenMP SIMD function
	timer.start();
	for ( int i=0; i<num_iterations; ++i ) {
		#pragma omp simd
		for ( int k=0; k<num_complex3; ++k ) {
			fma_Complex3(result2[k], a2[k], b2[k], c2[k]);
		}
	}
	timer.stop();
	std::cout << "openmp-simd(function): " << timer.get_duration_ms() << " ms"
	          << " (" << timer.get_duration_us() << " us)" << std::endl;


	//----- Linearize to Simple Datatypes -----//

	// With OpenMP SIMD, linearized to Complex = std::complex
	std::vector<Complex> a_linear(num_complex, a_dim), b_linear(num_complex, b_dim),
	                     c_linear(num_complex, c_dim), result_linear(num_complex, a_dim);
	timer.start();
	for ( int i=0; i<num_iterations; ++i ) {
		#pragma omp simd
		for ( int k=0; k<num_complex; ++k ) {
			result_linear[k] = a_linear[k] + b_linear[k]*c_linear[k];
		}
	}
	timer.stop();
	std::cout << "openmp-simd(Complex): " << timer.get_duration_ms() << " ms"
	          << " (" << timer.get_duration_us() << " us)" << std::endl;

	// With OpenMP SIMD, linearized to double
	AlignedVector<double> a_double(num_doubles, 1.0), b_double(num_doubles, 2.0),
	                    c_double(num_doubles, -1.0), result_double(num_doubles, 0.0);
	AlignedVector<double> d_double(num_doubles, 3.0), e_double(num_doubles, 5.0);
	timer.start();
	for ( int i=0; i<num_iterations; ++i ) {
		#pragma omp simd
		for ( int k=0; k<num_doubles; ++k ) {
			// (b + ic)*(d + ie) = (b*d - c*e)
			result_double[k] = a_double[k] + (b_double[k]*d_double[k] - c_double[k]*e_double[k]);
		}
	}
	timer.stop();
	std::cout << "openmp-simd(double): " << timer.get_duration_ms() << " ms"
	          << " (" << timer.get_duration_us() << " us)" << std::endl;

	//aligned(a_double.data(), b_double, c_double, d_double, e_double, result_double:64)

	return 0;
}
*/
