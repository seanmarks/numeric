#include "main.h"

template<typename T>
using AlignedAlloc = numeric::aligned::CacheAlignedAllocator<T>;

int main(int argc, char* argv[])
{
	//-----------------//
	//----- Setup -----//
	//-----------------//

	using Real = double;
	constexpr int         DIM = 3;
	constexpr std::size_t N   = DIM;

	// numeric
	using namespace numeric;
	using RealN          = Array<Real,N>;
	using RealVector     = std::vector<Real, AlignedAlloc<Real>>;
	using ComplexVector  = aligned::ComplexVector<Real>;
	using ComplexNVector = aligned::ComplexVector<Real>; // FIXME

	// STL
	using Complex           = std::complex<Real>;
	using ComplexN          = std::array<Complex,N>;
	using StdRealVector     = std::vector<Real>;
	using StdComplexVector  = std::vector<Complex>;
	using StdComplexNVector = std::vector<ComplexN>;

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
	}


	//------------------------//
	//----- SIMD Testing -----//
	//------------------------//

	std::cout << "//--------------------------//" << std::endl;
	std::cout << "//----- ComplexNVector -----//" << std::endl;
	std::cout << "//--------------------------//" << std::endl;
	std::cout  << std::endl;

	// Number of values
	int num_neigh = 5;
	int l = 6;
	int len = num_neigh*l;  // default
	len = 5000;

	// Number of iterations to perform
	int num_iterations = 100000;

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

	// Scalars
	Real a = 2.0, b = -3.0;
	Complex alpha = {{ 1.0,  2.0 }};
	Complex beta  = {{ 3.0, -4.0 }};
	Complex one   = {{ 1.0,  0.0 }};
	Complex unit  = {{ 0.0, -1.0 }};

	// Complex vectors
	const int len_TMP = N*len;
	StdComplexVector u_old(len, alpha);  ComplexVector u_new(len_TMP, alpha);
	StdComplexVector v_old(len, beta);   ComplexVector v_new(len_TMP, beta);
	StdComplexVector output_old(len);    ComplexVector output_new(len_TMP);

	// Arrays of scalars
	RealN da(a), db(b);
	ComplexN dalpha; dalpha.fill(alpha);
	ComplexN dbeta;  dbeta.fill(beta);

	StdComplexNVector du_old(len, dalpha);  ComplexNVector du_new(len_TMP, alpha);
	StdComplexNVector dv_old(len, dbeta);   ComplexNVector dv_new(len_TMP, beta);
	StdComplexNVector doutput_old(len);     ComplexNVector doutput_new(len_TMP);

	header = "TEST:  chain rule,  d(alpha*u) = d(alpha)*u + alpha*du  ('alpha' complex)";

	timer_old.start();
	for ( int k=0; k<num_iterations; ++k ) {
		for ( int i=0; i<len; ++i ) {
			for ( int d=0; d<DIM; ++d ) {
				doutput_old[i][d] = dalpha[d]*u_old[i] + alpha*du_old[i][d];
				// 'a'
				//doutput_old[i][d] = da[d]*u_old[i] + a*du_old[i][d];
			}
		}
	}
	timer_old.stop();
	//std::cout << "  u_old[0] = " << u_old[0]          << std::endl;

	// Group by m-value (packed by DIM)
	{
		ComplexNVector du_new(len_TMP, alpha);
		ComplexNVector dv_new(len_TMP, beta);
		ComplexNVector doutput_new(len_TMP);

		timer_new.start();
		for ( int k=0; k<num_iterations; ++k ) {
			const Real* u  = u_new.data();        // ComplexVector
			const Real* du = du_new.data();       // ComplexNVector
			Real* doutput  = doutput_new.data();  // ComplexNVector

			//int num_complex = DIM*len;

			int stride     = 2*DIM;  // number of Reals per ComplexN
			int num_values = stride*len;  // (2*DIM)*len
			#pragma omp simd aligned(u, du, doutput: CACHE_LINE_SIZE)
			for ( int i=0; i<len; ++i ) {
				for ( int d=0; d<DIM; ++d ) {
					int re_u = 2*i;
					int im_u = re_u + 1;
					int re = i*stride + 2*d;
					int im = re + 1;
					doutput[re] = (dalpha[d].real()*u[re_u] - dalpha[d].imag()*u[im_u]) + (alpha.real()*du[re] - alpha.imag()*du[im]);
					doutput[im] = (dalpha[d].real()*u[im_u] + dalpha[d].imag()*u[re_u]) + (alpha.real()*du[im] + alpha.imag()*du[re]);

					// 'a'
					//doutput[re] = da[d]*u[re_u] + a*du[re];
					//doutput[im] = da[d]*u[im_u] + a*du[im];
				}
			}
		}
		timer_new.stop();
		std::cout << "Group by m-value (i-value)" << std::endl;
		comparePerformance(header, -1.0, timer_new, timer_old);
	}


	// Group by dimension (array of N separate ComplexVectors)
	{
		using ComplexVectorN = std::array<ComplexVector, N>;
		ComplexVectorN du_new;       du_new.fill(u_new);
		ComplexVectorN dv_new;       dv_new.fill(v_new);
		ComplexVectorN doutput_new;  doutput_new.fill(output_new);

		timer_new.start();
		for ( int k=0; k<num_iterations; ++k ) {
			const Real* u = u_new.data();  // ComplexVector
			const int num_values = 2*len;

			for ( int d=0; d<DIM; ++d ) {
				const Real* du       = du_new[d].data();       // ComplexVector
				Real*       doutput  = doutput_new[d].data();  // ComplexVector

				const Real alpha_re  = alpha.real();
				const Real alpha_im  = alpha.imag();
				const Real dalpha_re = dalpha[d].real();
				const Real dalpha_im = dalpha[d].imag();

				//int im;
				#pragma omp simd aligned(u, du, doutput: CACHE_LINE_SIZE) //private(im)
				for ( int re=0; re<num_values; re+=2 ) {
					int im = re + 1;
					doutput[re] = (dalpha_re*u[re] - dalpha_im*u[im]) + (alpha_re*du[re] - alpha_im*du[im]);
					doutput[im] = (dalpha_re*u[im] + dalpha_im*u[re]) + (alpha_re*du[im] + alpha_im*du[re]);
					//doutput[re] = (dalpha[d].real()*u[re] - dalpha[d].imag()*u[im]) + (alpha.real()*du[re] - alpha.imag()*du[im]);
					//doutput[im] = (dalpha[d].real()*u[im] + dalpha[d].imag()*u[re]) + (alpha.real()*du[im] + alpha.imag()*du[re]);
				}
				/*
				// 'a'
				const Real* du       = du_new[d].data();       // ComplexVector
				Real*       doutput  = doutput_new[d].data();  // ComplexVector
				#pragma omp simd aligned(du, doutput: CACHE_LINE_SIZE)
				for ( int k=0; k<num_values; ++k ) {
					doutput[k] = da[d]*u[k] + a*du[k];
				}
				*/
			}
		}
		timer_new.stop();

		//check(output_new, output_old);
		std::cout << "Group by DIM" << std::endl;
		comparePerformance(header, -1.0, timer_new, timer_old);
	}



	/*
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
	*/

	return 0;
}
