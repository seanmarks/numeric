// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "main.h"

template<typename T>
using Vector = CommonTypes::Vector<T>;

template<typename T>
using Matrix = SphericalHarmonics::Matrix<T>;

// Bundle together outputs
struct SphericalHarmonicsOutput
{
	SphericalHarmonics::VectorComplex  Y_l;
	SphericalHarmonics::VectorComplex3 derivs_Y_l;
};


int main(int argc, char* argv[])
{
	// GPTL: start
	GPTL::GlobalOptions options;  // sets default options
	GPTL::initialize();           // Initialize GPTL
	GPTL::Timer gptl_timer("main");
	gptl_timer.start();


	//-----------------//
	//----- Setup -----//
	//-----------------//

	// CommonTypes
	static constexpr int X_DIM = CommonTypes::X_DIM;
	static constexpr int Y_DIM = CommonTypes::Y_DIM;
	static constexpr int Z_DIM = CommonTypes::Z_DIM;
	static constexpr int N_DIM = CommonTypes::N_DIM;
	using Real     = CommonTypes::Real;
	using Complex  = CommonTypes::Complex;
	using Real3    = CommonTypes::Real3;
	using Complex3 = CommonTypes::Complex3;

	using VectorReal3    = CommonTypes::VectorReal3;
	using VectorComplex  = CommonTypes::VectorComplex;
	using VectorComplex3 = CommonTypes::VectorComplex3;

	// numeric
	//using ComplexVector = aligned::ComplexVector<Real>;
	// STL
	using StdComplex       = std::complex<Real>;
	using StdVectorComplex = std::vector<StdComplex>;

	// Convert 'argv'
	std::vector<std::string> args(argc);
	for ( int c=0; c<argc; ++c ) {
		args[c] = std::string( argv[c] );
	}


	//-------------------//
	//----- Testing -----//
	//-------------------//

	std::cout << "//------------------------------//" << std::endl;
	std::cout << "//----- SphericalHarmonics -----//" << std::endl;
	std::cout << "//------------------------------//" << std::endl;
	std::cout  << std::endl;

	// Number of values
	int num_samples_per_angle = 100;

	// Number of iterations to perform
	int num_iterations = 100;
	//int num_iterations = 1000;
	//int num_iterations = 200000;
	//int num_iterations = 1000000;

	// Harmonic index
	int l = 6;

	// Input handling
	if ( argc > 1 ) {
		num_samples_per_angle = std::stoi( args[1] );
	}
	FANCY_ASSERT(num_samples_per_angle >= 1, "invalid length");
	if ( argc > 2 ) {
		num_iterations = std::stoi( args[2] );
	}
	FANCY_ASSERT(num_iterations >= 1, "invalid length");

	int num_points = num_samples_per_angle*num_samples_per_angle;

	std::cout << "Setup\n"
	          << "  num. points: " << num_points << "\n"
	          << "   iterations: " << num_iterations << "\n"
	          << std::endl;

	Timer timer_old, timer_new;
	std::string header;


	//----- Generate a Set of Points -----//

	Vector<double> phi(num_samples_per_angle);    // planar angle
	Vector<double> theta(num_samples_per_angle);  // zenith angle

	double delta_phi   = 2.0*constants::pi/static_cast<double>(num_samples_per_angle);
	double delta_theta = constants::pi/static_cast<double>(num_samples_per_angle);
	
	for ( int i=0; i<num_samples_per_angle; ++i ) {
		phi[i]   = i*delta_phi;
		theta[i] = i*delta_theta;
	} 

	// Same magnitude for all points
	double r = 1.0;
	Vector<double> norms(num_points, r);

	VectorReal3 points(num_points);
	for ( int i=0; i<num_samples_per_angle; ++i ) {
		for ( int j=0; j<num_samples_per_angle; ++j ) {
			int n = i*num_samples_per_angle + j;
			points[n][X_DIM] = norms[n]*cos(phi[i])*sin(theta[j]);
			points[n][Y_DIM] = norms[n]*sin(phi[i])*sin(theta[j]);
			points[n][Z_DIM] = norms[n]*cos(theta[j]);
		}
	}


	//----- Run Benchmark -----//

	bool do_fast_harmonics = true;
	bool need_derivatives  = true;
	SphericalHarmonics sph_harmonics(l, do_fast_harmonics);

	// Old approach
	std::vector<SphericalHarmonicsOutput> old_outputs(num_points);
	timer_old.start();
	for ( int k=0; k<num_iterations; ++k ) {
		/*
		if ( k % 10 == 0 ) {
			std::cout << "old iteration: " << k << std::endl;  // FIXME DEBUG
		}
		*/
		for ( int i=0; i<num_points; ++i ) {
			sph_harmonics.calculate_Y_l( points[i], norms[i], need_derivatives,
			                             old_outputs[i].Y_l, old_outputs[i].derivs_Y_l );
		}
	}
	timer_old.stop();
	std::cout << "Done old" << std::endl;  // FIXME DEBUG

	// New approach
	Matrix<Complex>  Y_l(num_points, l+1);
	Matrix<Complex3> derivs_Y_l(num_points, l+1);
	timer_new.start();
	for ( int k=0; k<num_iterations; ++k ) {
		/*
		if ( k % 10 == 0 ) {
			std::cout << "new iteration: " << k << std::endl;  // FIXME DEBUG
		}
		*/
		sph_harmonics.calculate( points, norms, need_derivatives,
		                         Y_l, derivs_Y_l );
	}
	timer_new.stop();
	std::cout << "Done new" << std::endl;  // FIXME DEBUG

	comparePerformance(header, 0.0, timer_new, timer_old);


	// GPTL: done
	gptl_timer.stop();
	GPTL::print("gptl.log");
	GPTL::finalize();

	return 0;
}
