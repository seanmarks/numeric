// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "main.h"

template<typename T>
using Vector = Legendre::Vector<T>;

template<typename T>
using Matrix = Legendre::Matrix<T>;


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

	std::cout << "//--------------------//" << std::endl;
	std::cout << "//----- Legendre -----//" << std::endl;
	std::cout << "//--------------------//" << std::endl;
	std::cout  << std::endl;

	// Convert 'argv'
	std::vector<std::string> args(argc);
	for ( int c=0; c<argc; ++c ) {
		args[c] = std::string( argv[c] );
	}

	// Number of values
	int num_points = 10000;

	// Number of iterations to perform
	int num_iterations = 5000;

	// Input handling
	if ( argc > 1 ) {
		num_points = std::stoi( args[1] );
	}
	FANCY_ASSERT(num_points >= 1, "invalid length");
	if ( argc > 2 ) {
		num_iterations = std::stoi( args[2] );
	}
	FANCY_ASSERT(num_iterations >= 1, "invalid iteration count");

	std::cout << "Setup\n"
	          << "  num. points: " << num_points << "\n"
	          << "   iterations: " << num_iterations << "\n"
	          << std::endl;

	Timer timer_old, timer_new;
	std::string header;

	// Set of points to test
	Vector<double> points(num_points);
	for ( int i=0; i<num_points; ++i ) {
		points[i] = static_cast<double>( i % 10 );
	}

	// Set up Legendre objects to test
	std::vector<int> l_values;
	std::vector<std::unique_ptr<Legendre>> legendre_ptrs_;
	std::unique_ptr<Legendre> tmp_ptr = nullptr;

	l_values.push_back(3);
	tmp_ptr = std::unique_ptr<Legendre>( new LegendreP3() );
	legendre_ptrs_.push_back( std::move(tmp_ptr) );

	l_values.push_back(4);
	tmp_ptr = std::unique_ptr<Legendre>( new LegendreP4() );
	legendre_ptrs_.push_back( std::move(tmp_ptr) );

	l_values.push_back(6);
	tmp_ptr = std::unique_ptr<Legendre>( new LegendreP6() );
	legendre_ptrs_.push_back( std::move(tmp_ptr) );


	//----- Run Benchmark -----//

	Vector<Vector<double>> p_vecs(num_points);
	Matrix<double> p_matrix, p_matrix_T;

	int num_l_values = l_values.size();
	for ( int n=0; n<num_l_values; ++n ) {
		std::cout << "//----- Testing l = " << l_values[n] << "-----//" << "\n"
		          << std::endl;

		// Old approach: individual points
		timer_old.start();
		for ( int k=0; k<num_iterations; ++k ) {
			for ( int i=0; i<num_points; ++i ) {
				legendre_ptrs_[n]->calculate(points[i], p_vecs[i]);
			}
		}
		timer_old.stop();

		// New approach #1
		timer_new.start();
		for ( int k=0; k<num_iterations; ++k ) {
			legendre_ptrs_[n]->calculate(points, p_matrix);
		}
		timer_new.stop();
		comparePerformance(header, 0.0, timer_new, timer_old);

		// New approach #2
		timer_new.start();
		for ( int k=0; k<num_iterations; ++k ) {
			legendre_ptrs_[n]->calculate_T(points, p_matrix_T);
		}
		timer_new.stop();
		comparePerformance(header, 0.0, timer_new, timer_old);
	}


	//----- Cleanup -----//

	// GPTL: done
	gptl_timer.stop();
	GPTL::print("gptl.log");
	GPTL::finalize();

	return 0;
}
