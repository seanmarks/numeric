#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

#include "numeric/Aligned.h"
#include "numeric/Assert.h"
#include "numeric/CommonTypes.h"
#include "numeric/CompareNumbers.h"
#include "numeric/ComplexVector.h"
#include "numeric/Constants.h"
#include "numeric/GptlWrappers.h"
#include "numeric/Matrix.h"
#include "numeric/SphericalHarmonics.h"
#include "numeric/Timer.h"

using namespace CommonTypes;

void matrix_ops()
{
	using Matrix = numeric::Matrix<double>;
	using Vector = CommonTypes::Vector<double>;

	int n = 4000;

	Matrix A(n, n, 1.0);
	Vector x(n, 2.0), b(n, 3.0);

	Vector result(n, 0.0);

	std::vector<Timer> timers;

	int num_iterations = 100;

	// Approach #1: no cache blocking
	timers.push_back( Timer() );
	timers.back().start();
	for ( int k=0; k<num_iterations; ++k ) {
		for ( int i=0; i<n; ++i ) {
			for ( int j=0; j<n; ++j ) {
					result[i] += A(i,j)*x[i];
			}
		}
	}
	timers.back().stop();

	// Approach #2: cache blocking
	timers.push_back( Timer() );
	timers.back().start();
	for ( int k=0; k<num_iterations; ++k ) {
		int block_size = 128;

		for ( int j_outer = 0; j_outer < n; j_outer += std::min(block_size, n-j_outer) ) {
			int max_j = std::min(n, j_outer + block_size);
			for ( int i=0; i<n; ++i ) {
				for ( int j_inner = j_outer; j_inner < max_j; ++j_inner ) {
						result[i] += A(i, j_inner)*x[i];
				}
			}
		}
	}
	timers.back().stop();

	// Print performance comparison
	int num_timers = timers.size();
	for ( int t=0; t<num_timers; ++t ) {
		double dt = timers[t].get_duration_ms();
		std::cout << "Method " << t+1 << ": " << dt << " ms";

		if ( t > 0 ) {
			double dt_ref  = timers[0].get_duration_ms();
			double speedup = 100.0*(dt_ref/dt - 1.0);
			std::cout << " (" << speedup << " speedup)";
		}
		std::cout << std::endl;
	}
}


int main(int argc, char* argv[])
{
	//----- Ylm-like calculation -----//

	int l = 6;
	int num_m_values = 2*l + 1;

	int num_points = 100;

	using Matrix = numeric::Matrix<double>;
	Matrix Y(num_m_values, num_points);
	Y.assign(1.0);

	// Gradient
	std::array<Matrix, N_DIM> dY;
	dY.fill(Y);

	// Weight functions (and derivatives)
	Vector<double> h(num_points, -1.0);
	VectorReal3    dh(num_points, {{1.0, 2.0, 3.0}});

	// Test vectorization
	auto result = dY;
	for ( int d=0; d<N_DIM; ++d ) {
		for ( int m=0; m<num_m_values; ++m ) {
			#pragma omp parallel for simd schedule(static)
			for ( int i=0; i<num_points; ++i ) {
				result[d](m,i) = h[i]*dY[d](m,i) + dh[i][d]*Y(m,i);
			}
		}
	}

	matrix_ops();
}