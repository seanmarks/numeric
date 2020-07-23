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

int main(int argc, char* argv[])
{
	int l = 6;
	int len = 100;

	using Matrix = numeric::Matrix<double>;
	Matrix Y(l, len);
	Y.assign(1.0);

	std::array<Matrix, N_DIM> dY;
	dY.fill(Y);

	Vector<double> h(len, -1.0);
	VectorReal3    dh(len, {{1.0, 2.0, 3.0}});

	auto result = dY;
	for ( int d=0; d<N_DIM; ++d ) {
		for ( int m=0; m<l; ++m ) {
			#pragma omp parallel for simd schedule(static)
			for ( int i=0; i<len; ++i ) {
				result[d](m,i) = h[i]*dY[d](m,i) + dh[i][d]*Y(m,i);
			}
		}
	}
}
