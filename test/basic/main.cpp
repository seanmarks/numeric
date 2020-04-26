#include "main.h"

#include <iostream>

#include <numeric/OpenMP.h>
#include <numeric/Array.h>


int main(int argc, char* argv[])
{
	std::cout << "Hello" << "\n";

	const int DIM = 3;
	numeric::Array<double, DIM> x(5.0);

	numeric::Array<float, DIM> y(2.0);
	x /= y;

	numeric::Array<double, DIM> y_d = 4.0;
	x = y;


	std::cout << "x=\n";
	for ( int i=0; i<DIM; ++i ) {
		std::cout << x[i] << "\n";
	}

	std::cout << "y=\n";
	for ( auto it = y.cbegin(); it != y.cend(); ++it ) {
		std::cout << *it << "\n";
	}

	std::cout << "y + y_d = " << y + y_d << "\n";

	return 0;
}
