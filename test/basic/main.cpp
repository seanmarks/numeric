
#include <iostream>
#include <memory>

#include "main.h"

#include "numeric/Array.h"
#include "numeric/OpenMP.h"
#include "numeric/CArrayRef.h"

static constexpr int DIM = 3;

int main(int argc, char* argv[])
{
	std::cout << "Hello" << "\n";

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


	//----- Test ref container -----//

	using rvec = double[DIM];
	int len = 4;

	std::unique_ptr<rvec> coords = std::unique_ptr<rvec>( new rvec[len] );
	rvec* coords_ptr = coords.get();

	CArrayRef<rvec> coords_ref( &coords_ptr, &len );
	std::cout << "coords.size() [ref] = " << coords_ref.size() << "\n";
	int i = 0;
	for ( auto it = coords_ref.begin(); it != coords_ref.end(); ++it, ++i ) {
		std::cout << "i: ";

		for ( int d=0; d<DIM; ++d ) {
			(*it)[d] = 4+3;

			std::cout << " " << coords_ptr[i][d];
		}
		std::cout << "\n";
	}


	const rvec* coords_const_ptr = coords.get();
	CArrayConstRef<rvec> coords_const_ref( &coords_const_ptr, &len );

	std::cout << coords_const_ref[2][0] << std::endl;

	return 0;
}
