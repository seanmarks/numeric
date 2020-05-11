#include <iostream>
#include <limits>

#include "numeric/CompareNumbers.h"

int main(const int argc, const char* argv [])
{
	// Two numbers, one ULP apart
	float a = 67329.234;
	float b = 67329.242;
	std::cout << "a = " << a << ", b = " << b << "\n"
	          << "  almost_equal(ulps) = " << numeric::almost_equal_ulps(a,b) << "\n";
}
