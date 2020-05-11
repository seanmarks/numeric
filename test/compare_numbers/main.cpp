#include <iostream>
#include <limits>

#include "numeric/Assert.h"
#include "numeric/CompareNumbers.h"

int main(const int argc, const char* argv [])
{
	// Two numbers, one ULP apart
	float a = 67329.234;
	float b = 67329.242;

	numeric::AlmostEqualUlps<float> almost_equal_ulps_f;
	FANCY_ASSERT( almost_equal_ulps_f(a,b), "numbers should be within 1 float-ULP" );
	/*

	using numeric::Float_t;
	Float_t rep_a(a);  
	Float_t rep_b(b);
	std::cout << "a = " << a << ", b = " << b << "\n"
	          << "  int(a) = " << rep_a.integer_rep() << "\n"
	          << "  int(b) = " << rep_b.integer_rep() << "\n"
	          << "    int(b) - int(a) = " << rep_b.integer_rep() - rep_a.integer_rep() << "\n" 
	          << "  almost_equal(ulps) = " << almost_equal_ulps_f(a,b) << "\n";
	*/

	// Two numbers (further apart)
	double c = 0.0;
	double d = 1.0e-8;
	numeric::AlmostEqualUlps<double> almost_equal_ulps_d;
	FANCY_ASSERT( ! almost_equal_ulps_d(c,d), "numbers should be far apart in double-ULPs" );

	/*
	using numeric::Double_t;
	Double_t rep_c(c);
	Double_t rep_d(d);

	std::cout << "c = " << c << ", d = " << d << "\n"
	          << "  int(c) = " << rep_c.integer_rep() << "\n"
	          << "  int(d) = " << rep_d.integer_rep() << "\n"
	          << "    int(d) - int(c) = " << rep_d.integer_rep() - rep_c.integer_rep() << "\n" 
	          << "  almost_equal(ulps) = " << almost_equal_ulps_d(c,d) << "\n";
	          //<< "  almost_equal(ulps) = " << numeric::almost_equal_ulps(c,d) << "\n";
	*/
}
