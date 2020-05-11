#include <iostream>
#include <limits>

#include "numeric/CompareNumbers.h"

int main(const int argc, const char* argv [])
{
	// Two numbers, one ULP apart
	float a = 67329.234;
	float b = 67329.242;

	using numeric::Float_t;
	Float_t rep_a(a);  
	Float_t rep_b(b);

	numeric::AlmostEqualUlps<double> almost_equal_ulps_f;

	std::cout << "a = " << a << ", b = " << b << "\n"
	          << "  int(a) = " << rep_a.integer_rep() << "\n"
	          << "  int(b) = " << rep_b.integer_rep() << "\n"
	          << "    int(b) - int(a) = " << rep_b.integer_rep() - rep_a.integer_rep() << "\n" 
	          << "  almost_equal(ulps) = " << almost_equal_ulps_f(a,b) << "\n";

	// Two numbers (further apart)
	double c = 0.0;
	double d = 1.0e-8;

	using numeric::Double_t;
	Double_t rep_c(c);
	Double_t rep_d(d);

	numeric::AlmostEqualUlps<double> almost_equal_ulps_d;

	std::cout << "c = " << c << ", d = " << d << "\n"
	          << "  int(c) = " << rep_c.integer_rep() << "\n"
	          << "  int(d) = " << rep_d.integer_rep() << "\n"
	          << "    int(d) - int(c) = " << rep_d.integer_rep() - rep_c.integer_rep() << "\n" 
	          << "  almost_equal(ulps) = " << almost_equal_ulps_d(c,d) << "\n";
	          //<< "  almost_equal(ulps) = " << numeric::almost_equal_ulps(c,d) << "\n";
}
