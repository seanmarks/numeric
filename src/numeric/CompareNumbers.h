/* 
 * Adapted from code written by Bruce Dawson
 * (https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition)
 * 
 * Key facts to note:
 * - 
 */

#include <cmath>
#include <cstdint>
#include <limits>

namespace numeric {

/*
// Original
union Float_t
{
	Float_t(float num = 0.0f) : f(num) {}
	// Portable extraction of components.
	bool    is_negative()  const { return (i >> 31) != 0;      }
	int32_t raw_mantissa() const { return i & ((1 << 23) - 1); }
	int32_t raw_exponent() const { return (i >> 23) & 0xFF;    }

	int32_t i;
	float f;
#ifdef _DEBUG
	struct
	{   // Bitfields for exploration. Do not use in production code.
		uint32_t mantissa : 23;
		uint32_t exponent : 8;
		uint32_t sign : 1;
	} parts;
#endif
};
*/


class Float_t
{
 public:
	using Int = std::int32_t;

	Float_t(const float value = 0.0f):
		value_(value),
		integer_rep_( *reinterpret_cast<Int*>(&value_) )
	{}

	bool is_negative() const {
		// TODO: portable sign bit check?
		return value_ < 0.0;
	}  
	Int integer_rep() const {
		return integer_rep_;
	}

 private:
	float value_;
	Int   integer_rep_;
};


// TODO: default max_diff calibration
bool almost_equal_ulps(
	float a, float b, float max_diff = 1.0e-5, int max_ulps_diff = 4
)
{
	// Check if the numbers are really close -- needed
	// when comparing numbers near zero.
	float abs_diff = std::fabs(a - b);
	if (abs_diff <= max_diff) {
		return true;
	}

	Float_t ua(a);
	Float_t ub(b);

	// Different signs means they do not match.
	if (ua.is_negative() != ub.is_negative()) {
		return false;
	}

	// Find the difference in ULPs.
	int ulps_diff = std::abs(ua.integer_rep() - ub.integer_rep());
	if (ulps_diff <= max_ulps_diff) {
		return true;
	}
	else {
		return false;
	}
}


// TODO: default max_diff calibration
bool almost_equal(
	float a, float b, float max_diff = 1.0e-5, float max_rel_diff = std::numeric_limits<float>::epsilon()
)
{
	// Check if the numbers are really close -- needed
	// when comparing numbers near zero.
	float diff = std::fabs(a - b);
	if (diff <= max_diff) {
		return true;
	}

	a = std::fabs(a);
	b = std::fabs(b);
	float largest = (b > a) ? b : a;

	if (diff <= largest * max_rel_diff) {
		return true;
	}
	else {
		return false;
	}
}

}
