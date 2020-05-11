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
#include <type_traits>

namespace numeric {


// Type traits struct with a signed Int type of the same size (in bytes)
// as floating-point type T
template<typename T>
struct FloatingPointNumberTraits;


// Stores different representations a given floating-point number
// - TODO: Store a dissection of the number (e.g. sign, mantissa, exponent)
template<typename T>
class FloatingPointNumber
{
 public:
	using Traits = FloatingPointNumberTraits<T>;
	using Int    = typename Traits::Int;

	static_assert(std::is_floating_point<T>::value, "floating-point type required");
	static_assert(std::is_integral<Int>::value,     "integral type required");
	// TODO require (sizeof(T) == sizeof(Int)) and Int be signed

	FloatingPointNumber(const T value = 0.0):
		value_(value),
		integer_rep_( *reinterpret_cast<Int*>(&value_) )
	{}

	bool is_negative() const {
		// TODO: sign bit check?
		return value_ < 0.0;
	}  
	Int integer_rep() const {
		return integer_rep_;
	}

 private:
	T   value_;
	Int integer_rep_;
};


template<>
struct FloatingPointNumberTraits<float>
{
	using Int = std::int32_t;
};

template<>
struct FloatingPointNumberTraits<double>
{
	using Int = std::int64_t;
};

using Float_t  = FloatingPointNumber<float>;
using Double_t = FloatingPointNumber<double>;


template<typename T>
class AlmostEqualUlps
{
 public:
	using Traits = FloatingPointNumberTraits<T>;
	using Int    = typename Traits::Int;

	static_assert(std::is_floating_point<T>::value, "floating-point type required");
	static_assert(std::is_integral<Int>::value,     "integral type required");

	// TODO calibrate default max_diff
	AlmostEqualUlps(
		const T   max_diff      = /*std::sqrt(*/ std::numeric_limits<T>::epsilon() /*)*/,
		const Int max_ulps_diff = 4
	):
		max_diff_(max_diff), max_ulps_diff_(max_ulps_diff)
	{}

	bool operator()(const T& a, const T& b) {
		return are_almost_equal(a, b);
	}

	bool are_almost_equal(const T& a, const T& b) {
		// Check if the numbers are really close -- needed when comparing numbers near zero.
		T abs_diff = std::fabs(a - b);
		if (abs_diff <= max_diff_) {
			return true;
		}

		FloatingPointNumber<T> rep_a(a), rep_b(b);

		// Different signs means they do not match.
		if (rep_a.is_negative() != rep_b.is_negative()) {
			return false;
		}

		// Find the difference in ULPs.
		Int ulps_diff = std::abs(rep_a.integer_rep() - rep_b.integer_rep());
		if (ulps_diff <= max_ulps_diff_) {
			return true;
		}
		else {
			return false;
		}
	}

 private:
	T   max_diff_;
	Int max_ulps_diff_;
};


// TODO: calibrate default max_diff
template<typename T>
bool almost_equal(
	T a, T b, T max_diff = 1.0e-5, T max_rel_diff = std::numeric_limits<T>::epsilon()
)
{
	// Check if the numbers are really close -- needed
	// when comparing numbers near zero.
	T diff = std::fabs(a - b);
	if (diff <= max_diff) {
		return true;
	}

	a = std::fabs(a);
	b = std::fabs(b);
	T largest = (b > a) ? b : a;

	if (diff <= largest * max_rel_diff) {
		return true;
	}
	else {
		return false;
	}
}


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

} // end namespace numeric
