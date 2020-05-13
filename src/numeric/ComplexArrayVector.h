/*
 * SIMD-friendly vector of complex-valued arrays
 * - Emulates an N-dimensional complex vector:
 *     std::vector< std::array<std::complex<T>,N> >
 *
*/

#pragma once
#ifndef COMPLEX_ARRAY_VECTOR_H
#define COMPLEX_ARRAY_VECTOR_H

#include <array>
#include <cassert>
#include <cstddef>
#include <complex>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Aligned.h"
#include "ComplexVector.h"

namespace numeric {
namespace aligned {

template<typename T>
using DefaultAlignedVector = std::vector<T, CacheAlignedAllocator<T>>;

template<typename T>
using DefaultAlignedComplexVector = ComplexVector<T, CacheAlignedAllocator<T>>;

template<
	typename    T,
	std::size_t N,
	typename    AlignedComplexVector = DefaultAlignedComplexVector<T>
>
class ComplexArrayVector
{
 public:
	// Sanity checks
	//static_assert(std::is_same<T, typename DefaultAlignedComplexVector::value_type>::value, "type mismatch");
	static_assert(N > 0, "empty arrays are not allowed");

	using size_type = std::size_t;

	/*
	// TODO: 
	// - Need to make an ElementView class and/or similar to serve as value/ref/ptr type
	//   since its underlying storate is not std::complex<T>
	using value_type = T;
	using allocator_type = AlignedVector::allocator_type;

	using reference       = value_type&;
	using const_reference = const value_type&;
	using pointer         = value_type*;
	using const_pointer   = const value_type*;

	using iterator               = pointer;
	using const_iterator         = const_pointer;
	using reverse_iterator       = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;
	*/

	using Complex = std::complex<T>;

	ComplexArrayVector() {
		for ( unsigned n=0; n<N; ++n ) {
			data_[n].clear();
		}
	}

	ComplexArrayVector(const size_type size) {
		resize(size);
	}
	ComplexArrayVector(const size_type size, const Complex& value) {
		assign(size, value);
	}
	ComplexArrayVector(const size_type size, const T& value) {
		assign(size, value);
	}

	/*
	template<typename OtherComplexArrayVector,
	         typename std::enable_if< std::is_convertible<typename OtherComplexArrayVector::value_type, Complex>::value >::type* = nullptr
	>
	ComplexArrayVector(const OtherComplexArrayVector& vec) {
		unsigned len = vec.size();
		resize(len);
		for ( unsigned i=0; i<len; ++i ) {
			this->real(i) = vec[i].real();
			this->imag(i) = vec[i].imag();
		}
	}
	*/

	// TODO:
	// - initializer list construction and filling
	// - construction from / conversion to std::vector<std::complex<T>>

	// Manage size
	size_type size() const {
		return data_[0].size();
	}
	void resize(const size_type size) {
		for ( unsigned n=0; n<N; ++n ) {
			data_[n].resize(size);
		}
	}
	void reset() {
		resize(0);
	}
	void resetSize() {
		reset();
	}
	void clear() {
		for ( unsigned n=0; n<N; ++n ) {
			data_[n].clear();
		}
	}
	size_type max_size() const noexcept {
		return std::numeric_limits<size_type>::max()/2;
	}

	// Capacity
	size_type capacity() const {
		size_type min_capacity = data_[0].capacity();
		for ( unsigned n=1; n<N; ++n ) {
			min_capacity = std::min(min_capacity, data_[n].capacity());
		}
		return min_capacity;
	}
	void reserve(const size_type size) {
		for ( unsigned n=0; n<N; ++n ) {
			data_[n].reserve(size);
		}
	}

	// Set size and contents
	// - TODO: set using iterators
	void assign(const size_type size, const Complex& value) {
		for ( unsigned n=0; n<N; ++n ) {
			data_[n].assign(size, value);
		}
	}
	void assign(const size_type size, const T& value) {
		assign(size, {{value, 0.0}});
	}

	// Set value
	template<typename U>
	void fill(const U& value) {
		for ( unsigned n=0; n<N; ++n ) {
			data_[n].fill(value);
		}
	}
	/*
	// TODO overload for particular types/initializer lists?
	void fill(const Complex& value) {
		assign(size_, value);
	}
	*/

	// Special values
	// - All zeros
	void zeros() {
		fill(0.0);
	}
	void zero() {
		zero();
	}
	// - Note that ones() means (1,0) and not (1,1)
	void ones() {
		fill(1.0);
	}
	void one() {
		one();
	}


	//------------------------------------//
	//----- Element/Component Access -----//
	//------------------------------------//

	// Access real/imag components of elements
	T&       real(const size_type i, const size_type n)       { return data_[n].real(i); }
	T&       imag(const size_type i, const size_type n)       { return data_[n].imag(i); }
	const T& real(const size_type i, const size_type n) const { return data_[n].real(i); }
	const T& imag(const size_type i, const size_type n) const { return data_[n].imag(i); }

	// Object that represents a reference to the 'i'th element of
	// a ComplexArrayVector, in the sense of 'vec[i]'
	// - TODO: Rebinding/copying/moving/swap?
	//   - Goal: make ComplexArrayVector compatible with STL algorithms
	// - TODO: ElementConstRef?
	class ElementRef {
	 public:
		ElementRef(ComplexArrayVector& vec, const size_type i):
			vec_(vec), i_(i) 
		{}

		/*
		// TODO

		// Change the parent ComplexArrayVector
		ElementRef& operator=(const Complex& value) {
			real_.get() = value.real();
			imag_.get() = value.imag();
			return *this;
		}

		// Convert to Complex
		operator Complex() const {
			return Complex(real_, imag_);
		}

		// Access components
		T&       real()       { return real_; }
		T&       imag()       { return imag_; }
		const T& real() const { return real_; }
		const T& imag() const { return imag_; }
		*/

	 private:
		std::reference_wrapper<ComplexArrayVector> vec_;
		size_type i_;
	};

	// Element access
	// - By reference (mutable): allows assignment via 'vec[i] = ...'
	/*
	ElementRef operator()(const size_type i) {
		return ElementRef(*this, i);
	}
	ElementRef operator[](const size_type i) {
		return (*this)(i);
	}
	// - By value
	Complex operator()(const size_type i) const {
		return {{ real(i), imag(i) }};
	}
	Complex operator[](const size_type i) const {
		return (*this)(i);
	}
	*/


	//-----------------------//
	//----- Data Access -----//
	//-----------------------//

	/*
	TODO
	// Access underlying block of memory by ptr
	T*       data()       { return data_.data(); }
	const T* data() const { return data_.data(); }

	// Access underlying vector (less common)
	AlignedVector&       vector()       { return data_; }
	const AlignedVector& vector() const { return data_; }

	// Output
	template<typename U, typename V>
	friend std::ostream& operator<<(std::ostream& os, const ComplexArrayVector<U,V>& vec);
	*/


	//----------------------//
	//----- Arithmetic -----//
	//----------------------//

	// TODO: Template expressions?
	ComplexArrayVector operator+(const ComplexArrayVector& other) {
		ComplexArrayVector output(other.size());
		//simd::complex::add( this->data(), other.data(), size_, output.data() );
		return output;
	}
	ComplexArrayVector operator-(const ComplexArrayVector& other) {
		ComplexArrayVector output(other.size());
		//simd::complex::subtract( this->data(), other.data(), size_, output.data() );
		return output;
	}
	ComplexArrayVector operator*(const ComplexArrayVector& other) {
		ComplexArrayVector output(other.size());
		//simd::complex::multiply( this->data(), other.data(), size_, output.data() );
		return output;
	}
	ComplexArrayVector operator/(const ComplexArrayVector& other) {
		ComplexArrayVector output(other.size());
		//simd::complex::divide( this->data(), other.data(), size_, output.data() );
		return output;
	}

	ComplexArrayVector& operator+=(const ComplexArrayVector& other) {
#ifndef NDEBUG
		// TODO: General-purpose consistency checking
		// - Check lengths
		// - Variadic template for checking that none of a set of ptrs alias one another?
		FANCY_ASSERT( this->data() != other.data(), "illegal aliasing" );
#endif
		//simd::complex::add_in_place( other.data(), size_, this->data() );
		return *this;
	}
	ComplexArrayVector& operator-=(const ComplexArrayVector& other) {
		//simd::complex::subtract_in_place( other.data(), size_, this->data() );
		return *this;
	}
	ComplexArrayVector& operator*=(const ComplexArrayVector& other) {
		//simd::real::multiply_in_place( other.data(), 2*size_, this->data() );
		return *this;
	}
	ComplexArrayVector& operator/=(const ComplexArrayVector& other) {
		//simd::real::divide_in_place( other.data(), 2*size_, this->data() );
		return *this;
	}

	// TODO: more operators

 private:
	// Complex values are stored as sequential pairs.
	// - Ex. for a ComplexArrayVector with two complex-value entries, z1 and z1,
	//     data_ = {{ Re(z1), Im(z1),
	//                Re(z2), Im(z2) }}
	std::array<AlignedComplexVector,N> data_;
};


template<typename T, std::size_t N, typename AlignedComplexVector = DefaultAlignedComplexVector<T>>
std::ostream& operator<<(std::ostream& os, const ComplexArrayVector<T,N,AlignedComplexVector>& vec)
{
	unsigned len = vec.size();
	os << "[";
	for ( unsigned i=0; i<len; ++i ) {
		if ( i > 0 ) {
			os << "\n ";
		}
		for ( unsigned n=0; n<N; ++n ) {
			os << "(" << vec.real(i,n) << ", " << vec.imag(i,n) << ")";
		}
	}
	os << "]\n";

	return os;
}

} // end namespace aligned
} // end namespace numeric

//#include "ComplexArrayVectorKernels.h"

#endif // ifndef COMPLEX_ARRAY_VECTOR_H
