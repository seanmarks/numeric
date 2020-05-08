/*
 * SIMD-friendly complex vector
 *
 * - Note: operators that return a ComplexVector rely on copy elision
 *   (NRVO) for performance
*/

#pragma once
#ifndef COMPLEX_VECTOR_H
#define COMPLEX_VECTOR_H

#include <cassert>
#include <cstddef>
#include <complex>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Aligned.h"

namespace aligned {

template<typename T>
using DefaultAlignedVector = std::vector<T, CacheAlignedAllocator<T>>;

template<
	typename T,
	typename AlignedVector = DefaultAlignedVector<T>
>
class ComplexVector
{
 public:
	// Sanity checks
	static_assert(std::is_same<T, typename AlignedVector::value_type>::value, "type mismatch");

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

	ComplexVector(): data_(0), size_(0) {}

	ComplexVector(const size_type size) {
		resize(size);
	}
	ComplexVector(const size_type size, const Complex& value) {
		assign(size, value);
	}
	ComplexVector(const size_type size, const T& value) {
		assign(size, value);
	}

	// TODO: initializer list construction and filling

	// Manage size
	void resize(const size_type size) {
		if ( size_ != size ) {
			size_ = size;
			data_.resize(2*size_);
		}
	}
	void resetSize() {
		resize(0);
	}
	void clear() {
		data_.clear();
	}
	size_type max_size() const noexcept {
		return std::numeric_limits<size_type>::max()/2;
	}

	size_type size() const {
		return size_;
	}
	size_type capacity() const {
		return data_.capacity()/2;
	}

	// Capacity
	void reserve(const size_type size) {
		data_.reserve(2*size);
	}

	//
	void assign(const size_type size, const Complex& value) {
		resize(size);
		unsigned len = data_.size();
		for ( unsigned k=0; k<len; k+=2 ) {
			data_[k]   = value.real();
			data_[k+1] = value.imag();
		}
	}
	void assign(const size_type size, const T& value) {
		assign(size, {{value, 0.0}});
	}

	// Set value
	template<typename U>
	void fill(const U& value) {
		assign(size_, value);
	}
	/*
	// TODO overload for particular types/initializer lists?
	void fill(const T& value) {
		assign(size_, value);
	}
	*/

	// Special values
	void zeros() {
		data_.assign( data_.size(), 0.0 );
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

	// Copy value
	Complex operator()(const size_type i) const {
		unsigned k = 2*i;
		return {{ data_[k], data_[k+1] }};
	}

	// Access real/imag components of elements
	T&       real(const size_type i)       { return data_[2*i];   }
	T&       imag(const size_type i)       { return data_[2*i+1]; }
	const T& real(const size_type i) const { return data_[2*i];   }
	const T& imag(const size_type i) const { return data_[2*i+1]; }

	// Output
	template<typename U, typename V>
	friend std::ostream& operator<<(std::ostream& os, const ComplexVector<U,V>& vec);

 private:
	// Complex values are stored as sequential pairs.
	// - Ex. for a ComplexVector with two complex-value entries, z1 and z1,
	//     data_ = {{ Re(z1), Im(z1),
	//                Re(z2), Im(z2) }}
	AlignedVector data_;
	size_type     size_ = 0;  // number of complex elements = data_.size()/2
};


template<typename T, typename AlignedVector = DefaultAlignedVector<T>>
std::ostream& operator<<(std::ostream& os, const ComplexVector<T,AlignedVector>& vec)
{
	unsigned len = vec.data_.size();
	os << "[";
	for ( unsigned k=0; k<len; k+=2 ) {
		if ( k > 0 ) {
			os << "\n ";
		}
		os << "(" << vec.data_[k] << ", " << vec.data_[k+1] << ")";
	}
	os << "]\n";

	return os;
}

} // end namespace aligned

#endif // ifndef COMPLEX_VECTOR_H
