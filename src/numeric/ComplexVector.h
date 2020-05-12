/*
 * SIMD-friendly complex vector
 *
*/

#pragma once
#ifndef COMPLEX_VECTOR_H
#define COMPLEX_VECTOR_H

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

namespace numeric {
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

	// TODO require OtherComplexVector::value_type convertible to std::complex<T>?
	template<typename OtherComplexVector>
	ComplexVector(const OtherComplexVector& vec) {
		unsigned len = vec.size();
		resize(len);
		for ( unsigned i=0; i<len; ++i ) {
			this->real(i) = vec[i].real();
			this->imag(i) = vec[i].imag();
		}
	}

	// TODO:
	// - initializer list construction and filling
	// - construction from / conversion to std::vector<std::complex<T>>

	// Manage size
	size_type size() const {
		return size_;
	}
	void resize(const size_type size) {
		if ( size_ != size ) {
			size_ = size;
			data_.resize(2*size_);
		}
	}
	void reset() {
		resize(0);
	}
	void resetSize() {
		reset();
	}
	void clear() {
		data_.clear();
	}
	size_type max_size() const noexcept {
		return std::numeric_limits<size_type>::max()/2;
	}

	// Capacity
	size_type capacity() const {
		return data_.capacity()/2;
	}
	void reserve(const size_type size) {
		data_.reserve(2*size);
	}

	// Set size and contents
	// - TODO: set using iterators
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
	void fill(const Complex& value) {
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


	//------------------------------------//
	//----- Element/Component Access -----//
	//------------------------------------//

	// Access real/imag components of elements
	T&       real(const size_type i)       { return data_[2*i];   }
	T&       imag(const size_type i)       { return data_[2*i+1]; }
	const T& real(const size_type i) const { return data_[2*i];   }
	const T& imag(const size_type i) const { return data_[2*i+1]; }

	// Object that represents a reference to the 'i'th element of
	// a ComplexVector, in the sense of 'vec[i]'
	// - TODO: Rebinding/copying/moving/swap?
	//   - Goal: make ComplexVector compatible with STL algorithms
	// - TODO: ElementConstRef?
	class ElementRef {
	 public:
		ElementRef(ComplexVector& vec, const size_type i):
			/*vec_(vec),*/ real_(vec.real(i)), imag_(vec.imag(i)) {}

		// Change the parent ComplexVector
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

	 private:
		//std::reference_wrapper<ComplexVector> vec_;
		std::reference_wrapper<T> real_, imag_;
	};

	// Element access
	// - By reference (mutable): allows assignment via 'vec[i] = ...'
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


	//-----------------------//
	//----- Data Access -----//
	//-----------------------//

	// Access underlying block of memory by ptr
	T*       data()       { return data_.data(); }
	const T* data() const { return data_.data(); }

	// Access underlying vector (less common)
	AlignedVector&       vector()       { return data_; }
	const AlignedVector& vector() const { return data_; }

	// Output
	template<typename U, typename V>
	friend std::ostream& operator<<(std::ostream& os, const ComplexVector<U,V>& vec);


	//----------------------//
	//----- Arithmetic -----//
	//----------------------//

	// TODO: Template expressions?
	ComplexVector operator+(const ComplexVector& other) {
		ComplexVector output(size_);
		simd::complex::add( this->data(), other.data(), size_, output.data() );
		return output;
	}
	ComplexVector operator-(const ComplexVector& other) {
		ComplexVector output(size_);
		simd::complex::subtract( this->data(), other.data(), size_, output.data() );
		return output;
	}
	ComplexVector operator*(const ComplexVector& other) {
		ComplexVector output(size_);
		simd::complex::multiply( this->data(), other.data(), size_, output.data() );
		return output;
	}
	ComplexVector operator/(const ComplexVector& other) {
		ComplexVector output(size_);
		simd::complex::divide( this->data(), other.data(), size_, output.data() );
		return output;
	}

	ComplexVector& operator+=(const ComplexVector& other) {
#ifndef NDEBUG
		// TODO: General-purpose consistency checking
		// - Check lengths
		// - Variadic template for checking that none of a set of ptrs alias one another?
		FANCY_ASSERT( this->data() != other.data(), "illegal aliasing" );
#endif
		simd::complex::add_in_place( other.data(), size_, this->data() );
		return *this;
	}
	ComplexVector& operator-=(const ComplexVector& other) {
		simd::complex::subtract_in_place( other.data(), size_, this->data() );
		return *this;
	}
	/*
	ComplexVector& operator*=(const ComplexVector& other) {
		simd::real::multiply_in_place( other.data(), 2*size_, this->data() );
		return *this;
	}
	ComplexVector& operator/=(const ComplexVector& other) {
		simd::real::divide_in_place( other.data(), 2*size_, this->data() );
		return *this;
	}
	*/

	// TODO: more operators

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
} // end namespace numeric

#include "ComplexVectorKernels.h"

#endif // ifndef COMPLEX_VECTOR_H
