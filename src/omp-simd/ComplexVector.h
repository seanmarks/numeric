// SIMD-friendly complex vector
//
// - Note: operators that return a ComplexVector rely on copy elision
//   (NRVO) for performance

#ifndef COMPLEX_VECTOR_H
#define COMPLEX_VECTOR_H

#include <cassert>
#include <cstddef>
#include <complex>
#include <exception>
#include <stdexcept>
#include <vector>

#include "Aligned.h"
#include "AlignedSIMD.h"
#include "CommonTypes.h"

namespace Aligned {

template<
	typename T,
	template<typename, typename> class VectorType = std::vector,
	template<typename> class AlignedAllocType = CacheAlignedAllocator
>
class ComplexVector
{
 public:
	// TODO: Remove?
	template<typename U, std::size_t N, template<typename,typename> class VectorU, template<typename> class AllocU>
	friend class VectorComplexN;

	template<typename U = T> using Complex = std::complex<U>;
	template<typename U = T> using Alloc   = AlignedAllocType<U>;
	template<typename U = T> using Vector  = VectorType<U, AlignedAllocType<U>>;

	// Shorthand
	using ComplexVectorType = ComplexVector<T, VectorType, AlignedAllocType>;

	ComplexVector() {
		check();
	}

	ComplexVector(const unsigned size) {
		resize(size);
		check();
	}

	template<typename U>
	ComplexVector(const unsigned size, const Complex<U>& value) {
		resize(size);
		for ( unsigned i=0; i<size; ++i ) {
			setValue(i, value);
		}
		check();
	}

	// Returns 'i'th entry by value
	template<typename U = T>
	Complex<U> operator()(const unsigned i) const {
		return Complex<U>( real_[i], imag_[i] );
	}
	template<typename U>
	Complex<U> operator[](const unsigned i) const {
		return (*this)(i);  // use operator() to avoid duplication
	}


	// TODO Allow setting by reference
	template<typename U = T>
	class ElementRef {
	 public:
		ElementRef(const int i, ComplexVectorType& v):
			real_(v.real_[i]), imag_(v.imag_[i]) 
		{}

		// TODO Conversion to Complex
		template<typename V>
		ElementRef& operator=(const Complex<V>& value) {
			real_ = value.real();
			imag_ = value.imag();
			return *this;
		}
	 private:
		U& real_;
		U& imag_;
	};

	//ElementRef<T> operator()

	// Sets the 'i'th entry to 'value'
	template<typename U>
	void setValue(const unsigned i, const Complex<U> value) {
		real_[i] = value.real();
		imag_[i] = value.imag();
	}

	// Operations between vectors
	// - Syntax:  output = x <op> y
	template<typename U, template<typename, typename> class V, template<typename> class A>
	friend void add(const ComplexVector<U,V,A>& x, const ComplexVector<U,V,A>& y, ComplexVector<U,V,A>& output);
	template<typename U, template<typename, typename> class V, template<typename> class A>
	friend void subtract(const ComplexVector<U,V,A>& x, const ComplexVector<U,V,A>& y, ComplexVector<U,V,A>& output);
	template<typename U, template<typename, typename> class V, template<typename> class A>
	friend void multiply(const ComplexVector<U,V,A>& x, const ComplexVector<U,V,A>& y, ComplexVector<U,V,A>& output);
	template<typename U, template<typename, typename> class V, template<typename> class A>
	friend void divide(const ComplexVector<U,V,A>& x, const ComplexVector<U,V,A>& y, ComplexVector<U,V,A>& output);

	// TODO other operators, axpy, ...
	ComplexVectorType& operator+=(const ComplexVectorType& other) {
		add(*this, other, *this);
		return *this;
	}
	ComplexVectorType& operator-=(const ComplexVectorType& other) {
		subtract(*this, other, *this);
		return *this;
	}
	ComplexVectorType& operator*=(const ComplexVectorType& other) {
		multiply(*this, other, *this);
		return *this;
	}
	ComplexVectorType& operator/=(const ComplexVectorType& other) {
		divide(*this, other, *this);
		return *this;
	}

	// Binary ops	
	// - Note: this paradigm seems to have the best performance out of many options,
	//   but these binary ops (e.g. +, -) are still slower than in-place ops (e.g. +=, -=)
	//   because they require creation of a temporary
	friend ComplexVectorType operator+(ComplexVectorType x, const ComplexVectorType& y) {
		x += y;
		return x;
	}
	friend ComplexVectorType operator-(ComplexVectorType x, const ComplexVectorType& y) {
		x -= y;
		return x;
	}
	friend ComplexVectorType operator*(ComplexVectorType x, const ComplexVectorType& y) {
		x *= y;
		return x;
	}
	friend ComplexVectorType operator/(ComplexVectorType x, const ComplexVectorType& y) {
		x /= y;
		return x;
	}

	// Operations between a vector and a real-valued scalar
	// - Syntax:  output = x <op> value
	template<typename U, template<typename, typename> class V, template<typename> class A>
	friend void add(const ComplexVector<U,V,A>& x, const U value, ComplexVector<U,V,A>& output);

	// TODO
	ComplexVectorType& operator+=(const T value);
	ComplexVectorType& operator-=(const T value);
	ComplexVectorType& operator*=(const T value);
	ComplexVectorType& operator/=(const T value);

	// Operations between a vector and a complex-valued scalar
	// - Syntax:  output = x <op> value
	template<typename U, template<typename, typename> class V, template<typename> class A>
	friend void add(const ComplexVector<U,V,A>& x, const Complex<U>& value, ComplexVector<U,V,A>& output);
	template<typename U, template<typename, typename> class V, template<typename> class A>
	friend void multiply(const ComplexVector<U,V,A>& x, const Complex<U>& value, ComplexVector<U,V,A>& output);

	// TODO multiply on the left
	/*
	template<typename U, template<typename, typename> class V, template<typename> class A>
	friend void multiply(const Complex<U>& value, const ComplexVector<U,V,A>& x, ComplexVector<U,V,A>& output);
	*/

	ComplexVectorType& operator+=(const Complex<T>& value) {
		add(*this, value, *this);
		return *this;
	}
	ComplexVectorType& operator-=(const Complex<T>& value) {
		subtract(*this, value, *this);
		return *this;
	}
	ComplexVectorType& operator*=(const Complex<T>& value) {
		multiply(*this, value, *this);
		return *this;
	}
	ComplexVectorType& operator/=(const Complex<T>& value);  // TODO

	// TODO construct from std::vector<std::complex<U>>
	// TODO convert to std::vector<std::complex<U>> and return it
	// TODO operator[] that returns Complex<T>

	template<class OtherVectorType>
	OtherVectorType convert() const {
		int len = length();
		OtherVectorType other_vec;
		other_vec.resize(len);
		for ( int i=0; i<len; ++i ) {
			other_vec[i] = copyElement(i);
		}

		return other_vec;
	}

	void resize(const unsigned size) {
		real_.resize(size);
		imag_.resize(size);
	}

	int size() const {
		return real_.size();
	}
	int length() const {
		return size();
	}

	void debug_print_data_ptr() const {
		// TODO
		//std::cout << "DEBUG:  real_.data() = " << real_.data() << "\n";
		return;
	}

	//----- Access underlying data -----//

	// TODO
	/*
	template<typename U>
	class Accessor {
	 public:
		Accessor(ComplexVector<U>& vec, const unsigned i):
			vec_(vec), i_(i), real_(vec.real(i)), imag_(vec.imag(i))
		{}

	 private:
		ComplexVector<U>& vec_;
		unsigned i_;
		U& real_, imag_;
	}
	*/

	Vector<T>& real()  { return real_; }
	Vector<T>& imag()  { return imag_; }
	const Vector<T>& real() const { return real_; }
	const Vector<T>& imag() const { return imag_; }

	T& real(const unsigned i) { return real_[i]; }
	T& imag(const unsigned i) { return imag_[i]; }
	const T& real(const unsigned i) const { return real_[i]; }
	const T& imag(const unsigned i) const { return imag_[i]; }

 private:
	Vector<T> real_, imag_;

	template<typename U = T>
	Complex<U> copyElement(const unsigned i) const {
		return Complex<U>( real_[i], imag_[i] );
	}

	//----- Debugging -----//

	void check() const {
#ifdef DEBUG
		assert( real_.size() == imag_.size() );
		checkAlignment();
#endif //ifdef DEBUG
	}

	// Checks the alignment of the allocators at runtime
	// TODO Would be nice to be able to do this at compile time...
	void checkAlignment() const {
#ifdef DEBUG
		const int re_alignment = this->real_.get_allocator().get_alignment();
		const int im_alignment = this->imag_.get_allocator().get_alignment();
		assert( re_alignment == CACHE_LINE_SIZE );
		assert( im_alignment == CACHE_LINE_SIZE );
		assert( re_alignment == im_alignment );
#endif //ifdef DEBUG
	}

	// TODO make custom exception for ComplexVector class?
	void checkIndex(const unsigned i) const {
		if ( i >= size() ) {
			std::stringstream err_ss;
			err_ss << "index " << i << " is out of bounds (size: " << size() << ")\n";
			throw std::runtime_error( err_ss.str() );
		}
	}
};


//----- Vectors -----//

template<typename T, template<typename, typename> class V, template<typename> class A>
void add(const ComplexVector<T,V,A>& x, const ComplexVector<T,V,A>& y, ComplexVector<T,V,A>& output)
{
	const int size = x.size();
	output.resize(size);

	SIMD::Vector::add( x.real_.data(), y.real_.data(), output.real_.data(), size );
	SIMD::Vector::add( x.imag_.data(), y.imag_.data(), output.imag_.data(), size );
}


template<typename T, template<typename, typename> class V, template<typename> class A>
void subtract(const ComplexVector<T,V,A>& x, const ComplexVector<T,V,A>& y, ComplexVector<T,V,A>& output)
{
	const int size = x.size();
	output.resize(size);

	SIMD::Vector::subtract( x.real_.data(), y.real_.data(), output.real_.data(), size );
	SIMD::Vector::subtract( x.imag_.data(), y.imag_.data(), output.imag_.data(), size );
}


template<typename T, template<typename, typename> class V, template<typename> class A>
void multiply(const ComplexVector<T,V,A>& x, const ComplexVector<T,V,A>& y, ComplexVector<T,V,A>& output)
{
	const int size = x.size();
	output.resize(size);

	// Extract ptrs so that OpenMP can see that the data is cache-aligned (below)
	const T* re_x = x.real_.data();
	const T* im_x = x.imag_.data();
	const T* re_y = y.real_.data();
	const T* im_y = y.imag_.data();
	T* re_output = output.real_.data();
	T* im_output = output.imag_.data();

	#pragma omp simd aligned(re_x, im_x, re_y, im_y, re_output, im_output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		SIMD::Complex::multiply(re_x, im_x, re_y, im_y, re_output, im_output, i);
	}
}


template<typename T, template<typename, typename> class V, template<typename> class A>
void divide(const ComplexVector<T,V,A>& x, const ComplexVector<T,V,A>& y, ComplexVector<T,V,A>& output)
{
	const int size = x.size();
	output.resize(size);

	// Extract ptrs so that OpenMP can see that the data is cache-aligned (below)
	const T* re_x = x.real_.data();
	const T* im_x = x.imag_.data();
	const T* re_y = y.real_.data();
	const T* im_y = y.imag_.data();
	T* re_output = output.real_.data();
	T* im_output = output.imag_.data();

	#pragma omp simd aligned(re_x, im_x, re_y, im_y, re_output, im_output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		SIMD::Complex::divide(re_x, im_x, re_y, im_y, re_output, im_output, i);
	}

	/*
	double re_tmp, inv_norm_y_i_sq;
	#pragma omp simd aligned(re_x, im_x, re_y, im_y, re_output, im_output: CACHE_LINE_SIZE) private(re_tmp, inv_norm_y_i_sq)
	for ( int i=0; i<size; ++i ) {
		inv_norm_y_i_sq = 1.0/(re_y[i]*re_y[i] + im_y[i]*im_y[i]);
		re_tmp       = inv_norm_y_i_sq*(re_x[i]*re_y[i] + im_x[i]*im_y[i]);
		im_output[i] = inv_norm_y_i_sq*(im_x[i]*re_y[i] - re_x[i]*im_y[i]);
		re_output[i] = re_tmp;
	}
	*/
}


//----- ComplexVector and Complex -----//

// TODO general complex type
template<typename T, template<typename,typename> class V, template<typename> class A>
void add(const ComplexVector<T,V,A>& x, const std::complex<T>& value, ComplexVector<T,V,A>& output)
{
	const int size = x.size();
	output.resize(size);

	SIMD::Vector::add( x.real_.data(), value.real(), output.real_.data(), size );
	SIMD::Vector::add( x.imag_.data(), value.imag(), output.imag_.data(), size );
}


// TODO general complex type
template<typename T, template<typename,typename> class V, template<typename> class A>
void subtract(const ComplexVector<T,V,A>& x, const std::complex<T>& value, ComplexVector<T,V,A>& output)
{
	const int size = x.size();
	output.resize(size);

	SIMD::Vector::subtract( x.real_.data(), value.real(), output.real_.data(), size );
	SIMD::Vector::subtract( x.imag_.data(), value.imag(), output.imag_.data(), size );
}


// Right-multiply
// TODO general complex type
template<typename T, template<typename, typename> class V, template<typename> class A>
void multiply(const ComplexVector<T,V,A>& x, const std::complex<T>& value, ComplexVector<T,V,A>& output)
{
	const int size = x.size();
	output.resize(size);

	// Extract ptrs so that OpenMP can see that the data is cache-aligned (below)
	const T* re_x = x.real_.data();
	const T* im_x = x.imag_.data();
	T* re_output = output.real_.data();
	T* im_output = output.imag_.data();

	const T re_value = value.real();
	const T im_value = value.imag();

	#pragma omp simd aligned(re_x, im_x, re_output, im_output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		SIMD::Complex::multiply(re_x, im_x, re_value, im_value, re_output, im_output, i);
	}
}


// TODO general complex type
template<typename T, template<typename, typename> class V, template<typename> class A>
void divide(const ComplexVector<T,V,A>& x, const std::complex<T>& value, ComplexVector<T,V,A>& output)
{
	const int size = x.size();
	output.resize(size);

	// Extract ptrs so that OpenMP can see that the data is cache-aligned (below)
	const T* re_x = x.real_.data();
	const T* im_x = x.imag_.data();
	T* re_output = output.real_.data();
	T* im_output = output.imag_.data();

	const T re_value = value.real();
	const T im_value = value.imag();

	#pragma omp simd aligned(re_x, im_x, re_output, im_output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		SIMD::Complex::divide(re_x, im_x, re_value, im_value, re_output, im_output, i);
	}
}

} // end namespace Aligned

#endif // ifndef COMPLEX_VECTOR_H
