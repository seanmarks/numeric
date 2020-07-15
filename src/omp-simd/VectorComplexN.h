// SIMD-friendly vector of complex-valued arrays (of fixed length)
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// - Note: operators that return a VectorComplexN rely on copy elision
//   (NRVO) for performance

#ifndef VECTOR_COMPLEX_N_H
#define VECTOR_COMPLEX_N_H

#include <cassert>
#include <cstddef>
#include <complex>
#include <exception>
#include <stdexcept>

#include "CommonTypes.h"
#include "ComplexVector.h"

namespace Aligned {

template<
	typename T,
	std::size_t N,
	template<typename, typename> class VectorType = std::vector,
	template<typename> class AlignedAllocType = CacheAlignedAllocator
>
class VectorComplexN
{
 public:
	template<typename U = T> using Complex  = std::complex<U>;
	template<typename U = T> using ComplexN = std::array<Complex<U>, N>;
	
	template<typename U = T> using Alloc         = AlignedAllocType<U>;
	template<typename U = T> using Vector        = VectorType<U, AlignedAllocType<U>>;
	template<typename U = T> using ComplexVector = ComplexVector<U, VectorType, AlignedAllocType>;

	using VectorComplexNType = VectorComplexN<T, N, VectorType, AlignedAllocType>;

	VectorComplexN():
		size_(0), data_(0) {}

	VectorComplexN(const unsigned size) {
		resize(size);
	}

	// TODO assign(...)
	template<typename U = T>
	VectorComplexN(const unsigned size, const ComplexN<U>& value) {
		resize(size);
		for ( unsigned i=0; i<size; ++i ) {
			setArray(i, value);
		}
	}

	// Get/set size
	void resize(const unsigned size) {
		size_ = size;
		data_.resize( N*size_ );
	}
	unsigned size() const {
		return size_;
	}
	unsigned length() const {
		return size();
	}

	// Sets 'i'th entry
	template<typename U = T>
	void setArray(const unsigned i, const ComplexN<U>& value) {
		for ( unsigned n=0; n<N; ++n ) {
			data_.setValue(i, value[n]);
		}
	}

	// Returns 'i'th array
	template<typename U = T>
	ComplexN<U> operator()(const unsigned i) const {
		return copyArray(i);
	}
	template<typename U = T>
	ComplexN<U> operator[](const unsigned i) const {
		return copyArray(i);
	}

	// Returns 'd'th value of 'i'th array
	template<typename U = T>
	Complex<U> operator()(const unsigned i, const unsigned n) const {
		return copyValue(i,n);
	}

	// TODO use an iterator? but how would it handle returning references...?
	// wrapper object to data_? deref returns copy?
	// - Iterator idea
	//   - Contains vec iterators to real/imag components, with built-in stride
	//     - Stride through values: ++it  -->  real_it,imag_it += 1
	//     - Stride through arrays: ++it  -->  real_it,imag_it += N
	template<typename U = T>
	VectorComplexNType& multiplyChunk(const ComplexVector<U>& vec, const int offset = 0) {
		T* re_this = data_.real_.data();
		T* im_this = data_.imag_.data();
		const U* re_vec = vec.real_.data();
		const U* im_vec = vec.imag_.data();

		// Working variables
		int    j, j_first = N_*offset, j_end = j_first + N_;
		double re_vec_i, im_vec_i, re_this_tmp;
	
		const int vec_size    = vec.size();   // num values (type U) in 'vec'
		for ( int i=0; i<vec_size; ++i ) {
			// 'i'th element of 'vec'
			const double re_vec_i = re_vec[i];
			const double im_vec_i = im_vec[i];

			// (a + ib)*(c + id) = (a*c - b*d) + i(a*d + b*c)
			// - pragma on the inner loop seems to produce the best speedup
			#pragma omp simd aligned(re_this, im_this : CACHE_LINE_SIZE) private(j_first, j_end, re_this_tmp)
			for ( j=j_first; j<j_end; ++j ) {
				re_this_tmp = re_this[j]*re_vec_i - im_this[j]*im_vec_i;  // (a*c - b*d)
				im_this[j]  = re_this[j]*im_vec_i + im_this[j]*re_vec_i;  // (a*d + b*c)
				re_this[j]  = re_this_tmp;
			}
			j_first =  j_end;
			j_end   += N_;
		}

		return *this;
	}

 private:
	unsigned size_;          // = data_.size()/N: number of ComplexN in unflattened vector
	ComplexVector<T> data_;  // flattened array: size = N*size_

	static constexpr int N_ = static_cast<int>(N);  // OpenMP prefers signed ints

	// Returns the 'i'th array
	template<typename U = T>
	ComplexN<U> copyArray(const unsigned i) const {
		ComplexN<U> array;
		for ( unsigned n=0; n<N; ++n ) {
			array[n] = copyValue(i, n);
		}
		return array;
	}

	// Returns the 'n'th value of the 'i'th array
	template<typename U = T>
	Complex<U> copyValue(const unsigned i, const unsigned n) const {
		return data_(N*i + n);
	}
};

} // end namespace Aligned

#endif // ifndef VECTOR_COMPLEX_N_H
