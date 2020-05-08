
#pragma once
#ifndef ALIGNED_SIMD_REAL_H
#define ALIGNED_SIMD_REAL_H

#include <cmath>
#include <cstddef>

namespace aligned {  // operations on cache-aligned data
namespace simd {     // SIMD kernels
namespace real {     // real arrays

// TODO: choose based on 
//  (1) available instruction set
//  (2) sizeof(T)
//  ... or allow the compiler to figure it out
static constexpr int SIMD_LEN = 4;

//----- Vector Inputs -----//

// output = x <op> y
template<typename T> inline
void add(const T* x, const T* y, const int size, T* output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] + y[i];
	}
}

template<typename T> inline
void subtract(const T* x, const T* y, const int size, T* output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] - y[i];
	}
}

template<typename T> inline
void multiply(const T* x, const T* y, const int size, T* output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i]*y[i];
	}
}

template<typename T> inline
void divide(const T* x, const T* y, const int size, T* output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i]/y[i];
	}
}


//----- Vector and Scalar -----//

// output = x <op> a

template<typename T> inline
void add(const T* x, const T a, const int size, T* output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] + a;
	}
}

template<typename T> inline
void subtract(const T* x, const T a, const int size, T* output)
{
	add(x, -a, size, output);
}

template<typename T> inline
void multiply(const T* x, const T a, const int size, T* output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = a*x[i];
	}
}

template<typename T> inline
void divide(const T* x, const T a, const int size, T* output)
{
	multiply(x, 1.0/a, size, output);
}


//----- Mixed inputs -----//

// output = a*x + y
template<typename T> inline
void axpy(const T a, const T* x, const T* y, const int size, T* output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = a*x[i] + y[i];
	}
}

// output = a*x + b*y
template<typename T> inline
void axpby(const T a, const T* x, const T b, const T* y, const int size, T* output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = a*x[i] + b*y[i];
	}
}


} // end namespace real
} // end namespace simd
} // end namespace aligned

#endif // ALIGNED_SIMD_REAL_H
