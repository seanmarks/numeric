
#pragma once
#ifndef ALIGNED_SIMD_REAL_H
#define ALIGNED_SIMD_REAL_H

#include <cmath>
#include <cstddef>

namespace numeric {  // numeric arrays only
namespace aligned {  // operations on cache-aligned data
namespace simd {     // SIMD kernels
namespace real {     // real arrays

// TODO: choose based on 
//  (1) available instruction set
//  (2) sizeof(T)
//  ... or allow the compiler to figure it out
//static constexpr int SIMD_LEN = 4;


//-------------------------//
//----- Vector Inputs -----//
//-------------------------//

// output = x <op> y
template<typename T> inline
void add(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] + y[i];
	}
}

template<typename T> inline
void subtract(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] - y[i];
	}
}

template<typename T> inline
void multiply(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i]*y[i];
	}
}

template<typename T> inline
void divide(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i]/y[i];
	}
}


// In-place:  output <op>= x

template<typename T> inline
void add_inplace(const T* CXX_RESTRICT x, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] += x[i];
	}
}

template<typename T> inline
void subtract_inplace(const T* CXX_RESTRICT x, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] -= x[i];
	}
}

template<typename T> inline
void multiply_inplace(const T* CXX_RESTRICT x, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] *= x[i];
	}
}

template<typename T> inline
void divide_inplace(const T* CXX_RESTRICT x, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] /= x[i];
	}
}


//-----------------------------//
//----- Vector and Scalar -----//
//-----------------------------//

// output = x <op> a

template<typename T> inline
void add(const T* CXX_RESTRICT x, const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] + a;
	}
}

template<typename T> inline
void subtract(const T* CXX_RESTRICT x, const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] - a;
	}
}

template<typename T> inline
void multiply(const T* x, const T a, const int size, T* output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = a*x[i];
	}
}

template<typename T> inline
void divide(const T* x, const T a, const int size, T* output)
{
	multiply(x, 1.0/a, size, output);
}


// In-place:  output <op>= a

template<typename T> inline
void add_inplace(const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] += a;
	}
}

template<typename T> inline
void subtract_inplace(const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] -= a;
	}
}

template<typename T> inline
void multiply_inplace(const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] *= a;
	}
}

template<typename T> inline
void divide_inplace(const T a, const int size, T* CXX_RESTRICT output)
{
	multiply_inplace(1.0/a, size, output);
}


//----- Mixed inputs -----//

// output = a*x + y
template<typename T> inline
void axpy(const T a, const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = a*x[i] + y[i];
	}
}

// output = a*x + b*y
template<typename T> inline
void axpby(const T a, const T* CXX_RESTRICT  x, const T b, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = a*x[i] + b*y[i];
	}
}


} // end namespace real
} // end namespace simd
} // end namespace aligned
} // end namespace numeric

#endif // ALIGNED_SIMD_REAL_H
