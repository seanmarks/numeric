// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef ALIGNED_SIMD_REAL_H
#define ALIGNED_SIMD_REAL_H

#include <cmath>
#include <cstddef>

// TODO: choose based on 
//  (1) available instruction set
//  (2) sizeof(T)
//  ... or allow the compiler to figure it out
//static constexpr int SIMD_LEN = 4;

namespace numeric {  // numeric arrays only
namespace aligned {  // operations on cache-aligned data

namespace simd {  // SIMD kernels
namespace real {  // real vectors


//-------------------------//
//----- Vector Inputs -----//
//-------------------------//

// output = x + y
template<typename T> inline
void add(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] + y[i];
	}
}

// output = x - y
template<typename T> inline
void subtract(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] - y[i];
	}
}

// output = x * y (element-wise)
template<typename T> inline
void multiply(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i]*y[i];
	}
}

// output = x / y (element-wise)
template<typename T> inline
void divide(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE) //simdlen(SIMD_LEN)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i]/y[i];
	}
}

// output += x
template<typename T> inline
void add_in_place(const T* CXX_RESTRICT x, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] += x[i];
	}
}

// output -= x
template<typename T> inline
void subtract_in_place(const T* CXX_RESTRICT x, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] -= x[i];
	}
}

// output *= x  (element-wise)
template<typename T> inline
void multiply_in_place(const T* CXX_RESTRICT x, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] *= x[i];
	}
}

// output *= x  (element-wise)
template<typename T> inline
void divide_in_place(const T* CXX_RESTRICT x, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] /= x[i];
	}
}


//-----------------------------//
//----- Vector and Scalar -----//
//-----------------------------//

// output = x + a
template<typename T> inline
void add(const T* CXX_RESTRICT x, const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] + a;
	}
}

// output = x - a
template<typename T> inline
void subtract(const T* CXX_RESTRICT x, const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i] - a;
	}
}

// output = a*x
template<typename T> inline
void multiply(const T* CXX_RESTRICT x, const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = a*x[i];
	}
}

// output = x/a
template<typename T> inline
void divide(const T* CXX_RESTRICT x, const T a, const int size, T* CXX_RESTRICT output)
{
	multiply(x, 1.0/a, size, output);
}

// output += a
template<typename T> inline
void add_in_place(const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] += a;
	}
}

// output -= a
template<typename T> inline
void subtract_in_place(const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] -= a;
	}
}

// output *= a
template<typename T> inline
void multiply_in_place(const T a, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] *= a;
	}
}

// output /= a
template<typename T> inline
void divide_in_place(const T a, const int size, T* CXX_RESTRICT output)
{
	multiply_in_place(1.0/a, size, output);
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
