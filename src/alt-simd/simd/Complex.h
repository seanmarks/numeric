// SIMD kernels and helper functions for complex vectors
// - Use OpenMP SIMD
// - Assume cache-aligned data
//
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)


#pragma once
#ifndef ALIGNED_SIMD_COMPLEX_H
#define ALIGNED_SIMD_COMPLEX_H

#include <cmath>
#include <cstddef>

namespace numeric {  // numeric arrays only
namespace aligned {  // operations on cache-aligned data

namespace simd    {  // SIMD kernels
namespace complex {  // complex vectors


//-------------------------//
//----- Vector Inputs -----//
//-------------------------//

// output = u + v
template<typename T> inline
void add(
	const T* CXX_RESTRICT u_re, const T* CXX_RESTRICT u_im, const T* CXX_RESTRICT v_re, const T* CXX_RESTRICT v_im, const int size,
	T* CXX_RESTRICT output_re, T* CXX_RESTRICT output_im
)
{
	#pragma omp simd aligned(u_re, v_re, output_re: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output_re[i] = u_re[i] + v_re[i];
	}
	#pragma omp simd aligned(u_im, v_im, output_im: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output_im[i] = u_im[i] + v_im[i];
	}
}

// output = u - v
template<typename T> inline
void subtract(
	const T* CXX_RESTRICT u_re, const T* CXX_RESTRICT u_im, const T* CXX_RESTRICT v_re, const T* CXX_RESTRICT v_im, const int size,
	T* CXX_RESTRICT output_re, T* CXX_RESTRICT output_im
)
{
	#pragma omp simd aligned(u_re, v_re, output_re: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output_re[i] = u_re[i] - v_re[i];
	}
	#pragma omp simd aligned(u_im, v_im, output_im: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output_im[i] = u_im[i] - v_im[i];
	}
}

// output = u * v (element-wise)
template<typename T> inline
void multiply(
	const T* CXX_RESTRICT u_re, const T* CXX_RESTRICT u_im, const T* CXX_RESTRICT v_re, const T* CXX_RESTRICT v_im, const int size,
	T* CXX_RESTRICT output_re, T* CXX_RESTRICT output_im
)
{
	#pragma omp simd aligned(u_re, u_im, v_re, v_im, output_re, output_im: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output_re[i] = u_re[i]*v_re[i] - u_im[i]*v_im[i];
		output_im[i] = u_re[i]*v_im[i] + u_im[i]*v_re[i];
	}
}

// output = u / v (element-wise)
template<typename T> inline
void divide(const T* CXX_RESTRICT u, const T* CXX_RESTRICT v, const int size, T* CXX_RESTRICT output)
{
	/*
	const int num_values = 2*size;
	int im = 1;
	T inv_norm_v_sq;
	#pragma omp simd aligned(u, v, output: CACHE_LINE_SIZE) private(im, inv_norm_v_sq)
	for ( int re=0; re<num_values; re+=2 ) {
		im = re + 1;
		inv_norm_v_sq = 1.0/(v[re]*v[re] + v[im]*v[im]);  // 1/|v|^2
		output[re] = ( u[re]*v[re] + u[im]*v[im] )*inv_norm_v_sq;
		output[im] = ( u[im]*v[re] - u[re]*v[im] )*inv_norm_v_sq;
	}
	*/
}

// TODO more in-place ops
// output += u
template<typename T> inline
void add_in_place(const T* CXX_RESTRICT u, const int size, T* CXX_RESTRICT output)
{
	//simd::real::add_in_place(u, 2*size, output);
}

// output -= u
template<typename T> inline
void subtract_in_place(const T* CXX_RESTRICT u, const int size, T* CXX_RESTRICT output)
{
	//simd::real::subtract_in_place(u, 2*size, output);
}



//-----------------------------//
//----- Vector and Scalar -----//
//-----------------------------//

// output = alpha*output
template<typename T> inline
void left_multiply_in_place(const T alpha_re, const T alpha_im, const int size, T* CXX_RESTRICT output)
{
	/*
	const int num_values = 2*size;
	int im = 1;
	T   tmp_re;
	#pragma omp simd aligned(output: CACHE_LINE_SIZE) private(im, tmp_re)
	for ( int re=0; re<num_values; re+=2 ) {
		im = re + 1; 
		tmp_re     = alpha_re*output[re] - alpha_im*output[im];
		output[im] = alpha_re*output[im] + alpha_im*output[re];
		output[re] = tmp_re;
	}
	*/
}


/*
// output = x * y (element-wise)
template<typename T> inline
void multiply(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i]*y[i];
	}
}

// output = x / y (element-wise)
template<typename T> inline
void divide(const T* CXX_RESTRICT x, const T* CXX_RESTRICT y, const int size, T* CXX_RESTRICT output)
{
	#pragma omp simd aligned(x, y, output: CACHE_LINE_SIZE)
	for ( int i=0; i<size; ++i ) {
		output[i] = x[i]/y[i];
	}
}
*/


} // end namespace complex
} // end namespace simd

} // end namespace aligned
} // end namespace numeric

#endif // define ALIGNED_SIMD_COMPLEX_H
