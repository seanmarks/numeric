/*
 * SIMD kernels and helper functions for complex vectors
 * - Use OpenMP SIMD
 * - Assume cache-aligned data
 * - Complex vectors have the following layout
 *   - size       = length of complex vectors = number of complex values of the form a+ib
 *   - num_values = 2*size = total number of real numbers required to represent it
 */


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
void add(const T* CXX_RESTRICT u, const T* CXX_RESTRICT v, const int size, T* CXX_RESTRICT output)
{
	simd::real::add(u, v, 2*size, output);
}

// output = u - v
template<typename T> inline
void subtract(const T* CXX_RESTRICT u, const T* CXX_RESTRICT v, const int size, T* CXX_RESTRICT output)
{
	simd::real::subtract(u, v, 2*size, output);
}

// output = u * v (element-wise)
template<typename T> inline
void multiply(const T* CXX_RESTRICT u, const T* CXX_RESTRICT v, const int size, T* CXX_RESTRICT output)
{
	const int num_values = 2*size;
	int im = 1;
	#pragma omp simd aligned(u, v, output: CACHE_LINE_SIZE) private(im)
	for ( int re=0; re<num_values; re+=2 ) {
		im = re + 1; 
		output[re] = u[re]*v[re] - u[im]*v[im];
		output[im] = u[re]*v[im] + u[im]*v[re];
	}
}

// output = u / v (element-wise)
template<typename T> inline
void divide(const T* CXX_RESTRICT u, const T* CXX_RESTRICT v, const int size, T* CXX_RESTRICT output)
{
	const int num_values = 2*size;
	int im = 1;
	T inv_norm_v_sq;
	#pragma omp simd aligned(u, v, output: CACHE_LINE_SIZE) private(im, inv_norm_v_sq)
	for ( int re=0; re<num_values; re+=2 ) {
		im = re + 1;
		inv_norm_v_sq = 1.0/(v[re]*v[re] + v[im]*v[im]);  // 1/|v|^2
		output[re] = ( u[re]*v[re] + u[im]*v[im] )*inv_norm_v_sq;
		output[im] = ( u[im]*v[re] - u[re]*v[im] )*inv_norm_v_sq;

		//output[re] = ( u[re]*v[re] - u[im]*v[im] )*inv_norm_v_sq;
		//output[im] = ( u[im]*v[re] - u[re]*v[im] )*inv_norm_v_sq;
	}
}

// TODO in-place


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




// TODO
// - Rename functions based on specifiers?
//   - e.g. 'multiply'  -->  multiply_simple ('uniform(...) linear(i: 1) notinbranch')
//                           multiply_contiguous

/*
// output = x*y
#pragma omp declare simd \
	aligned(re_x, im_x, re_y, im_y, re_output, im_output: CACHE_LINE_SIZE) \
	uniform(re_x, im_x, re_y, im_y, re_output, im_output) \
	linear(i: 1) simdlen(ALIGNED_SIMD_LEN) notinbranch
template<typename T> inline
void multiply(const T* re_x, const T* im_x, const T* re_y, const T* im_y,  
              T* re_output, T* im_output, const int i)
{
	T re_tmp     = re_x[i]*re_y[i] - im_x[i]*im_y[i];
	im_output[i] = re_x[i]*im_y[i] + im_x[i]*re_y[i];
	re_output[i] = re_tmp;
}

// output = x/y
#pragma omp declare simd \
	aligned(re_x, im_x, re_y, im_y, re_output, im_output: CACHE_LINE_SIZE) \
	uniform(re_x, im_x, re_y, im_y, re_output, im_output) \
	linear(i: 1) simdlen(ALIGNED_SIMD_LEN) notinbranch
template<typename T> inline
void divide(const T* re_x, const T* im_x, const T* re_y, const T* im_y,  
              T* re_output, T* im_output, const int i)
{
	T inv_norm_y_i_sq = 1.0/(re_y[i]*re_y[i] + im_y[i]*im_y[i]);
	T re_tmp     = inv_norm_y_i_sq*(re_x[i]*re_y[i] + im_x[i]*im_y[i]);
	im_output[i] = inv_norm_y_i_sq*(im_x[i]*re_y[i] - re_x[i]*im_y[i]);
	re_output[i] = re_tmp;
}

// Right-multiply: output = x*a
#pragma omp declare simd \
	aligned(re_x, im_x, re_output, im_output: CACHE_LINE_SIZE) \
	uniform(re_x, im_x, re_output, im_output) \
	linear(i: 1) simdlen(ALIGNED_SIMD_LEN) notinbranch
template<typename T> inline
void multiply(const T* re_x, const T* im_x, const T re_a, const T im_a,
              T* re_output, T* im_output, const int i)
{
	T re_tmp     = re_x[i]*re_a - im_x[i]*im_a;
	im_output[i] = re_x[i]*im_a + im_x[i]*re_a;
	re_output[i] = re_tmp;
}

// Left-multiply: output = a*x
#pragma omp declare simd \
	aligned(re_x, im_x, re_output, im_output: CACHE_LINE_SIZE) \
	uniform(re_x, im_x, re_output, im_output) \
	linear(i: 1) simdlen(ALIGNED_SIMD_LEN) notinbranch
template<typename T> inline
void multiply(const T re_a, const T im_a, const T* re_x, const T* im_x, 
              T* re_output, T* im_output, const int i)
{
	T re_tmp     = re_a*re_x[i] - im_a*im_x[i];
	im_output[i] = re_a*im_x[i] + im_a*re_x[i];
	re_output[i] = re_tmp;
}
 
// Division: output = x/a = x*(1/a)
#pragma omp declare simd \
	aligned(re_x, im_x, re_output, im_output: CACHE_LINE_SIZE) \
	uniform(re_x, im_x, re_output, im_output) \
	linear(i: 1) simdlen(ALIGNED_SIMD_LEN) notinbranch
template<typename T> inline
void divide(const T* re_x, const T* im_x, const T re_a, const T im_a,
            T* re_output, T* im_output, const int i)
{
	T inv_norm_a_sq = 1.0/(re_a*re_a + im_a*im_a); // TODO: option to pass as input param? wasteful to recompute
	T re_tmp     = inv_norm_a_sq*(re_x[i]*re_a + im_x[i]*im_a);
	im_output[i] = inv_norm_a_sq*(im_x[i]*re_a - re_x[i]*im_a);
	re_output[i] = re_tmp;
}
*/

} // end namespace complex
} // end namespace simd

} // end namespace aligned
} // end namespace numeric

#endif // define ALIGNED_SIMD_COMPLEX_H
