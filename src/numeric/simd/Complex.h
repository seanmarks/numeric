// SIMD kernels and helper functions
// - Uses OpenMP SIMD
// - Functions generally assume cache-aligned data

#pragma once
#ifndef ALIGNED_SIMD_H
#define ALIGNED_SIMD_H

#include <cmath>
#include <cstddef>

// TODO how best to choose?
#ifndef ALIGNED_SIMD_LEN 
#  define ALIGNED_SIMD_LEN 16
#endif

namespace numeric {
namespace aligned {
namespace simd    {

namespace complex {

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

#endif // define ALIGNED_SIMD_H
