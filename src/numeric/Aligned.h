/* Aligned.h
 * !!! README !!!
 * - Convenience header for creating and performing arithmetic operations
 *   on cache-aligned vectors
 * - The 'aligned' namespace **assumes** that all vector-like objects (i.e. variable-length storage)
 *   are allocatedto cache-aligned blocks of memory. This is for optimal performance.
 *   - If you use non-aligned vector-like objects with the SIMD kernels, the results will be **WRONG!**
 *   - The OpenMD SIMD kernals are written with the explicit requirement of cache alignment
 *   - In short: you should use aligned:: objects with aligned::simd kernels
 */

#pragma once
#ifndef ALIGNED_H
#define ALIGNED_H

#include "AlignedAllocator.h"
#include "AlignedSIMD.h"

#endif // ifndef ALIGNED_H
