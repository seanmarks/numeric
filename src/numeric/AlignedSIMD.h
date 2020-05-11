/* Main header for including SIMD kernels and helper functions
 * - !!! See Aligned.h !!! about using these functions safely and correctly
 * - Uses OpenMP SIMD
 * - Functions assume cache-aligned data
 */

#pragma once
#ifndef ALIGNED_SIMD_H
#define ALIGNED_SIMD_H

#include <cmath>
#include <cstddef>

#include "Aligned.h"

// Look for a valid restrict keyword
// - If it's not available, use a dummy
// - Make sure to undef it before leaving file scope.
#ifndef CXX_RESTRICT
#  define ALIGNED_SIMD_RESTRICT
#  define CXX_RESTRICT ""
#endif

#include "simd/Real.h"
#include "simd/Complex.h"

// Undefine local macros
#ifdef ALIGNED_SIMD_RESTRICT
#  undef ALIGNED_SIMD_RESTRICT
#  undef CXX_RESTRICT
#endif

#endif // define ALIGNED_SIMD_H
