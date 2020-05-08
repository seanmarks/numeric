/* Main header for including SIMD kernels and helper functions
 * - Uses OpenMP SIMD
 * - Functions generally assume cache-aligned data
 */

#pragma once
#ifndef ALIGNED_SIMD_H
#define ALIGNED_SIMD_H

#include <cmath>
#include <cstddef>

#include "Aligned.h"

#include "simd/Real.h"
#include "simd/Complex.h"

#endif // define ALIGNED_SIMD_H
