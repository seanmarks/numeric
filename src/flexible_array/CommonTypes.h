/* CommonTypes
 *
 * Namespace with commonly-used types
 *
 */

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <array>
#include <complex>
#include <vector>

namespace CommonTypes {

//----- Scalar Types -----//

using Real    = double;
using Complex = std::complex<Real>;


//----- Vector Types -----//

// std::vector with default allocator
template<typename ValueType>
using Vector = std::vector<ValueType>;

// (x,y,z) components
static const int DIM_  = 3;
static const int X_DIM = 0;
static const int Y_DIM = 1;
static const int Z_DIM = 2;
static const std::array<std::string, DIM_> axis_names = {{ "x", "y", "z" }};

// 1-D range
using Real2 = std::array<Real,2>;
using Range = Real2;

// 3-D range
using Bounds3D = std::array<Range, DIM_>;

/*
// Raw C arrays of floats (for interfacing with GROMACS-related codes)
#ifndef COMMON_TYPES_PLUMED_MODE
using Rvec      = float[DIM_];
#using RvecArray = SimpleVector<Rvec>;
#else
using Rvec      = PLMD::Vector;
using RvecArray = std::vector<Rvec>;  // std::vector of PLMD::vector
#endif
*/

// STL-derived array types of length 3
using Real3    = std::array<Real,DIM_>;
using Int3     = std::array<int,DIM_>;
using Complex3 = std::array<Complex, DIM_>;

// Vectors of indefinite length
//using VectorReal     = Vector<Real>;
//using VectorReal3    = Vector<Real3>;
//using VectorComplex  = Vector<Complex>;
//using VectorComplex3 = Vector<Complex3>;


//----- Matrix types -----//

using Matrix3x3 = std::array<Real3, DIM_>;
using Matrix    = Matrix3x3;

// Box matrix for a triclinic cell
// - "GROMACS" convention: box vectors arranged as [ a^T; b^T; c^T ]
using Box = Matrix3x3;


//----- Zeros -----//

const Real3 zero_3 = {{ 0.0, 0.0, 0.0 }};

const Matrix zero_matrix = {{ zero_3, zero_3, zero_3 }};

} // namespace CommonTypes

#endif // ifndef COMMON_TYPES_H
