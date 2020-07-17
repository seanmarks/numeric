// CommonTypes
// - AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// - Namespace with commonly-used types and some related constants

#pragma once
#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <array>
#include <complex>
#include <iostream>
#include <vector>

#include "AlignedAllocator.h"
#include "Constants.h"

namespace CommonTypes {


//----- Scalar Types -----//

using Real    = double;
using Complex = std::complex<Real>;


//----- Array Types -----//

template<typename T, std::size_t N>
using Array = std::array<T,N>;

// (x,y,z) components
static constexpr int N_DIM  = 3;
static constexpr int X_DIM = 0;
static constexpr int Y_DIM = 1;
static constexpr int Z_DIM = 2;
static const std::array<std::string, N_DIM> axis_names = {{ "x", "y", "z" }};


using Real3    = Array<Real,3>;
using Complex3 = Array<Complex, N_DIM>;
using Int3     = Array<int,N_DIM>;

// Ranges
using Real2    = Array<Real,2>;
using Range    = Real2;
using Bounds3D = Array<Range, N_DIM>;

// Misc.
using Int2 = std::array<int,2>;


//----- Vector Types -----//

template<typename T>
using CacheAlignedAllocator = numeric::aligned::CacheAlignedAllocator<T>;

template<typename T>
using Vector = std::vector<T,CacheAlignedAllocator<T>>;
//template<typename T> using Vector = std::vector<T>;

using VectorComplex = Vector<Complex>;
using VectorReal    = Vector<Real>;
using ComplexVector = VectorComplex;
using RealVector    = VectorReal;

using VectorReal3    = Vector<Real3>;
using VectorComplex3 = Vector<Complex3>;


//----- Matrix types -----//

using Matrix3x3 = std::array<Real3, N_DIM>;
using Matrix    = Matrix3x3;

// Box matrix for a triclinic cell
// - "GROMACS" convention: box vectors arranged as [ a^T; b^T; c^T ]
using Box = Matrix3x3;


//----- Interface types -----//

// XdrFile
//using Rvec      = XdrFile::Rvec;
//using RvecArray = XdrFile::RvecArray;

} // end namespace CommonTypes


template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T,N>& x)
{
	os << "[";
	for ( unsigned i=0; i<N; ++i ) {
		if ( i > 0 ) { os << ", "; }
		os << x[i];
	}
	os << "]";
	return os;
}


template<typename T, std::size_t N, std::size_t M>
std::ostream& operator<<(std::ostream& os, const std::array<std::array<T,M>,N>& x)
{
	os << "[";
	for ( unsigned i=0; i<N; ++i ) {
		if ( i > 0 ) {
			os << " ";
		}
		os << x[i] << "\n";
	}
	os << "]";
	return os;
}


template<typename T, typename A>
std::ostream& operator<<(std::ostream& os, const std::vector<T,A>& v)
{
	unsigned N = v.size();
	os << "[";
	for ( unsigned i=0; i<N; ++i ) {
		if ( i > 0 ) {
			os << "\n";
		}
		os << v[i];
	}
	os << "]\n";
	return os;
}

#endif // ifndef COMMON_TYPES_H
