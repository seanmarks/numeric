// RealSphericalHarmonics.h
//
// ABOUT: Calculate the real spherical harmonics, y_{l,m}(phi,eta)
//
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// NOTE: 
//	   - phi = planar angle (range [0,2*pi])
//     - eta = cos(theta), where theta is the polar angle (range [0,pi])
//     - "n" before a number in a function name means "-" 
//			- e.g. "n1" for "-1" ("negative 1")
//
// DEVELOPMENT
//	   - inline all "small" functions?

#ifndef REAL_SPHERICAL_HARMONICS_H
#define REAL_SPHERICAL_HARMONICS_H

#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Assert.h"
#include "CommonTypes.h"
#include "Constants.h"
#include "Legendre.h"
#include "Matrix.h"


class RealSphericalHarmonics
{
 public:
	//----- Constants and Typedefs -----//

  static constexpr int X_DIM = CommonTypes::X_DIM;
  static constexpr int Y_DIM = CommonTypes::Y_DIM;
  static constexpr int Z_DIM = CommonTypes::Z_DIM;
  static constexpr int N_DIM = CommonTypes::N_DIM;  // Dimensionality of simulation 

	template<typename T>
	using Vector = CommonTypes::Vector<T>;
	template<typename T>
	using Matrix = numeric::Matrix<T, Vector<T>>;

	// Scalars
	using Real    = CommonTypes::Real;
	using Complex = CommonTypes::Complex;

	// Arrays
	using Real3    = CommonTypes::Real3;

	// Vectors
	using VectorReal     = CommonTypes::VectorReal;
	using VectorComplex  = CommonTypes::VectorComplex;
	using VectorReal3    = CommonTypes::VectorReal3;
	using VectorComplex3 = CommonTypes::VectorComplex3;
	using RealVector     = CommonTypes::RealVector;
	using ComplexVector  = CommonTypes::ComplexVector;


	//----- Constructors -----//

	RealSphericalHarmonics(const int l);


	//----- Spherical harmonics: y_{l,m} -----//

	void calculate(
		const VectorReal3&    x,  // (num_points,3)
		const Vector<double>& r,  // (num_points)
		const bool            need_derivatives,
		Matrix<double>& y_l,        // (num_points, l+1)
		Matrix<Real3>&  derivs_y_l  // FIXME: 3 separate matrices?
	) const;


	//----- Get/Set Functions -----//

	void setHarmonicIndex(const int l);

	int getHarmonicIndex() {
		return l_;
	}

	int get_num_m_values() const {
		return num_m_values_;
	}


	//----- Vector Operations -----//

	/*
	// Note: These functions are not intended to operate on arbitrary complex vectors.
	//       Rather, they are meant to operate on vectors with the same layout as 
	//       RealSphericalHarmonics vectors
	// - i.e. length 2*l+1 or l+1 (fast_harmonics), with values for m<0 at the end

	// Norm of a complex vector: 
	//   sqrt(<u,u>) where < , > is the Hermitian inner product
	double norm(const VectorComplex& vec) const;

	// Hermitian inner product
	Complex innerProduct(const VectorComplex& x, const VectorComplex& y) const;

	// Hermitian inner product, where each element of 'x' multiplies the complex conjugate 
	// of each of the DIM entries in the corresponding element of 'dy'
	Complex3 innerProduct(const VectorComplex& x, const VectorComplex3& dy) const;
	Complex3 innerProduct(
		const VectorComplex& x, 
		const VectorComplex3::const_iterator& dy_first,
		const VectorComplex3::const_iterator& dy_last    // one past the last element
	) const;
	*/


	//----- Misc. -----//

	// Spherical coordinates, but using eta = cos(theta) instead of theta
	void calculateModifiedSphericalCoordinates(
		const Real3& x,
		// Output
		double& r, 
		double& phi, // Planar angle (range [0, 2*pi])
		double& eta  // cos(theta), where theta = azimuthal angle (range [0,pi])
	) const;


 protected:
	// Calculates n! recursively
	static constexpr std::size_t factorial(const std::size_t n) {
		return n > 0 ? n*factorial(n-1) : 1;
	}

 private:
	int l_;             // harmonic index
	int num_m_values_;  // = 2*l + 1 (normal) or = l + 1 (do_fast_harmonics)

	// Coefficients
	//   COEFF_Y_{l}_{m}  ==>  (-1)^m * sqrt((2*l + 1)/(4*pi) * (l - m)!/(l + m)!)
	Vector<double> coeff_y_l_;

	// Use to compute Legendre polynomials and their derivatives
	std::unique_ptr<Legendre> legendre_ptr_ = nullptr;

	// Buffers
	// - TODO: make use of these somehow?
	mutable Vector<double> inv_r_, xhat_, yhat_, eta_;
	//mutable Vector<double> re_zeta_,  im_zeta_;
	mutable std::array<Vector<double>, N_DIM> re_deriv_zeta_, im_deriv_zeta_, deriv_eta_;
	//mutable VectorReal3    re_deriv_zeta_, im_deriv_zeta_, deriv_eta_;

	mutable Matrix<double>  p_x_;  // Legendre polynomials (and derivatives)

	/*
	mutable Vector<Complex> zeta_;
	mutable Matrix<Complex> zeta_m_minus_1_;

	mutable VectorComplex3 deriv_zeta_;
	*/
};

#endif // REAL_SPHERICAL_HARMONICS_H
