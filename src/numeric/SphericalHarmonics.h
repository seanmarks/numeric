// SphericalHarmonics.h
//
// ABOUT: Calculate the spherical harmonics Y_l^(m) (eta,phi)
//
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// NOTE: 
//     - eta = cos(theta) where theta is the polar angle (range [0,pi])
//	   - phi = planar angle (range [0,2*pi])
//     - "n" before a number in a function name means "-" 
//			- e.g. "n1" for "-1" ("negative 1")
//
// DEVELOPMENT
//	   - inline all "small" functions?

#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

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
#include "Matrix.h"

class SphericalHarmonics
{
 public:
	SphericalHarmonics(
		const int  l,
		const bool do_fast_harmonics = true
	);


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
	using Complex3 = CommonTypes::Complex3;

	// Vectors
	using VectorReal     = CommonTypes::VectorReal;
	using VectorComplex  = CommonTypes::VectorComplex;
	using VectorReal3    = CommonTypes::VectorReal3;
	using VectorComplex3 = CommonTypes::VectorComplex3;
	using RealVector     = CommonTypes::RealVector;
	using ComplexVector  = CommonTypes::ComplexVector;


	//----- Spherical harmonics: Y_l_m -----//

	void calculate_Y_l(
		const Real3&    x, 
		const double    r, 
		const bool      need_derivatives,
		VectorComplex&  Y_l, 
		VectorComplex3& derivs_Y_l
	) const;

	void calculate(
		const VectorReal3&    x, 
		const Vector<double>& r, 
		const bool            need_derivatives,
		Matrix<Complex>&  Y_l, 
		Matrix<Complex3>& derivs_Y_l
	) const;

	// Compute Legendre polynomial P_l and first l derivatives
	// - Output format: p = { P_l, d/dx P_l, ... , d^l/dx^l P_l }
	void legendreP0(const double x, Vector<double>& p) const;
	void legendreP3(const double x, Vector<double>& p) const;
	void legendreP4(const double x, Vector<double>& p) const;
	void legendreP6(const double x, Vector<double>& p) const;

	void legendreP0(const Vector<double>& x, Matrix<double>& p) const;
	void legendreP3(const Vector<double>& x, Matrix<double>& p) const;
	void legendreP4(const Vector<double>& x, Matrix<double>& p) const;
	void legendreP6(const Vector<double>& x, Matrix<double>& p) const;


	//----- Get/Set Functions -----//

	void setHarmonicIndex(const int l);

	int getHarmonicIndex() {
		return l_;
	}

	int get_num_m_values() const {
		return num_m_values_;
	}
	
	bool do_fast_harmonics() const {
		return do_fast_harmonics_;
	}


	//----- Vector Operations -----//

	// Note: These functions are not intended to operate on arbitrary complex vectors.
	//       Rather, they are meant to operate on vectors with the same layout as 
	//       SphericalHarmonics vectors
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
	//----- Universal Constants -----//
	static constexpr double PI_ = constants::pi;

	//----- Coefficients -----//

	// Legendre polynomials 
	//   COEFF_P{l}       ==>  leading coefficient of P_l(x)
	//   COEFF_D{m}_P{l}  ==>  leading coefficient of d^m/dx^m P_l(x)

	//static constexpr double COEFF_P2  =  TODO

	static constexpr double COEFF_P3    = 1.0/2.0;
	static constexpr double COEFF_D1_P3 = 3.0/2.0;
	static constexpr double COEFF_D2_P3 = 15.0;
	static constexpr double COEFF_D3_P3 = 15.0;

	static constexpr double COEFF_P4    = 1.0/8.0;
	static constexpr double COEFF_D1_P4 = 5.0/2.0;
	static constexpr double COEFF_D2_P4 = 15.0/2.0;
	static constexpr double COEFF_D3_P4 = 105.0;
	static constexpr double COEFF_D4_P4 = 105.0;

	//static constexpr double COEFF_P5 = TODO

	static constexpr double COEFF_P6    = 1.0/16.0;
	static constexpr double COEFF_D1_P6 = 21.0/8.0;
	static constexpr double COEFF_D2_P6 = 105.0/8.0;
	static constexpr double COEFF_D3_P6 = 315.0/2.0;
	static constexpr double COEFF_D4_P6 = 945.0/2.0;
	static constexpr double COEFF_D5_P6 = 10395.0;
	static constexpr double COEFF_D6_P6 = 10395.0;


	// Spherical harmonics
	//   COEFF_Y_{l}_{m}  ==>  (-1)^m * sqrt((2*l + 1)/(4*pi) * (l - m)!/(l + m)!)
	// - TODO
	//   - Better way to go about this?
	//   - Could merge these coefficients with P_l ones
	const Vector<double> coeff_Y_0_ = {
		0.2820947917738780
	};

	const Vector<double> coeff_Y_3_ = {{
		0.7463526651802310,
		-0.2154534560761000,
		0.0681323650955522,
		-0.0278149215755189
	}};

	const Vector<double> coeff_Y_4_ = {{
		0.8462843753216340,
		-0.1892349391515120,
		0.0446031029038193,
		-0.0119206806752224,
		0.0042145970709046
	}};

	const Vector<double> coeff_Y_6_ = {{
		1.0171072362820500,
		-0.1569430538290060,
		0.0248148756521035,
		-0.0041358126086839,
		0.0007550926197968,
		-0.0001609862874555,
		0.0000464727381991
	}};


 private:
	int l_;             // harmonic index
	int num_m_values_;  // = 2*l + 1 (normal) or = l + 1 (do_fast_harmonics)

	// "Fast harmonics:" only compute Y_l,m for m >= 0
	bool do_fast_harmonics_;

	Vector<double> coeff_Y_l_;

	// Use to compute Legendre polynomials
	std::function< void(const double, Vector<double>&) > legendreP_function_;

	std::function< void(const Vector<double>&, Matrix<double>&) > legendreP_matrix_function_;

	// Buffers
	// - TODO: make use of these somehow?
	mutable Vector<double>  inv_r_, xhat_, yhat_, eta_;
	mutable Matrix<double>  p_x_;
	mutable Vector<Complex> zeta_;
	mutable Matrix<Complex> zeta_m_minus_1_;

	mutable Vector<double> x2_;

	mutable VectorReal3    deriv_eta_;
	mutable VectorComplex3 deriv_zeta_;
};

#endif // SPHERICAL_HARMONICS_H
