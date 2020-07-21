// Legendre polynomials and their derivatives
//
// NOTE
// - COEFF_D{m}_P{l}  ==>  leading coefficient of d^m/dx^m P_l(x)

#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <cmath>
#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

#include "Assert.h"
#include "CommonTypes.h"
#include "Matrix.h"


class Legendre
{
 public:
	template<typename T> using Vector = CommonTypes::Vector<T>;
	template<typename T> using Matrix = numeric::Matrix<T, Vector<T>>;

	Legendre() {}

	virtual
	int getHarmonicIndex() const = 0;

	// Compute Legendre polynomial P_l and first l derivatives
	// - Output format:  p = { P_l, d/dx P_l, ... , d^l/dx^l P_l }
	virtual
	void calculate(const double x, Vector<double>& p) const = 0;

	// Matrix version of the above function
	virtual
	void calculate(
		const Vector<double>& x,  // (num_points)
		Matrix<double>& p  // (num_points, l+1)
	) const = 0;

	// Calculates the transpose of 'p' from the function above
	virtual
	void calculate_T(
		const Vector<double>& x,  // (num_points)
		Matrix<double>& p  // (l+1, num_points)
	) const;
};


// P_0
class LegendreP0 : public Legendre
{
 public:
	LegendreP0(): Legendre() {}

	virtual
	void calculate(const double x, Vector<double>& p) const override;

	virtual
	void calculate(const Vector<double>& x, Matrix<double>& p) const override;

	virtual
	void calculate_T(const Vector<double>& x, Matrix<double>& p) const override;

	virtual
	int getHarmonicIndex() const override {
		return l_;
	}

 private:
	static constexpr int l_ = 0;

	static constexpr double COEFF_D0_P0 = 1.0/1.0;
};


// P_3
class LegendreP3 : public Legendre
{
 public:
	LegendreP3(): Legendre() {}

	virtual
	void calculate(const double x, Vector<double>& p) const override;

	virtual
	void calculate(const Vector<double>& x, Matrix<double>& p) const override;

	virtual
	void calculate_T(const Vector<double>& x, Matrix<double>& p) const override;

	virtual
	int getHarmonicIndex() const override {
		return l_;
	}

 private:
	static constexpr int l_ = 3;

	static constexpr double COEFF_D0_P3 = 1.0/2.0;
	static constexpr double COEFF_D1_P3 = 3.0/2.0;
	static constexpr double COEFF_D2_P3 = 15.0;
	static constexpr double COEFF_D3_P3 = 15.0;

	mutable Vector<double> x2_;
};


// P_4
class LegendreP4 : public Legendre
{
 public:
	LegendreP4(): Legendre() {}

	virtual
	void calculate(const double x, Vector<double>& p) const override;

	virtual
	void calculate(const Vector<double>& x, Matrix<double>& p) const override;

	virtual
	void calculate_T(const Vector<double>& x, Matrix<double>& p) const override;

	virtual
	int getHarmonicIndex() const override {
		return l_;
	}

 private:
	static constexpr int l_ = 4;

	static constexpr double COEFF_D0_P4 = 1.0/8.0;
	static constexpr double COEFF_D1_P4 = 5.0/2.0;
	static constexpr double COEFF_D2_P4 = 15.0/2.0;
	static constexpr double COEFF_D3_P4 = 105.0;
	static constexpr double COEFF_D4_P4 = 105.0;

	mutable Vector<double> x2_;
};


// P_6
class LegendreP6 : public Legendre
{
 public:
	LegendreP6(): Legendre() {}

	virtual
	void calculate(const double x, Vector<double>& p) const override;

	virtual
	void calculate(const Vector<double>& x, Matrix<double>& p) const override;

	virtual
	void calculate_T(const Vector<double>& x, Matrix<double>& p) const override;

	virtual
	int getHarmonicIndex() const override {
		return l_;
	}

 private:
	static constexpr int l_ = 6;

	static constexpr double COEFF_D0_P6 = 1.0/16.0;
	static constexpr double COEFF_D1_P6 = 21.0/8.0;
	static constexpr double COEFF_D2_P6 = 105.0/8.0;
	static constexpr double COEFF_D3_P6 = 315.0/2.0;
	static constexpr double COEFF_D4_P6 = 945.0/2.0;
	static constexpr double COEFF_D5_P6 = 10395.0;
	static constexpr double COEFF_D6_P6 = 10395.0;

	mutable Vector<double> x2_;
};

#endif // ifndef LEGENDRE_H
