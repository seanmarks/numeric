#include "Legendre.h"


void LegendreP0::calculate(const double x, Vector<double>& p) const
{
	p.assign(1, COEFF_D0_P0);
}


void LegendreP0::calculate(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	p.assign({{num_points, 1}}, COEFF_D0_P0);
}




void LegendreP3::calculate(const double x, Vector<double>& p) const
{
	p.resize(4);

	double x2 = x*x;
	p[0] = COEFF_D0_P3*x*(5.0*x2 - 3.0);
	p[1] = COEFF_D1_P3*(5.0*x2 - 1.0);
	p[2] = COEFF_D2_P3*x;
	p[3] = COEFF_D3_P3;
}


void LegendreP3::calculate(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	static constexpr int num_cols = 4;
	p.resize(num_points, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];
		p(i,0) = COEFF_D0_P3*x[i]*(5.0*x2 - 3.0);
		p(i,1) = COEFF_D1_P3*(5.0*x2 - 1.0);
		p(i,2) = COEFF_D2_P3*x[i];
		p(i,3) = COEFF_D3_P3;
	}
}




void LegendreP4::calculate(const double x, Vector<double>& p) const
{
	p.resize(5);

	double x2 = x*x;
	p[0] = COEFF_D0_P4*(x2*(35.0*x2 - 30.0) + 3.0);
	p[1] = COEFF_D1_P4*x*(7.0*x2 - 3.0);
	p[2] = COEFF_D2_P4*(7.0*x2 - 1.0);
	p[3] = COEFF_D3_P4*x;
	p[4] = COEFF_D4_P4;
}


void LegendreP4::calculate(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	static constexpr int num_cols = 5;
	p.resize(num_points, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];
		p(i,0) = COEFF_D0_P4*(x2*(35.0*x2 - 30.0) + 3.0);
		p(i,1) = COEFF_D1_P4*x[i]*(7.0*x2 - 3.0);
		p(i,2) = COEFF_D2_P4*(7.0*x2 - 1.0);
		p(i,3) = COEFF_D3_P4*x[i];
		p(i,4) = COEFF_D4_P4;
	}
}




void LegendreP6::calculate(const double x, Vector<double>& p) const
{
	p.resize(7);

	double x2 = x*x;
	p[0] = COEFF_D0_P6*(x2*(x2*(231.0*x2 - 315.0) + 105.0) - 5.0);
	p[1] = COEFF_D1_P6*x*(x2*(33.0*x2 - 30.0) + 5.0);
	p[2] = COEFF_D2_P6*(x2*(33.0*x2 - 18.0) + 1.0);
	p[3] = COEFF_D3_P6*x*(11.0*x2 - 3.0);
	p[4] = COEFF_D4_P6*(11.0*x2 - 1.0);
	p[5] = COEFF_D5_P6*x;
	p[6] = COEFF_D6_P6;
}


void LegendreP6::calculate(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	static constexpr int num_cols = 7;
	p.resize(num_points, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];
		p(i,0) = COEFF_D0_P6*(x2*(x2*(231.0*x2 - 315.0) + 105.0) - 5.0);
		p(i,1) = COEFF_D1_P6*x[i]*(x2*(33.0*x2 - 30.0) + 5.0);
		p(i,2) = COEFF_D2_P6*(x2*(33.0*x2 - 18.0) + 1.0);
		p(i,3) = COEFF_D3_P6*x[i]*(11.0*x2 - 3.0);
		p(i,4) = COEFF_D4_P6*(11.0*x2 - 1.0);
		p(i,5) = COEFF_D5_P6*x[i];
		p(i,6) = COEFF_D6_P6;
	}
}
