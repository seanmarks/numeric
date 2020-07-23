// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "Legendre.h"

void Legendre::calculate_T(const Vector<double>& x, Matrix<double>& p) const
{
	Matrix<double> p_tmp;
	calculate(x, p_tmp);

	int num_points   = x.size();
	int num_m_values = getHarmonicIndex() + 1;
	p.resize(num_m_values, num_points);
	for ( int m=0; m<num_m_values; ++m ) {
		for ( int i=0; i<num_points; ++i ) {
			p(m,i) = p_tmp(i,m);
		}
	}
}


// GCC 7.5.0 posts link errors if this definition is not included
constexpr double LegendreP0::COEFF_D0_P0;


void LegendreP0::calculate(const double x, Vector<double>& p) const
{
	p.assign(l_ + 1, COEFF_D0_P0);
}


void LegendreP0::calculate(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	p.assign({{num_points, l_ + 1}}, COEFF_D0_P0);
}


void LegendreP0::calculate_T(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	p.assign({{l_ + 1, num_points}}, COEFF_D0_P0);
}




void LegendreP3::calculate(const double x, Vector<double>& p) const  // MISSED?
{
	p.resize(l_ + 1);

	double x2 = x*x;
	p[0] = COEFF_D0_P3*x*(5.0*x2 - 3.0);
	p[1] = COEFF_D1_P3*(5.0*x2 - 1.0);
	p[2] = COEFF_D2_P3*x;
	p[3] = COEFF_D3_P3;
}  // MISSED?


void LegendreP3::calculate(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	static constexpr int num_cols = l_ + 1;
	p.resize(num_points, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for simd schedule(static)  // MISSED?
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];  // OPT
		p(i,0) = COEFF_D0_P3*x[i]*(5.0*x2 - 3.0);
		p(i,1) = COEFF_D1_P3*(5.0*x2 - 1.0);
		p(i,2) = COEFF_D2_P3*x[i];
		p(i,3) = COEFF_D3_P3;  // MISSED?
	}
}


void LegendreP3::calculate_T(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	static constexpr int num_rows = l_ + 1;
	p.resize(num_rows, num_points);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for simd schedule(static)  // MISSED?
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];  // OPT
		p(0,i) = COEFF_D0_P3*x[i]*(5.0*x2 - 3.0);
		p(1,i) = COEFF_D1_P3*(5.0*x2 - 1.0);
		p(2,i) = COEFF_D2_P3*x[i];
		p(3,i) = COEFF_D3_P3;  // MISSED?
	}

	/*
	//x2_.resize(num_points);

	#pragma omp parallel  // NOT VEC
	{
		// TODO OPENMP CHUNKSIZE
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			const double x2 = x[i]*x[i];  // OPT
			p(0,i) = COEFF_D0_P3*x[i]*(5.0*x2 - 3.0);  // MISSED?
			p(1,i) = COEFF_D1_P3*(5.0*x2 - 1.0);  // OPT
			p(2,i) = COEFF_D2_P3*x[i];  // OPT
			p(3,i) = COEFF_D3_P3;  // OPT
		}
		// TODO OPENMP CHUNKSIZE
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			x2_[i] = x[i]*x[i];  // OPT
			p(0,i) = COEFF_D0_P3*x[i]*(5.0*x2_[i] - 3.0);  // MISSED?
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(1,i) = COEFF_D1_P3*(5.0*x2_[i] - 1.0);  // OPT
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(2,i) = COEFF_D2_P3*x[i];  // OPT
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(3,i) = COEFF_D3_P3;  // OPT
		}
	}
	*/
}



void LegendreP4::calculate(const double x, Vector<double>& p) const  // MISSED?
{
	p.resize(l_ + 1);

	double x2 = x*x;
	p[0] = COEFF_D0_P4*(x2*(35.0*x2 - 30.0) + 3.0);
	p[1] = COEFF_D1_P4*x*(7.0*x2 - 3.0);
	p[2] = COEFF_D2_P4*(7.0*x2 - 1.0);
	p[3] = COEFF_D3_P4*x;
	p[4] = COEFF_D4_P4;
}  // MISSED?


void LegendreP4::calculate(const Vector<double>& x, Matrix<double>& p) const  // MISSED?
{
	const int num_points = x.size();
	static constexpr int num_cols = l_ + 1;
	p.resize(num_points, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for simd schedule(static)  // MISSED?
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];  // OPT
		p(i,0) = COEFF_D0_P4*(x2*(35.0*x2 - 30.0) + 3.0);
		p(i,1) = COEFF_D1_P4*x[i]*(7.0*x2 - 3.0);
		p(i,2) = COEFF_D2_P4*(7.0*x2 - 1.0);
		p(i,3) = COEFF_D3_P4*x[i];
		p(i,4) = COEFF_D4_P4;  // MISSED?
	}
}  // MISSED?


void LegendreP4::calculate_T(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	static constexpr int num_rows = l_ + 1;
	p.resize(num_rows, num_points);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for simd schedule(static)  // MISSED?
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];  // OPT
		p(0,i) = COEFF_D0_P4*(x2*(35.0*x2 - 30.0) + 3.0);
		p(1,i) = COEFF_D1_P4*x[i]*(7.0*x2 - 3.0);
		p(2,i) = COEFF_D2_P4*(7.0*x2 - 1.0);
		p(3,i) = COEFF_D3_P4*x[i];
		p(4,i) = COEFF_D4_P4;  // MISSED?
	}

	/*
	//x2_.resize(num_points);
	#pragma omp parallel
	{
		// TODO OPENMP CHUNKSIZE
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			x2_[i] = x[i]*x[i];  // OPT
			p(0,i) = COEFF_D0_P4*(x2_[i]*(35.0*x2_[i] - 30.0) + 3.0);  // MISSED?
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(1,i) = COEFF_D1_P4*x[i]*(7.0*x2_[i] - 3.0);  // OPT
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(2,i) = COEFF_D2_P4*(7.0*x2_[i] - 1.0);  // OPT
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(3,i) = COEFF_D3_P4*x[i];  // OPT
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(4,i) = COEFF_D4_P4;  // OPT
		}
	}
	*/
}




void LegendreP6::calculate(const double x, Vector<double>& p) const
{
	p.resize(l_ + 1);

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
	static constexpr int num_cols = l_ + 1;
	p.resize(num_points, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for simd schedule(static)  // MISSED?
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];  // OPT
		p(i,0) = COEFF_D0_P6*(x2*(x2*(231.0*x2 - 315.0) + 105.0) - 5.0);
		p(i,1) = COEFF_D1_P6*x[i]*(x2*(33.0*x2 - 30.0) + 5.0);
		p(i,2) = COEFF_D2_P6*(x2*(33.0*x2 - 18.0) + 1.0);
		p(i,3) = COEFF_D3_P6*x[i]*(11.0*x2 - 3.0);
		p(i,4) = COEFF_D4_P6*(11.0*x2 - 1.0);
		p(i,5) = COEFF_D5_P6*x[i];
		p(i,6) = COEFF_D6_P6;  // MISSED?
	}
}


void LegendreP6::calculate_T(const Vector<double>& x, Matrix<double>& p) const
{
	const int num_points = x.size();
	static constexpr int num_rows = l_ + 1;
	p.resize(num_rows, num_points);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for simd schedule(static)  // MISSED?
	for ( int i=0; i<num_points; ++i ) {
		const double x2 = x[i]*x[i];  // OPT
		p(0,i) = COEFF_D0_P6*(x2*(x2*(231.0*x2 - 315.0) + 105.0) - 5.0);
		p(1,i) = COEFF_D1_P6*x[i]*(x2*(33.0*x2 - 30.0) + 5.0);
		p(2,i) = COEFF_D2_P6*(x2*(33.0*x2 - 18.0) + 1.0);
		p(3,i) = COEFF_D3_P6*x[i]*(11.0*x2 - 3.0);
		p(4,i) = COEFF_D4_P6*(11.0*x2 - 1.0);
		p(5,i) = COEFF_D5_P6*x[i];
		p(6,i) = COEFF_D6_P6;  // MISSED?
	}

	/*
	//x2_.resize(num_points);

	#pragma omp parallel
	{
		// TODO OPENMP CHUNKSIZE
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			const double x2 = x[i]*x[i];
			p(0,i) = COEFF_D0_P6*(x2*(x2*(231.0*x2 - 315.0) + 105.0) - 5.0);
			p(1,i) = COEFF_D1_P6*x[i]*(x2*(33.0*x2 - 30.0) + 5.0);
			p(2,i) = COEFF_D2_P6*(x2*(33.0*x2 - 18.0) + 1.0);
			p(3,i) = COEFF_D3_P6*x[i]*(11.0*x2 - 3.0);
			p(4,i) = COEFF_D4_P6*(11.0*x2 - 1.0);
			p(5,i) = COEFF_D5_P6*x[i];
			p(6,i) = COEFF_D6_P6;

			// ALT
			x2_[i] = x[i]*x[i];
			p(0,i) = COEFF_D0_P6*(x2_[i]*(x2_[i]*(231.0*x2_[i] - 315.0) + 105.0) - 5.0);
			p(1,i) = COEFF_D1_P6*x[i]*(x2_[i]*(33.0*x2_[i] - 30.0) + 5.0);
			p(2,i) = COEFF_D2_P6*(x2_[i]*(33.0*x2_[i] - 18.0) + 1.0);
			p(3,i) = COEFF_D3_P6*x[i]*(11.0*x2_[i] - 3.0);
			p(4,i) = COEFF_D4_P6*(11.0*x2_[i] - 1.0);
			p(5,i) = COEFF_D5_P6*x[i];
			p(6,i) = COEFF_D6_P6;
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(1,i) = COEFF_D1_P6*x[i]*(x2_[i]*(33.0*x2_[i] - 30.0) + 5.0);
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(2,i) = COEFF_D2_P6*(x2_[i]*(33.0*x2_[i] - 18.0) + 1.0);
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(3,i) = COEFF_D3_P6*x[i]*(11.0*x2_[i] - 3.0);
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(4,i) = COEFF_D4_P6*(11.0*x2_[i] - 1.0);
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(5,i) = COEFF_D5_P6*x[i];
		}
		#pragma omp for simd schedule(static)
		for ( int i=0; i<num_points; ++i ) {
			p(6,i) = COEFF_D6_P6;
		}
	}
	*/
}
