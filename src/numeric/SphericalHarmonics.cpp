// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "SphericalHarmonics.h"


SphericalHarmonics::SphericalHarmonics(const int l, const bool do_fast_harmonics) 
 : do_fast_harmonics_(do_fast_harmonics)
{
	setHarmonicIndex(l);
}


void SphericalHarmonics::setHarmonicIndex(const int l) 
{
	FANCY_ASSERT( l >= 0, "invalid harmonic index: " << l );

	// Set harmonic index and number of m-values
	l_ = l;
	if ( do_fast_harmonics_ ) {
		num_m_values_ = l_ + 1;
	}
	else {
		num_m_values_ = 2*l_ + 1;
	}

	// Set coefficients and necessary legendreP function
	// - TODO: Register and set using a map?
	if ( l == 0 ) {
		coeff_Y_l_ = coeff_Y_0_;
		legendreP_function_        = [this](const double x, Vector<double>& p)          { legendreP0(x, p); };
		legendreP_matrix_function_ = [this](const Vector<double>& x, Matrix<double>& p) { legendreP0(x, p); };
	}
	else if ( l == 3 ) {
		coeff_Y_l_ = coeff_Y_3_;
		legendreP_function_        = [this](const double x, Vector<double>& p)          { legendreP3(x, p); };
		legendreP_matrix_function_ = [this](const Vector<double>& x, Matrix<double>& p) { legendreP3(x, p); };
	}
	else if ( l == 4 ) {
		coeff_Y_l_ = coeff_Y_4_;
		legendreP_function_        = [this](const double x, Vector<double>& p)          { legendreP4(x, p); };
		legendreP_matrix_function_ = [this](const Vector<double>& x, Matrix<double>& p) { legendreP4(x, p); };
	}
	else if ( l == 6 ) {
		coeff_Y_l_ = coeff_Y_6_;
		legendreP_function_        = [this](const double x, Vector<double>& p)          { legendreP6(x, p); };
		legendreP_matrix_function_ = [this](const Vector<double>& x, Matrix<double>& p) { legendreP6(x, p); };
	}
	else {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  Harmonic index l=" << l_ << " is not supported.\n";
		throw std::runtime_error( err_ss.str() );
	}
}


// Calculate Y_l,m for m = -l, ... , l
// - Output format: Y_l: [m=0, 1, ..., l, -1, -2, ... -l] (similarly for derivatives)
void SphericalHarmonics::calculate_Y_l(
		const Real3& x, const double r, const bool need_derivatives,
		VectorComplex& Y_l, VectorComplex3& derivs_Y_l
) const
{
	// Input checks
	FANCY_ASSERT( r > 0.0, "invalid input: r = " << r );
	
	// Ensure enough memory has been allocated
	Y_l.resize(num_m_values_);
	if ( need_derivatives ) { derivs_Y_l.resize(num_m_values_); }


	//----- Precompute useful quantities -----//

	// Imaginary unit
	const Complex i = Complex(0.0, 1.0);

	// Floating-point multiplication is faster than floating-point division
	double inv_r = 1.0/r;

	// Unit vector in direction of position vector x
	Real3 xhat;
	for ( int d=0; d<N_DIM; ++d ) {
		xhat[d] = x[d]*inv_r;
	}

	// zeta = (x + i*y)/r = xhat + i*yhat
	Complex zeta = Complex(xhat[X_DIM], xhat[Y_DIM]);

	// eta = cos(theta) = z/r
	double eta = xhat[Z_DIM];

	// Derivatives (if requested) 
	Real3    deriv_eta;
	Complex3 deriv_zeta;
	double   zhat_times_inv_r = xhat[Z_DIM]*inv_r;
	if ( need_derivatives ) {
		deriv_eta[X_DIM] = -xhat[X_DIM]*zhat_times_inv_r;
		deriv_eta[Y_DIM] = -xhat[Y_DIM]*zhat_times_inv_r;
		deriv_eta[Z_DIM] = (1.0 - eta*eta)*inv_r;

		deriv_zeta[X_DIM] = (1.0 - xhat[X_DIM]*zeta)*inv_r;
		deriv_zeta[Y_DIM] = (i - xhat[Y_DIM]*zeta)*inv_r;
		deriv_zeta[Z_DIM] = -zhat_times_inv_r*zeta;
	}


	//----- Compute Y_{l,m} -----//

	// First, compute the Legendre polynomials and their derivatives
	Vector<double> p_x(num_m_values_);
	legendreP_function_(eta, p_x);

	// m = 0
	int m = 0;
	Y_l[m] = coeff_Y_l_[m]*p_x[m];
	if ( need_derivatives ) {
		if ( l_ > 0 ) {
			for ( int d=0; d<N_DIM; ++d ) {
				derivs_Y_l[m][d] = coeff_Y_l_[m]*p_x[m+1]*deriv_eta[d];
			}
		}
		else {
			// Y_{0,0} is constant
			for ( int d=0; d<N_DIM; ++d ) {
				derivs_Y_l[m][d] = 0.0;
			}
		}
	}

	// m > 0
	Complex zeta_m_minus_1 = Complex(1.0, 0.0);
	for ( int m=1; m<=l_; ++m ) {

		Y_l[m] = coeff_Y_l_[m]*zeta*zeta_m_minus_1*p_x[m];

		if ( need_derivatives ) {
			for ( int d=0; d<N_DIM; ++d ) {
				derivs_Y_l[m][d] = (m*p_x[m])*deriv_zeta[d];
				if ( m < l_ ) {
					derivs_Y_l[m][d] += (deriv_eta[d]*p_x[m+1])*zeta;
				}
				derivs_Y_l[m][d] *= coeff_Y_l_[m]*zeta_m_minus_1;
			}
		}

		zeta_m_minus_1 *= zeta;
	}

	// m < 0
	// - Note: Y_{l,-m} = (-1)^m * complex_conjugate(Y_{l,m})
	if ( not do_fast_harmonics_ ) {
		int index;
		for ( int m=1; m<=l_; m++ ) {
			index = l_ + m;
			if ( m % 2 == 0 ) {
				Y_l[index] = std::conj( Y_l[m] );
				if ( need_derivatives ) {
					for ( int d=0; d<N_DIM; ++d ) {
						derivs_Y_l[index][d] = std::conj( derivs_Y_l[m][d] );
					}
				}
			}
			else {
				Y_l[index] = std::conj( -Y_l[m] );
				if ( need_derivatives ) {
					for ( int d=0; d<N_DIM; ++d ) {
						derivs_Y_l[index][d] = std::conj( -derivs_Y_l[m][d] );
					}
				}
			}
		}
	}
}


void SphericalHarmonics::calculate(
	const VectorReal3& x, const Vector<double>& r, const bool need_derivatives,
	Matrix<Complex>& Y_l, Matrix<Complex3>& derivs_Y_l
) const
{
	// Input checks
	//FANCY_ASSERT( r > 0.0, "invalid input: r = " << r );
	
	// Ensure enough memory has been allocated
	const int len = x.size();
	Y_l.resize(len, num_m_values_);
	if ( need_derivatives ) { derivs_Y_l.resize(len, num_m_values_); }


	//----- Precompute useful quantities -----//

	// Imaginary unit
	const Complex i = Complex(0.0, 1.0);

	#pragma omp parallel
	{
		// TODO: chunk size
		//int chunk_size = 
		#pragma omp for 
		for ( int j=0; j<len; ++j ) {
			// Floating-point multiplication is faster than floating-point division
			double inv_r = 1.0/r[j];

			// Unit vector in direction of position vector x
			Real3 xhat;
			for ( int d=0; d<N_DIM; ++d ) {
				xhat[d] = x[j][d]*inv_r;
			}

			// zeta = (x + i*y)/r = xhat + i*yhat
			Complex zeta = Complex(xhat[X_DIM], xhat[Y_DIM]);

			// eta = cos(theta) = z/r
			double eta = xhat[Z_DIM];

			// Derivatives (if requested) 
			Real3    deriv_eta;
			Complex3 deriv_zeta;
			double   zhat_times_inv_r = xhat[Z_DIM]*inv_r;
			if ( need_derivatives ) {
				deriv_eta[X_DIM] = -xhat[X_DIM]*zhat_times_inv_r;
				deriv_eta[Y_DIM] = -xhat[Y_DIM]*zhat_times_inv_r;
				deriv_eta[Z_DIM] = (1.0 - eta*eta)*inv_r;

				deriv_zeta[X_DIM] = (1.0 - xhat[X_DIM]*zeta)*inv_r;
				deriv_zeta[Y_DIM] = (i - xhat[Y_DIM]*zeta)*inv_r;
				deriv_zeta[Z_DIM] = -zhat_times_inv_r*zeta;
			}


			//----- Compute Y_{l,m} -----//

			// First, compute the Legendre polynomials and their derivatives
			// - TODO: is matrix version notably faster?
			Vector<double> p_x(num_m_values_);
			legendreP_function_(eta, p_x);

			// m = 0
			int m = 0;
			Y_l(j,m) = coeff_Y_l_[m]*p_x[m];
			if ( need_derivatives ) {
				if ( l_ > 0 ) {
					for ( int d=0; d<N_DIM; ++d ) {
						derivs_Y_l(j,m)[d] = coeff_Y_l_[m]*p_x[m+1]*deriv_eta[d];
					}
				}
				else {
					// Y_{0,0} is constant
					for ( int d=0; d<N_DIM; ++d ) {
						derivs_Y_l(j,m)[d] = 0.0;
					}
				}
			}

			// m > 0
			Complex zeta_m_minus_1 = Complex(1.0, 0.0);
			for ( int m=1; m<=l_; ++m ) {
				Y_l(j,m) = coeff_Y_l_[m]*zeta*zeta_m_minus_1*p_x[m];

				if ( need_derivatives ) {
					for ( int d=0; d<N_DIM; ++d ) {
						derivs_Y_l(j,m)[d] = (m*p_x[m])*deriv_zeta[d];
						if ( m < l_ ) {
							derivs_Y_l(j,m)[d] += (deriv_eta[d]*p_x[m+1])*zeta;
						}
						derivs_Y_l(j,m)[d] *= coeff_Y_l_[m]*zeta_m_minus_1;
					}
				}

				zeta_m_minus_1 *= zeta;
			}
		}  // end loop over positions

		// m < 0
		// - Note: Y_{l,-m} = (-1)^m * complex_conjugate(Y_{l,m})
		if ( ! do_fast_harmonics_ ) {
			#pragma omp for
			for ( int j=0; j<len; ++j ) {
				for ( int m=1; m<=l_; m++ ) {
					int index = l_ + m;
					if ( m % 2 == 0 ) {
						Y_l(j,index) = std::conj( Y_l(j,m) );
						if ( need_derivatives ) {
							for ( int d=0; d<N_DIM; ++d ) {
								derivs_Y_l(j,index)[d] = std::conj( derivs_Y_l(j,m)[d] );
							}
						}
					}
					else {
						Y_l(j,index) = std::conj( -Y_l(j,m) );
						if ( need_derivatives ) {
							for ( int d=0; d<N_DIM; ++d ) {
								derivs_Y_l(j,index)[d] = std::conj( -derivs_Y_l(j,m)[d] );
							}
						}
					}
				}
			}
		}
	} // end pragma omp parallel
}


void SphericalHarmonics::legendreP0(const double x, Vector<double>& p) const
{
	p.assign(1, coeff_Y_0_[0]);
}


void SphericalHarmonics::legendreP3(const double x, Vector<double>& p) const
{
	p.resize(4);

	double x2 = x*x;
	p[0] = COEFF_P3*x*(5.0*x2 - 3.0);
	p[1] = COEFF_D1_P3*(5.0*x2 - 1.0);
	p[2] = COEFF_D2_P3*x;
	p[3] = COEFF_D3_P3;
}


void SphericalHarmonics::legendreP4(const double x, Vector<double>& p) const
{
	p.resize(5);

	double x2 = x*x;
	p[0] = COEFF_P4*(x2*(35.0*x2 - 30.0) + 3.0);
	p[1] = COEFF_D1_P4*x*(7.0*x2 - 3.0);
	p[2] = COEFF_D2_P4*(7.0*x2 - 1.0);
	p[3] = COEFF_D3_P4*x;
	p[4] = COEFF_D4_P4;
}


void SphericalHarmonics::legendreP6(const double x, Vector<double>& p) const
{
	p.resize(7);

	double x2 = x*x;
	double x2_11 = 11.0*x2;
	double x2_33 = 33.0*x2;
	p[0] = COEFF_P6*(x2*(x2*(231.0*x2 - 315.0) + 105.0) - 5.0);
	p[1] = COEFF_D1_P6*x*(x2*(x2_33 - 30.0) + 5.0);
	p[2] = COEFF_D2_P6*(x2*(x2_33 - 18.0) + 1.0);
	p[3] = COEFF_D3_P6*x*(x2_11 - 3.0);
	p[4] = COEFF_D4_P6*(x2_11 - 1.0);
	p[5] = COEFF_D5_P6*x;
	p[6] = COEFF_D6_P6;
}


void SphericalHarmonics::legendreP0(const Vector<double>& x, Matrix<double>& p) const
{
	const int len = x.size();
	p.assign( {{len, 1}}, coeff_Y_0_[0] );
}


void SphericalHarmonics::legendreP3(const Vector<double>& x, Matrix<double>& p) const
{
	const int len = x.size();
	static constexpr int num_cols = 4;
	p.resize(len, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for
	for ( int i=0; i<len; ++i ) {
		const double x2 = x[i]*x[i];
		p(i,0) = COEFF_P3*x[i]*(5.0*x2 - 3.0);
		p(i,1) = COEFF_D1_P3*(5.0*x2 - 1.0);
		p(i,2) = COEFF_D2_P3*x[i];
		p(i,3) = COEFF_D3_P3;
	}
}


void SphericalHarmonics::legendreP4(const Vector<double>& x, Matrix<double>& p) const
{
	const int len = x.size();
	static constexpr int num_cols = 5;
	p.resize(len, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for
	for ( int i=0; i<len; ++i ) {
		const double x2 = x[i]*x[i];
		p(i,0) = COEFF_P4*(x2*(35.0*x2 - 30.0) + 3.0);
		p(i,1) = COEFF_D1_P4*x[i]*(7.0*x2 - 3.0);
		p(i,2) = COEFF_D2_P4*(7.0*x2 - 1.0);
		p(i,3) = COEFF_D3_P4*x[i];
		p(i,4) = COEFF_D4_P4;
	}
}


void SphericalHarmonics::legendreP6(const Vector<double>& x, Matrix<double>& p) const
{
	const int len = x.size();
	static constexpr int num_cols = 7;
	p.resize(len, num_cols);

	// TODO OPENMP CHUNKSIZE
	#pragma omp parallel for
	for ( int i=0; i<len; ++i ) {
		const double x2 = x[i]*x[i];
		const double x2_11 = 11.0*x2;
		const double x2_33 = 33.0*x2;
		p(i,0) = COEFF_P6*(x2*(x2*(231.0*x2 - 315.0) + 105.0) - 5.0);
		p(i,1) = COEFF_D1_P6*x[i]*(x2*(x2_33 - 30.0) + 5.0);
		p(i,2) = COEFF_D2_P6*(x2*(x2_33 - 18.0) + 1.0);
		p(i,3) = COEFF_D3_P6*x[i]*(x2_11 - 3.0);
		p(i,4) = COEFF_D4_P6*(x2_11 - 1.0);
		p(i,5) = COEFF_D5_P6*x[i];
		p(i,6) = COEFF_D6_P6;
	}
}

 
void SphericalHarmonics::calculateModifiedSphericalCoordinates(
		const Real3& x, double& r, double& phi, double& eta) const
{
	r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

	// Angular variables
	phi = atan2(x[1], x[0]);  // phi = arctan(y/x)
	if (phi < 0.0) { phi += 2.0*PI_; }
	eta = x[2]/r;       		  // eta = cos(theta) = z/r

	// Make sure roundoff error doesn't make eta leave the range [-1,1]
	if      ( eta >  1.0 ) { eta =  1.0; }
	else if ( eta < -1.0 ) { eta = -1.0; }
}

 
double SphericalHarmonics::norm(const VectorComplex& vec) const
{
	double norm = 0.0;
	if ( do_fast_harmonics_ ) {
		// m>0 terms
		for ( int m=1; m<=l_; ++m ) {
			norm += (vec[m].real())*(vec[m].real()) + (vec[m].imag())*(vec[m].imag());
		}

		// Add term for m=0, and double sum of m>0 terms to account for m<0
		norm = (vec[0].real())*(vec[0].real()) + (vec[0].imag())*(vec[0].imag()) + 2.0*norm;
	}
	else {
		// Sum over all 2*l+1 values
		for ( int m=0; m<=num_m_values_; ++m ) {
			norm += (vec[m].real())*(vec[m].real()) + (vec[m].imag())*(vec[m].imag());
		}
	}

	norm = sqrt(norm);
	return norm;
}

 
SphericalHarmonics::Complex SphericalHarmonics::innerProduct(
		const VectorComplex& x, const VectorComplex& y) const
{
	Complex product(0.0, 0.0);
	if ( do_fast_harmonics_ ) {
		// m>0 terms
		for ( int m=1; m<=l_; ++m ) {
			product += x[m] * std::conj(y[m]);
		}

		// Add term for m=0, and double sum of m>0 terms to account for m<0
		product = x[0] * std::conj(y[0]) + 2.0*product;
	}
	else {
		// Sum over all 2*l+1 values
		for ( int m=0; m<=num_m_values_; ++m ) {
			product += x[m] * std::conj(y[m]);
		}
	}

	return product;
}

 
SphericalHarmonics::Complex3 SphericalHarmonics::innerProduct(
		const VectorComplex& x, const VectorComplex3& dy) const
{
	const auto dy_first = dy.begin();
	const auto dy_last  = std::next(dy.begin(), num_m_values_);

	return this->innerProduct(x, dy_first, dy_last);
}

 
SphericalHarmonics::Complex3 SphericalHarmonics::innerProduct(
		const VectorComplex& x,
		const VectorComplex3::const_iterator& dy_first, 
		const VectorComplex3::const_iterator& dy_last
) const
{
	Complex3 product;
	if ( do_fast_harmonics_ ) {
		for ( int d=0; d<N_DIM; ++d ) {
			product[d] = Complex(0.0, 0.0);

			// m>0
			auto dy_it = std::next(dy_first, 1);
			for ( int m=1; m<=l_; ++m ) {
				product[d] += x[m] * std::conj((*dy_it)[d]);
				++dy_it;
			}

			// Add term for m=0, and double terms for m>0 to account for m<0
			product[d] = x[0] * std::conj((*dy_first)[d]) + 2.0*product[d];
		}
	}
	else {
		for ( int d=0; d<N_DIM; ++d ) {
			product[d] = Complex(0.0, 0.0);

			// Sum over all terms
			auto dy_it = dy_first;
			for ( int m=0; m<=num_m_values_; ++m ) {
				product[d] += x[m] * std::conj((*dy_it)[d]);
				++dy_it;
			}
		}
	}

	return product;
}
