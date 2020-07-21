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

	l_ = l;
	if ( do_fast_harmonics_ ) {
		num_m_values_ = l_ + 1;
	}
	else {
		num_m_values_ = 2*l_ + 1;
	}

	// Use the appropriate Legendre specialization
	// - TODO: Register and set using a map?
	if ( l == 0 ) {
		legendre_ptr_ = std::unique_ptr<Legendre>( new LegendreP0 );
	}
	else if ( l == 3 ) {
		legendre_ptr_ = std::unique_ptr<Legendre>( new LegendreP3 );
	}
	else if ( l == 4 ) {
		legendre_ptr_ = std::unique_ptr<Legendre>( new LegendreP4 );
	}
	else if ( l == 6 ) {
		legendre_ptr_ = std::unique_ptr<Legendre>( new LegendreP6 );
	}
	else {
		std::stringstream err_ss;
		err_ss << "Error in " << FANCY_FUNCTION << "\n"
		       << "  Harmonic index l=" << l_ << " is not supported.\n";
		throw std::runtime_error( err_ss.str() );
	}

	// Coefficients
	coeff_Y_l_.resize(l_+1);
	const double prefac = (2.0*l_ + 1.0)/(4.0*constants::pi);
	for ( int m=0; m<=l_; ++m ) {
		coeff_Y_l_[m] = prefac * factorial(l-m)/static_cast<double>(factorial(l+m));
		if ( m > 0 ) {
			coeff_Y_l_[m] *= 2.0;  // extra factor not present with complex Y_{l,m}
		}
		coeff_Y_l_[m] = sqrt(coeff_Y_l_[m]);

		if ( m % 2 != 0 ) {
			coeff_Y_l_[m] = -coeff_Y_l_[m];  // Condon-Shortley phase
		}
	}
}


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
	legendre_ptr_->calculate(eta, p_x);

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
	for ( int m=1; m<=l_; ++m ) {  // TODO VEC CONTROL FLOW

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
	const int num_points = x.size();
	Y_l.resize(num_points, num_m_values_);
	if ( need_derivatives ) {
		derivs_Y_l.resize(num_points, num_m_values_);
	}


	//----- Precompute useful quantities -----//

	// Imaginary unit
	const Complex i = Complex(0.0, 1.0);

	inv_r_.resize(num_points);
	xhat_.resize(num_points);
	yhat_.resize(num_points);
	eta_.resize(num_points);
	zeta_.resize(num_points);
	zeta_m_minus_1_.resize(num_points, l_+1);
	if ( need_derivatives ) {
		deriv_eta_.resize(num_points);
		deriv_zeta_.resize(num_points);
	}

	#pragma omp parallel
	{
		// TODO: chunk size
		//int chunk_size = 
		#pragma omp for simd
		for ( int j=0; j<num_points; ++j ) {  // OPT
			// Floating-point multiplication is faster than floating-point division
			inv_r_[j] = 1.0/r[j];

			// Unit vector in direction of position vector x
			xhat_[j] = x[j][X_DIM]*inv_r_[j];
			yhat_[j] = x[j][Y_DIM]*inv_r_[j];
			eta_[j]  = x[j][Z_DIM]*inv_r_[j];  // eta = cos(theta) = z/r = zhat

			// zeta = (x + i*y)/r = xhat + i*yhat
			zeta_[j] = Complex(xhat_[j], yhat_[j]);  // OPT?
		}

		// FIXME VEC
		#pragma omp for
		for ( int j=0; j<num_points; ++j ) {  // OPT
			zeta_m_minus_1_(j,0) = Complex(1.0, 0.0);
			#pragma omp simd
			for ( int m=1; m<=l_; ++m ) {
				zeta_m_minus_1_(j,m) = zeta_m_minus_1_(j,m-1) * zeta_[j];
			}
		}

		// Derivatives (if requested) 
		if ( need_derivatives ) {
			#pragma omp for simd
			for ( int j=0; j<num_points; ++j ) {
				double zhat_times_inv_r = eta_[j]*inv_r_[j];

				deriv_eta_[j][X_DIM] = -xhat_[j]*zhat_times_inv_r;  // OPT
				deriv_eta_[j][Y_DIM] = -yhat_[j]*zhat_times_inv_r;
				deriv_eta_[j][Z_DIM] = (1.0 - eta_[j]*eta_[j])*inv_r_[j];

				deriv_zeta_[j][X_DIM] = (1.0 - xhat_[j]*zeta_[j])*inv_r_[j];   // MISSED
				deriv_zeta_[j][Y_DIM] = (i - yhat_[j]*zeta_[j])*inv_r_[j];
				deriv_zeta_[j][Z_DIM] = -zhat_times_inv_r*zeta_[j];
			}
		}
	}

	// Compute the Legendre polynomials and their derivatives
	// - TODO: is matrix version notably faster?
	p_x_.resize(num_points, num_m_values_);
	legendre_ptr_->calculate(eta_, p_x_);


	//----- Compute Y_{l,m} -----//

	#pragma omp parallel
	{
		#pragma omp for
		for ( int j=0; j<num_points; ++j ) {
			// m = 0
			int m = 0;
			Y_l(j,m) = coeff_Y_l_[m]*p_x_(j,m);
			if ( need_derivatives ) {
				for ( int d=0; d<N_DIM; ++d ) {
					derivs_Y_l(j,m)[d] = coeff_Y_l_[m]*p_x_(j,m+1)*deriv_eta_[j][d];
				}
			}

			// m > 0
			//Complex zeta_m_minus_1 = Complex(1.0, 0.0);
			for ( int m=1; m<=l_; ++m ) {  // TODO: VEC CONTROL FLOW
				Y_l(j,m) = coeff_Y_l_[m]*zeta_[j]*zeta_m_minus_1_(j,m-1)*p_x_(j,m);
				//zeta_m_minus_1 *= zeta_[j];  // TODO PRECOMPUTE
			}
			if ( need_derivatives ) {
				//Complex zeta_m_minus_1 = Complex(1.0, 0.0);
				for ( int m=1; m<=l_; ++m ) {  // TODO: VEC CONTROL FLOW
					for ( int d=0; d<N_DIM; ++d ) {
						derivs_Y_l(j,m)[d] = m*p_x_(j,m)*deriv_zeta_[j][d];
						if ( m < l_ ) {  // FIXME
							derivs_Y_l(j,m)[d] += deriv_eta_[j][d]*p_x_(j,m+1)*zeta_[j];
						}
						derivs_Y_l(j,m)[d] *= coeff_Y_l_[m]*zeta_m_minus_1_(j,m-1);
					}
					//zeta_m_minus_1 *= zeta_[j];
				}
			}
		}  // end loop over positions

		// m < 0
		// - Note: Y_{l,-m} = (-1)^m * complex_conjugate(Y_{l,m})
		if ( ! do_fast_harmonics_ ) {
			#pragma omp for
			for ( int j=0; j<num_points; ++j ) {
				for ( int m=1; m<=l_; m++ ) {  // TODO: VEC CONTROL FLOW
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

 
void SphericalHarmonics::calculateModifiedSphericalCoordinates(
		const Real3& x, double& r, double& phi, double& eta) const
{
	r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

	// Angular variables
	phi = atan2(x[1], x[0]);  // phi = arctan(y/x)
	if (phi < 0.0) { phi += 2.0*constants::pi; }
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
