// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "RealSphericalHarmonics.h"


RealSphericalHarmonics::RealSphericalHarmonics(const int l) 
{
	setHarmonicIndex(l);
}


void RealSphericalHarmonics::setHarmonicIndex(const int l) 
{
	FANCY_ASSERT( l >= 0, "invalid harmonic index: " << l );

	l_ = l;
	num_m_values_ = l_ + 1;

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
	coeff_y_l_.resize(num_m_values_);
	const double prefac = (2.0*l_ + 1.0)/(4.0*constants::pi);
	for ( int m=0; m<=l_; ++m ) {
		coeff_y_l_[m] = prefac * factorial(l-m)/static_cast<double>(factorial(l+m));
		if ( m > 0 ) {
			coeff_y_l_[m] *= 2.0;  // extra factor not present with complex Y_{l,m}
		}
		coeff_y_l_[m] = sqrt(coeff_y_l_[m]);

		if ( m % 2 != 0 ) {
			coeff_y_l_[m] = -coeff_y_l_[m];  // Condon-Shortley phase
		}
	}
}


void RealSphericalHarmonics::calculate(
		const VectorReal3& x, const Vector<double>& r, const bool need_derivatives,
		Matrix<double>& y_l, Matrix<Real3>& derivs_y_l
) const
{
	FANCY_ASSERT( legendre_ptr_ != nullptr, "missing object");
	FANCY_ASSERT( legendre_ptr_->getHarmonicIndex() == l_, "harmonic index mismatch" );

	// Ensure enough memory has been allocated
	const int num_points = x.size();
	y_l.resize(num_points, 2*l_+1);
	if ( need_derivatives ) {
		derivs_y_l.resize(num_points, 2*l_+1);
	}

	inv_r_.resize(num_points);
	xhat_.resize(num_points);
	yhat_.resize(num_points);
	eta_.resize(num_points);
	if ( need_derivatives ) {
		for ( int d=0; d<N_DIM; ++d ) {
			re_deriv_zeta_[d].resize(num_points);
			im_deriv_zeta_[d].resize(num_points);
			deriv_eta_[d].resize(num_points);
		}
		/*
		re_deriv_zeta_.resize(num_points);
		im_deriv_zeta_.resize(num_points);
		deriv_eta_.resize(num_points);
		*/
	}

	#pragma omp parallel
	{
		// TODO: CHUNK_SIZE
		const int chunk_size = 8;
		#pragma omp for simd schedule(static,chunk_size)
		for ( int j=0; j<num_points; ++j ) {  // OPT
			// Floating-point multiplication is faster than floating-point division
			inv_r_[j] = 1.0/r[j];

			// Unit vector in direction of position vector x
			xhat_[j] = x[j][X_DIM]*inv_r_[j];
			yhat_[j] = x[j][Y_DIM]*inv_r_[j];
			eta_[j]  = x[j][Z_DIM]*inv_r_[j];  // eta = cos(theta) = z/r = zhat
		}

		if ( need_derivatives ) {
			const int chunk_size = 8;
			#pragma omp for simd schedule(static,chunk_size)
			for ( int j=0; j<num_points; ++j ) {
				double zhat_times_inv_r = eta_[j]*inv_r_[j];

				deriv_eta_[X_DIM][j] = -xhat_[j]*zhat_times_inv_r;
				deriv_eta_[Y_DIM][j] = -yhat_[j]*zhat_times_inv_r;
				deriv_eta_[Z_DIM][j] = (1.0 - eta_[j]*eta_[j])*inv_r_[j];

				re_deriv_zeta_[X_DIM][j] = inv_r_[j]*(1.0 - xhat_[j]*xhat_[j]);
				re_deriv_zeta_[Y_DIM][j] = -inv_r_[j]*xhat_[j]*yhat_[j];
				re_deriv_zeta_[Z_DIM][j] = deriv_eta_[X_DIM][j];

				im_deriv_zeta_[X_DIM][j] = re_deriv_zeta_[Y_DIM][j];
				im_deriv_zeta_[Y_DIM][j] = inv_r_[j]*(1.0 - yhat_[j]*yhat_[j]);
				im_deriv_zeta_[Z_DIM][j] = deriv_eta_[Y_DIM][j];

				/*
				// FIXME: alternative layout

				deriv_eta_[j][X_DIM] = -xhat_[j]*zhat_times_inv_r;
				deriv_eta_[j][Y_DIM] = -yhat_[j]*zhat_times_inv_r;
				deriv_eta_[j][Z_DIM] = (1.0 - eta_[j]*eta_[j])*inv_r_[j];

				re_deriv_zeta_[j][X_DIM] = inv_r_[j]*(1.0 - xhat_[j]*xhat_[j]);
				re_deriv_zeta_[j][Y_DIM] = -inv_r_[j]*xhat_[j]*yhat_[j];
				re_deriv_zeta_[j][Z_DIM] = deriv_eta_[j][X_DIM];

				im_deriv_zeta_[j][X_DIM] = re_deriv_zeta_[j][Y_DIM];
				im_deriv_zeta_[j][Y_DIM] = inv_r_[j]*(1.0 - yhat_[j]*yhat_[j]);
				im_deriv_zeta_[j][Z_DIM] = deriv_eta_[j][Y_DIM];
				*/
			}
		}
	}

	// Compute the Legendre polynomials and their derivatives
	// - TODO: is matrix version notably faster?
	p_x_.resize(num_points, num_m_values_);
	legendre_ptr_->calculate(eta_, p_x_);

	#pragma omp parallel
	{
		#pragma omp for
		for ( int j=0; j<num_points; ++j ) {
			// m = 0
			int m = 0;
			y_l(j,m) = coeff_y_l_[m]*p_x_(j,m);
			if ( need_derivatives ) {
				#pragma omp simd
				for ( int d=0; d<N_DIM; ++d ) {
					derivs_y_l(j,m)[d] = coeff_y_l_[m]*p_x_(j,m+1)*deriv_eta_[d][j];
					//derivs_y_l(j,m)[d] = coeff_y_l_[m]*p_x_(j,m+1)*deriv_eta_[j][d];
				}
			}

			// FIXME
			// m > 0
			//Complex zeta_m_minus_1 = Complex(1.0, 0.0);
			double re_zeta_m_minus_1 = 1.0, im_zeta_m_minus_1 = 0.0;
			for ( int m=1; m<=l_; ++m ) {
				y_l(j,m) = coeff_y_l_[m]*xhat_[j]*xhat_[j]*p_x_(j,m);  // FIXME
				//y_l(j,m) = coeff_y_l_[m]*zeta_[j]*zeta_m_minus_1_(j,m-1)*p_x_(j,m);
			}
			if ( need_derivatives ) {
				for ( int m=1; m<=l_; ++m ) {
					for ( int d=0; d<N_DIM; ++d ) {
						derivs_y_l(j,m)[d] = m*p_x_(j,m)*re_deriv_zeta_[d][j];  // FIXME
						//derivs_y_l(j,m)[d] = m*p_x_(j,m)*deriv_zeta_[j][d];
						if ( m < l_ ) {  // FIXME
							derivs_y_l(j,m)[d] += deriv_eta_[d][j]*p_x_(j,m+1)*xhat_[j];  // FIXME
							//derivs_y_l(j,m)[d] += deriv_eta_[j][d]*p_x_(j,m+1)*zeta_[j];  
						}
						derivs_y_l(j,m)[d] *= coeff_y_l_[m]*xhat_[j];
						//derivs_y_l(j,m)[d] *= coeff_y_l_[m]*zeta_m_minus_1_(j,m-1);
					}
					//zeta_m_minus_1 *= zeta_[j];
				}
			}
		}
	}

/*
	//----- Compute Y_{l,m} -----//

	#pragma omp parallel
	{
		#pragma omp for
		for ( int j=0; j<num_points; ++j ) {
			// m = 0
			int m = 0;
			y_l(j,m) = coeff_y_l_[m]*p_x_(j,m);
			if ( need_derivatives ) {
				for ( int d=0; d<N_DIM; ++d ) {
					derivs_y_l(j,m)[d] = coeff_y_l_[m]*p_x_(j,m+1)*deriv_eta_[j][d];
				}
			}

			// m > 0
			//Complex zeta_m_minus_1 = Complex(1.0, 0.0);
			for ( int m=1; m<=l_; ++m ) {  // TODO: VEC CONTROL FLOW
				y_l(j,m) = coeff_y_l_[m]*zeta_[j]*zeta_m_minus_1_(j,m-1)*p_x_(j,m);
				//zeta_m_minus_1 *= zeta_[j];  // TODO PRECOMPUTE
			}
			if ( need_derivatives ) {
				//Complex zeta_m_minus_1 = Complex(1.0, 0.0);
				for ( int m=1; m<=l_; ++m ) {  // TODO: VEC CONTROL FLOW
					for ( int d=0; d<N_DIM; ++d ) {
						derivs_y_l(j,m)[d] = m*p_x_(j,m)*deriv_zeta_[j][d];
						if ( m < l_ ) {  // FIXME
							derivs_y_l(j,m)[d] += deriv_eta_[j][d]*p_x_(j,m+1)*zeta_[j];
						}
						derivs_y_l(j,m)[d] *= coeff_y_l_[m]*zeta_m_minus_1_(j,m-1);
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
						y_l(j,index) = std::conj( y_l(j,m) );
						if ( need_derivatives ) {
							for ( int d=0; d<N_DIM; ++d ) {
								derivs_y_l(j,index)[d] = std::conj( derivs_y_l(j,m)[d] );
							}
						}
					}
					else {
						y_l(j,index) = std::conj( -y_l(j,m) );
						if ( need_derivatives ) {
							for ( int d=0; d<N_DIM; ++d ) {
								derivs_y_l(j,index)[d] = std::conj( -derivs_y_l(j,m)[d] );
							}
						}
					}
				}
			}
		}
	} // end pragma omp parallel
*/
}


/* 
double RealSphericalHarmonics::norm(const VectorComplex& vec) const
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

 
RealSphericalHarmonics::Complex RealSphericalHarmonics::innerProduct(
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

 
RealSphericalHarmonics::Complex3 RealSphericalHarmonics::innerProduct(
		const VectorComplex& x, const VectorComplex3& dy) const
{
	const auto dy_first = dy.begin();
	const auto dy_last  = std::next(dy.begin(), num_m_values_);

	return this->innerProduct(x, dy_first, dy_last);
}

 
RealSphericalHarmonics::Complex3 RealSphericalHarmonics::innerProduct(
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
*/
