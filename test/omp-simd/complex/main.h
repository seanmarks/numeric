
#pragma once
#ifndef MAIN_H
#define MAIN_H

#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

#include "numeric/Aligned.h"
#include "numeric/Assert.h"
#include "numeric/CompareNumbers.h"
#include "numeric/ComplexVector.h"
#include "numeric/Timer.h"

void comparePerformance(std::string header, const double rmsd, const Timer& timer_new, const Timer& timer_old)
{
	double t_old_ms = timer_old.get_duration_ms();
	double t_new_ms = timer_new.get_duration_ms();

	double t_old_us = timer_old.get_duration_us();
	double t_new_us = timer_new.get_duration_us();
	double speedup = 100.0*(t_old_us/t_new_us - 1.0);

	if ( header.empty() ) {
		header = "Compare";
	}
	else if ( header.back() == '\n' ) {
		header.pop_back();
	}

	std::cout << header << "\n"
	          << "  (rmsd = " << rmsd << ")\n"
	          << "  old: " << t_old_ms << " ms (" << t_old_us << " us)\n"
	          << "  new: " << t_new_ms << " ms (" << t_new_us << " us)\n"
	          << "  speedup: " << speedup << "%\n"
	          << std::endl;
}


// TODO update
template<typename T>
double rmsd(
	const numeric::aligned::ComplexVector<T>& new_vec,
	const std::vector<std::complex<T>>&       std_vec
)
{
	double sum = 0.0;
	unsigned length = std_vec.size();
	assert( length == new_vec.size() );

	for ( unsigned i=0; i<length; ++i ) {
		sum += std::norm( std_vec[i] - new_vec[i] );  // std::norm(z) = |z|^2
	}
	return sqrt( sum/static_cast<double>(length) );
}


template<typename T>
void check(
	const numeric::aligned::ComplexVector<T>& new_vec,
	const std::vector<std::complex<T>>&       std_vec
)
{
	unsigned length = std_vec.size();
	FANCY_ASSERT( length == new_vec.size(), "length mismatch" );

	numeric::AlmostEqualUlps<T> almost_equal_ulps;
	almost_equal_ulps.setMaxUlpsDiff(10);  // TODO: tune this

	for ( unsigned i=0; i<length; ++i ) {
		FANCY_ASSERT( almost_equal_ulps(new_vec[i].real(), std_vec[i].real()),
		              "numbers should be close in ULPs"
		              << "(value = " << new_vec[i].real() << ", ref = " << std_vec[i].real()
		              << ", diff = " << new_vec[i].real() - std_vec[i].real() << ")" );

		FANCY_ASSERT( almost_equal_ulps(new_vec[i].imag(), std_vec[i].imag()),
		              "numbers should be close in ULPs"
		              << "(value = " << new_vec[i].imag() << ", ref = " << std_vec[i].imag()
		              << ", diff = " << new_vec[i].imag() - std_vec[i].imag() << ")" );
	}
}

#endif // ifndef MAIN_H
