
#pragma once
#ifndef MAIN_H
#define MAIN_H

#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

#include "numeric/Assert.h"
#include "numeric/Aligned.h"
#include "numeric/CompareNumbers.h"
#include "numeric/Timer.h"

void comparePerformance(std::string header, const double rmsd, const Timer& timer_old, const Timer& timer_new)
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


// Cache-aligned vector
template<typename T>
using AlignedAlloc = numeric::aligned::CacheAlignedAllocator<T>;


template<typename V, typename W>
double rmsd(const V& std_vec, const W& new_vec)
{
	double sum = 0.0;
	unsigned length = std_vec.size();
	assert( length == new_vec.size() );

	for ( unsigned i=0; i<length; ++i ) {
		double delta = std_vec[i] - new_vec[i];
		sum += delta*delta;
	}
	return sqrt( sum/static_cast<double>(length) );
}

#endif // ifndef MAIN_H
