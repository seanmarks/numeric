// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#pragma once
#ifndef MAIN_H
#define MAIN_H

#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

#include "numeric/Aligned.h"
#include "numeric/Assert.h"
#include "numeric/CommonTypes.h"
#include "numeric/CompareNumbers.h"
#include "numeric/ComplexVector.h"
#include "numeric/Constants.h"
#include "numeric/GptlWrappers.h"
#include "numeric/Matrix.h"
#include "numeric/Legendre.h"
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

#endif // ifndef MAIN_H
