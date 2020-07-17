// Wrappers around GPTL library calls
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// TODO:
// - Throw an exception with some functions (e.g. print_summary()) if
//   GPTL is not available/initialized?

#pragma once
#ifndef GPTL_WRAPPER_H
#define GPTL_WRAPPER_H

#include <iostream>
#include <string>

#include "OpenMP.h"

class MpiCommunicator;


namespace GPTL {

void initialize();

void print(const std::string& file);

//void print(const std::string& file, const MpiCommunicator& mpi_communicator);

//void printSummary(const std::string& file, const MpiCommunicator& mpi_communicator);

void finalize();

// Wrapper around GPTL timer functions
// - TODO: use a registrar for timers so that names are unique?
// - TODO: do anything with error codes?
class Timer
{
 public:
	Timer(const std::string& name);

	Timer();

	void start();
	void stop();

 private:
	std::string name_ = "default";

	// Note: passing '0' to GPTLinit_handle() makes it produce the hash index corresponding to
	// the given timer name
	int hash_index_ = 0;
};

// Wrapper around 'GPTLsetoption'
// - The default constructor sets default option values
class GlobalOptions
{
 public:
	GlobalOptions();

	void setDefaults();
};

// TODO Find a way to have a 'static'-like (shared) object, since GPTL options
//      are all shared globally
//static GlobalOptions global_options("test");

} // end namespace GPTL

#endif // ifndef GPTL_WRAPPER_H
