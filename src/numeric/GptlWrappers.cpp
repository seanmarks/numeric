// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#include "GptlWrappers.h"

// Parallelization
//#include "MpiCommunicator.h"
#include "OpenMP.h"

// GPTL headers (include here for best encapsulation)
#if HAVE_GPTL
#  include "gptl.h"
#  ifdef MPI_ENABLED 
#    include "gptlmpi.h"
#  endif
#endif


namespace GPTL {

//----- Loose Functions -----//

// This flag helps prevent multiple calls to GPTLinitialize
// - Please don't mess with it. Please.
// - TODO Move to a class for encapsulation?
static bool is_gptl_initialized_ = false;

void initialize()
{
#if HAVE_GPTL
	// If OpenMP is enabled, only the master thread should initialize GPTL
	// - Multiple initializations are an error
	// - It seems to be okay if MPI_Init has already been called, and multiple
	//   independent ranks call GPTLinitialize()
	//   - TODO: Verify with devs
	// - If the program is not in a parallel region when initialize() is called,
	//   the pragmas are ignored
	// - Note that 'master' does not impose an implicit barrier, so an explicit one is needed
	// - TODO Is it possible to make this threadsafe for use with other kinds of threads?
	//   - pthreads? C++ threads?
	#pragma omp master
	{
		if ( not is_gptl_initialized_ ) {
			GPTLinitialize();
			is_gptl_initialized_ = true;
		}
	}
	#pragma omp barrier
#endif // if HAVE_GPTL
}


void print(const std::string& file) 
{
#if HAVE_GPTL
	GPTLpr_file( file.c_str() );
#endif // if HAVE_GPTL
}


/*
void print(const std::string& file, const MpiCommunicator& mpi_communicator)
{
#if HAVE_GPTL
	int my_rank = mpi_communicator.getRank();
	std::string rank_file = file + "." + std::to_string(my_rank);
	GPTLpr_file( rank_file.c_str() );
#endif // if HAVE_GPTL
}
*/


/*
void printSummary(const std::string& file, const MpiCommunicator& mpi_communicator)
{
#if HAVE_GPTL
	GPTLpr_summary_file( mpi_communicator.getCommunicator(), file.c_str() );
#endif // if HAVE_GPTL
}
*/


void finalize()
{
#if HAVE_GPTL
#pragma omp master
	{
		if ( is_gptl_initialized_ ) {
			GPTLfinalize();
			is_gptl_initialized_ = false;
		}
	}
#pragma omp barrier
#endif // if HAVE_GPTL
}

//----- Timer -----//

Timer::Timer(const std::string& name):
	name_(name)
{
#if HAVE_GPTL
	GPTLinit_handle( name_.c_str(), &hash_index_ );
#endif
}


Timer::Timer():
	Timer("default")
{
}


void Timer::start()
{
#if HAVE_GPTL
	GPTLstart_handle( name_.c_str(), &hash_index_ );
	//GPTLstart( name_.c_str() );
#endif
}


void Timer::stop()
{
#if HAVE_GPTL
	GPTLstop_handle( name_.c_str(), &hash_index_ );
	//GPTLstop( name_.c_str() );
#endif // if HAVE_GPTL
}


//----- GlobalOptions -----//

GlobalOptions::GlobalOptions()
{	
	setDefaults();
}


void GlobalOptions::setDefaults()
{
#if HAVE_GPTL
	GPTLsetoption(GPTLpercent, 1);  // Print a column with % of first timer
	//GPTLsetoption(GPTLoverhead, 0);  // Don't print overhead estimate
#endif // if HAVE_GPTL
}

} // end namespace GPTL
