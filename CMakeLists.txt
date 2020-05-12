cmake_minimum_required(VERSION 3.15)

# Include directory with custom CMake modules
# - These contain useful helper functions
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(setIfUnset)
include(WriteCXXSource)

project(numeric)

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Default to 'Release' build type
setIfUnset(CMAKE_BUILD_TYPE Release)

# Verbose output
# - TODO: restrict to 'debug'?
set(CMAKE_VERBOSE_MAKEFILE ON)

# Put all binaries and libraries into subdirectories, instead of the main 'build' directory
# - Adapted from CMake tutorial (v3.15)
setIfUnset(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
setIfUnset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
setIfUnset(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")


#######################################
### Important variables and options ###
#######################################

# GPTL (for benchmarking)
option(GPTL_ENABLED "Enable GPTL" OFF)
set(GPTL_DIR "${HOME}/programs/gptl" CACHE STRING "Location of GPTL installation")

# Parallelization
option(OPENMP_ENABLED "Enable OpenMP" OFF)

# Output names
set(LIB_NUMERIC    "libnumeric"      CACHE STRING "Name of the library to be created")
set(LIB_OMP_SIMD   "libompsimd"      CACHE STRING "Name of the library to be created")

#set(NUMERIC_DRIVER "${PROJECT_NAME}" CACHE STRING "Name of the executable")


#################
### Libraries ###
#################

# TODO: function for finding libraries specified via command line
# - Generalize to finding static and shared versions
# - Use 'target_link_libraries' and 'target_include_directories' instead
#   for finer-grained control?

### GPTL ###

if(GPTL_ENABLED)
	# TODO: make target-wise? put in config file, that also gets included in the plumed patch?
	add_definitions(-DHAVE_GPTL)

	message("-- gptl enabled")
	message("-- gptl location: ${GPTL_DIR}")

	# Look for static library
	set(LIB_GPTL "gptl")
	find_library(LIB_GPTL_STATIC
							 NAMES lib${LIB_GPTL}.a ${LIB_GPTL}.a
							 HINTS "${GPTL_DIR}/lib"
							 )
	if(LIB_GPTL_STATIC)
		message("-- Found static library for ${LIB_GPTL}: ${LIB_GPTL_STATIC}")
	endif()
endif()


### MPI ###

#if(MPI_ENABLED)
#	message("-- MPI enabled")
#	# Add MPI
#	find_package(MPI)
#	if(MPI_FOUND)
#		# TODO Move to config.h?
#		add_definitions(-DMPI_ENABLED)
#		add_definitions(-DHAVE_LIBMPI)
#	else()
#    message(FATAL_ERROR "MPI could not be located")
#	endif()
#endif()


### OpenMP ###

if(OPENMP_ENABLED)
	message(STATUS "OpenMP enabled")
	# Add OpenMP
	find_package(OpenMP)
	if(NOT OpenMP_FOUND)
		if(NOT APPLE)
			message(FATAL_ERROR "OpenMP could not be located")
		else()
			# Mac is notoriously difficult: try to find Homebrew OpenMP semi-manually
			message(STATUS "Will try to find Homebrew OpenMP")
			include(findHomebrewOpenMP)
			findHomebrewOpenMP()
		endif()
	endif()

	# TODO: compiler-dependent option?
	target_compile_options(OpenMP::OpenMP_CXX INTERFACE "-fopenmp-simd")
endif()

# Find a restrict keyword for C++ (compiler-dependent)
include(FindCXXRestrict)
if(HAVE_CXX_RESTRICT)
	string(REPLACE "=" ";" FLAG_LIST ${CXX_RESTRICT_FLAG})
	list(GET FLAG_LIST 1 CXX_RESTRICT_FLAG_NAME)

	message(STATUS "Found restrict keyword for C++ - ${CXX_RESTRICT_FLAG_NAME} (${CXX_RESTRICT_FLAG})")
	message(STATUS "  Will replace the macro - ${CXX_RESTRICT_MACRO}")
else()
	message(STATUS "Could not find a restrict keyword equivalent for C++")
	set(CXX_RESTRICT_FLAG "-D${CXX_RESTRICT_MACRO}=")
endif()


######################
### Compiler Flags ###
######################

# Basic flags
add_compile_options(-g -Wall)
add_definitions(${CXX_RESTRICT_FLAG})
#target_compile_options(${PROJECT_NAME} PUBLIC -g -Wall)

# Optimizations
set(NUMERIC_OPT_DEBUG_FLAGS   -O0)
set(NUMERIC_OPT_RELEASE_FLAGS -O3 -ffast-math -march=native)

# FIXME
#list(APPEND NUMERIC_OPT_RELEASE_FLAGS -fopt-info-vec-missed)
list(APPEND NUMERIC_OPT_RELEASE_FLAGS -fopt-info-vec-optimized)
list(APPEND NUMERIC_OPT_RELEASE_FLAGS -ftree-vectorizer-verbose=1)

#if(GNU)
	# FIXME this does not seem to be working
	#string(APPEND NUMERIC_OPT_RELEASE_FLAGS -ftree-vectorize)
#endif()
add_compile_options("$<$<CONFIG:DEBUG>:${NUMERIC_OPT_DEBUG_FLAGS}>")
add_compile_options("$<$<CONFIG:RELEASE>:${NUMERIC_OPT_RELEASE_FLAGS}>")
#target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:${NUMERIC_OPT_DEBUG_FLAGS}>")
#target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:RELEASE>:${NUMERIC_OPT_RELEASE_FLAGS}>")


#############
### Build ###
#############

#add_definitions(-DNUMERIC_STANDALONE_MODE)
add_subdirectory(src)


############
### Test ###
############

enable_testing()
add_subdirectory(test)
