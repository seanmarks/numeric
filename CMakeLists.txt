# AUTHOR: Sean M. Marks (https://github.com/seanmarks)

cmake_minimum_required(VERSION 3.15)

# Include directory with custom CMake modules
# - These contain useful helper functions
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(setIfUnset)
include(WriteCXXSource)

project(Numeric)

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
set(LIB_OMP_SIMD  "libompsimd"  CACHE STRING "Name of the library to be created")
#set(LIB_ALT_SIMD  "libaltsimd"  CACHE STRING "Name of the library to be created")

#set(LIB_NUMERIC   "libnumeric"  CACHE STRING "Name of the library to be created")
#set(NUMERIC_DRIVER "${PROJECT_NAME}" CACHE STRING "Name of the executable")


#################
### Libraries ###
#################

# Link options
option(PREFER_STATIC_LIBS "Prefer static libraries when linking" ON)
if(BUILD_SHARED_LIBS)
	# FIXME: reasonable?
	set(PREFER_STATIC_LIBS OFF)
endif()
if(PREFER_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
	message(STATUS "Will prefer static libraries")
	set(xdrfile_PREFER_STATIC_LIBS ON)
	set(GPTL_PREFER_STATIC_LIBS    ON)
endif()

# GPTL (for benchmarking)
option(GPTL_ENABLED     "Enable GPTL" OFF)
option(GPTL_AUTOPROFILE "Enable GPTL autoprofiling" OFF)
if(GPTL_ENABLED)
	message(STATUS "GPTL enabled")
	include(FindGPTL)
	find_package(GPTL)
	if (GPTL_FOUND)
		message(STATUS "GPTL library: ${GPTL_LIBRARIES}")
		if(GPTL_AUTOPROFILE)
			message(STATUS "GPTL will autoprofile code")
		endif()
	endif()
endif()

# MPI
option(MPI_ENABLED "Enable MPI" OFF)
if(MPI_ENABLED)
	message(STATUS "MPI enabled")
	find_package(MPI REQUIRED)
endif()

# OpenMP
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
			include(FindHomebrewOpenMP)
			FindHomebrewOpenMP()
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
#list(APPEND NUMERIC_OPT_RELEASE_FLAGS -ftree-vectorize)

# FIXME DEBUG
# Control verbosity of feedback from auto-vectorizer
option(VERBOSE_VECTORIZE "" OFF)
set(VERBOSE_VECTORIZE_FLAGS "")
# FIXME: if(GNU) not being recognized
if(TRUE)
	list(APPEND VERBOSE_VECTORIZE_FLAGS -fopt-info-vec-optimized)
	list(APPEND VERBOSE_VECTORIZE_FLAGS -ftree-vectorizer-verbose=1)

	if(VERBOSE_VECTORIZE)
		list(APPEND NUMERIC_OPT_RELEASE_FLAGS ${VERBOSE_VECTORIZE_FLAGS})
	endif()
endif()


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
