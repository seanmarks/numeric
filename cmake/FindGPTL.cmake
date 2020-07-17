# FindGPTL - locates the GPTL library (https://jmrosinski.github.io/GPTL/)
#
# AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#
# INPUT: optional variables
#   GPTL_DIR                - hint for the install location
#   PREFER_STATIC_LIBS      - if set, also sets GPTL_PREFER_STATIC_LIBS
#   GPTL_PREFER_STATIC_LIBS - prefer to link to a static library
#
# OUTPUT: defines the following variables
#   GPTL_FOUND        - ON/OFF for success/failure
#   GPTL_LIBRARIES    - list of libraries required to link to GPTL
#   GPTL_INCLUDE_DIRS - list of directories that must be included
#
# SOURCES
# - Talk by Daniel Pfeifer at C++Now 2017
# - CMake tutorial (https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Find-Libraries)
# - FindBoost.cmake (cmake 3.15)

# TODO: depending on compilation settings, check for -fPIC

# Optionally, prefer static library
if(PREFER_STATIC_LIBS)
	set(GPTL_PREFER_STATIC_LIBS ON)
endif()
if(GPTL_PREFER_STATIC_LIBS)
	# Save CMake variable
	set(_GPTL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

	# Modify suffix preferences
	if(WIN32)
		list(INSERT CMAKE_FIND_LIBRARY_SUFFIXES 0 .lib .a)
	else()
		set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
	endif()
endif()

if(GPTL_DIR)
	# Use the hint provided
	list(APPEND CMAKE_PREFIX_PATH ${GPTL_DIR})
endif()

# Look for GPTL
find_path(GPTL_INCLUDE_DIR gptl.h
	HINTS ${GPTL_DIR}
)
find_library(GPTL_LIBRARY GPTL
	HINTS ${GPTL_DIR}
	NAMES libgptl gptl
)
mark_as_advanced(GPTL_INCLUDE_DIR GPTL_LIBRARY)

# Make sure everything works
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( GPTL
	REQUIRED_VARS GPTL_INCLUDE_DIR GPTL_LIBRARY
	)

if(GPTL_FOUND)
	# Important variables
	set(GPTL_LIBRARIES    ${GPTL_LIBRARY})
	set(GPTL_INCLUDE_DIRS ${GPTL_INCLUDE_DIR})
endif()

if(GPTL_FOUND AND NOT TARGET GPTL::GPTL)
	# Create a target under the appropriate "namespace"
	add_library(GPTL::GPTL UNKNOWN IMPORTED)
	set_target_properties(GPTL::GPTL PROPERTIES
		IMPORTED_LOCATION "${GPTL_LIBRARY}"
		IMPORTED_LINK_INTERFACE_LANGUAGES "C"
		INTERFACE_INCLUDE_DIRECTORIES "${GPTL_INCLUDE_DIR}"
	)
endif()

# Reset CMake variable
if(GPTL_PREFER_STATIC_LIBS)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${_GPTL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()
if(GPTL_DIR)
	list(REMOVE_AT CMAKE_PREFIX_PATH -1)
endif()
