### LIB_NUMERIC ###
#
# AUTHOR: Sean M. Marks (https://github.com/seanmarks)

# Builds main library: numeric

# Standard install locations
include(GNUInstallDirs)

# Declare library
if(BUILD_SHARED_LIBS)
	add_library(numeric SHARED "")
else()
	# Need position-independent code to be able to link to a shared library
	# - TODO: flag to override this depending on intended purpose (e.g. "LINK_ALL_STATIC")
	add_library(numeric STATIC "")
	set_target_properties(numeric PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# Add source files
file(GLOB_RECURSE LIB_NUMERIC_SRC "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")
target_sources(numeric PRIVATE ${LIB_NUMERIC_SRC})

# FIXME DEBUG
if(VERBOSE_VECTORIZE_FLAGS)
	# Verbose vectorizer output only for specific files
	set_source_files_properties(
			SphericalHarmonics.cpp
		PROPERTIES
			COMPILE_OPTIONS "${VERBOSE_VECTORIZE_FLAGS}"
	)
endif()

# GPTL
if(GPTL_ENABLED)
	target_compile_definitions(numeric PUBLIC HAVE_GPTL=1)
	target_link_libraries(numeric PUBLIC GPTL::GPTL)
	target_include_directories(numeric PUBLIC ${GPTL_INCLUDE_DIRS})

	if(GPTL_AUTOPROFILE)
		# TODO: flags on different systems
		# TODO: public or private?
		target_compile_options(numeric PRIVATE -finstrument-functions)
		#target_compile_options(numeric PUBLIC -finstrument-functions)
		target_link_options(numeric PUBLIC -finstrument-functions)
		target_link_options(numeric PUBLIC -rdynamic)
	endif()
endif()

# OpenMP
if(OPENMP_ENABLED)
	target_compile_definitions(numeric PUBLIC HAVE_OPENMP=1)
	target_link_libraries(numeric PUBLIC OpenMP::OpenMP_CXX)
endif()

# MPI
if(MPI_ENABLED)
	target_compile_definitions(numeric PUBLIC HAVE_MPI=1)
	target_compile_definitions(numeric PUBLIC MPI_ENABLED=1)
endif()

# Where the libraries should be installed
install(TARGETS numeric DESTINATION lib)
