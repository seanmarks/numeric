# Testing

#############
### Setup ###
#############

# Find bash
find_program (BASH_PROGRAM bash)
if(NOT BASH_PROGRAM)
	message(FATAL_ERROR "Unable to find bash for testing")
endif()

# Testing driver script (absolute path)
set(TEST_DRIVER_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/run_test.sh)

set(ECHO_FAILED_DIFFS "-diff" CACHE STRING "tests: whether to print diffs for tests that fail")

# Directory with testing files
# - This is inside <repo>/test, *not* the build directory
set(TEST_DIR_ROOT ${CMAKE_CURRENT_LIST_DIR})


#############
### Tests ###
#############

### Builtin Tests ###

# Root directory where testing binaries will be built
set(TEST_BIN_ROOT ${CMAKE_BINARY_DIR}/test/bin)

# Tests are differentiated by their directories

set(NUMERIC_TEST_NAMES "")
list(APPEND NUMERIC_TEST_NAMES "basic")
list(APPEND NUMERIC_TEST_NAMES "omp-simd/real")
list(APPEND NUMERIC_TEST_NAMES "omp-simd/complex")

foreach(TEST_NAME IN LISTS NUMERIC_TEST_NAMES)
	# TODO: move to function: build test driver with custom 'main' (link to library)

	string(REPLACE "/" "-" TEST_PROGRAM ${TEST_NAME})
	set(TEST_PROGRAM "test_${TEST_PROGRAM}")
	message("TEST_PROGRAM = ${TEST_PROGRAM}")

	set(TEST_DIR ${TEST_DIR_ROOT}/${TEST_NAME})   # directory with source files for this test
	set(TEST_BIN ${TEST_BIN_ROOT}/${TEST_NAME})   # where to build this testing binary
	add_executable(${TEST_PROGRAM} ${TEST_DIR}/main.cpp)
	set_target_properties(${TEST_PROGRAM}
												PROPERTIES
												RUNTIME_OUTPUT_DIRECTORY ${TEST_BIN}
												)
	target_link_libraries(${TEST_PROGRAM} PUBLIC ${LIB_NUMERIC})
	target_include_directories(${TEST_PROGRAM} PUBLIC ${CMAKE_SOURCE_DIR}/src)

	# - Force target 'check' (i.e. 'make check') to depend on building the test executable;
	#   otherwise, it won't get built
	#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TEST_PROGRAM})
	add_test(NAME ${TEST_NAME}
					 WORKING_DIRECTORY ${TEST_BIN}
					 COMMAND ${TEST_PROGRAM}
					 )
endforeach()


set(OMP_SIMD_TEST_NAMES "")
list(APPEND OMP_SIMD_TEST_NAMES "old-omp-simd")

foreach(TEST_NAME IN LISTS OMP_SIMD_TEST_NAMES)
	# TODO: move to function: build test driver with custom 'main' (link to library)

	string(REPLACE "/" "-" TEST_PROGRAM ${TEST_NAME})
	set(TEST_PROGRAM "test_${TEST_PROGRAM}")
	message("TEST_PROGRAM = ${TEST_PROGRAM}")

	set(TEST_DIR ${TEST_DIR_ROOT}/${TEST_NAME})   # directory with source files for this test
	set(TEST_BIN ${TEST_BIN_ROOT}/${TEST_NAME})   # where to build this testing binary
	add_executable(${TEST_PROGRAM} ${TEST_DIR}/main.cpp)
	set_target_properties(${TEST_PROGRAM}
												PROPERTIES
												RUNTIME_OUTPUT_DIRECTORY ${TEST_BIN}
												)
	target_link_libraries(${TEST_PROGRAM} PUBLIC ${LIB_OMP_SIMD})
	target_include_directories(${TEST_PROGRAM} PUBLIC ${CMAKE_SOURCE_DIR}/src)

	# - Force target 'check' (i.e. 'make check') to depend on building the test executable;
	#   otherwise, it won't get built
	#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TEST_PROGRAM})
	add_test(NAME ${TEST_NAME}
					 WORKING_DIRECTORY ${TEST_BIN}
					 COMMAND ${TEST_PROGRAM}
					 )
endforeach()
