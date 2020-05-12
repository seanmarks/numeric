### LIB_ALT_SIMD ###

# Builds main library: ${LIB_ALT_SIMD}

# Locate source files
file(GLOB_RECURSE LIB_ALT_SIMD_SRC "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

# Declare
add_library(${LIB_ALT_SIMD} STATIC ${LIB_ALT_SIMD_SRC})

# OpenMP
if(OPENMP_ENABLED)
	target_link_libraries(${LIB_ALT_SIMD} PUBLIC OpenMP::OpenMP_CXX)
endif()

# GPTL
if(GPTL_ENABLED)
	target_include_directories(${LIB_ALT_SIMD} PUBLIC "${GPTL_DIR}/include")

	# Prefer static linkage
	if(LIB_GPTL_STATIC)
		target_link_libraries(${LIB_ALT_SIMD} PUBLIC ${LIB_GPTL_STATIC})
	else()
		# Fall back on dynamic linkage of GPTL
		target_link_libraries(${LIB_ALT_SIMD} PUBLIC ${LIB_GPTL})
		target_link_directories(${LIB_ALT_SIMD} PUBLIC "${GPTL_DIR}/lib")
	endif()
endif()

# Where the libraries should be installed
install(TARGETS ${LIB_ALT_SIMD} DESTINATION lib)
