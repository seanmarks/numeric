### LIB_OMP_SIMD ###

# AUTHOR: Sean M. Marks (https://github.com/seanmarks)

# Builds main library: ${LIB_OMP_SIMD}

# Locate source files
file(GLOB_RECURSE LIB_OMP_SIMD_SRC "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

# Declare
add_library(${LIB_OMP_SIMD} OBJECT ${LIB_OMP_SIMD_SRC})

# OpenMP
if(OPENMP_ENABLED)
	target_link_libraries(${LIB_OMP_SIMD} PUBLIC OpenMP::OpenMP_CXX)
endif()

# GPTL
if(GPTL_ENABLED)
	target_include_directories(${LIB_OMP_SIMD} PUBLIC "${GPTL_DIR}/include")

	# Prefer static linkage
	if(LIB_GPTL_STATIC)
		target_link_libraries(${LIB_OMP_SIMD} PUBLIC ${LIB_GPTL_STATIC})
	else()
		# Fall back on dynamic linkage of GPTL
		target_link_libraries(${LIB_OMP_SIMD} PUBLIC ${LIB_GPTL})
		target_link_directories(${LIB_OMP_SIMD} PUBLIC "${GPTL_DIR}/lib")
	endif()
endif()

# Where the libraries should be installed
install(TARGETS ${LIB_OMP_SIMD} DESTINATION lib)



#function(targetLinkLibraryPreferStatic TARGET SCOPE LIBNAME)
#	if(LIB_${LIBNAME}_STATIC)
#		target_link_libraries(${TARGET} ${SCOPE} ${LIB_${LIBNAME}_STATIC})
#	else()
#		target_link_libraries(${TARGET} ${SCOPE} ${LIB_${LIBNAME}})
#		target_link_directories(${TARGET} ${SCOPE} "${${LIBNAME}_DIR}/lib")
#	endif()
#endfunction()
