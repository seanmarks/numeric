# AUTHOR: Sean M. Marks (https://github.com/seanmarks)

# Each subdirectory builds a library

add_subdirectory(numeric)
add_subdirectory(omp-simd)
#add_subdirectory(alt-simd)


### Driver ###

# Add the executable
#add_executable(${NUMERIC_DRIVER} driver/main.cpp driver/main.h)
#target_link_libraries(${NUMERIC_DRIVER} PUBLIC ${LIB_NUMERIC})
#target_include_directories(${NUMERIC_DRIVER} PUBLIC ${CMAKE_SOURCE_DIR}/src)

#if(APPLE)
#	# Run 'codesign' on executables built on Mac
#	# - If MPI is enabled, this will prevent the endless firewall prompts about allowing
#	#   incoming connections (it will only appear when the program is first run)
#	add_custom_command(TARGET ${NUMERIC_DRIVER}
#	                   POST_BUILD
#	                   WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
#	                   COMMAND codesign --force --deep --sign - ${NUMERIC_DRIVER}
#	                   )
#endif()

# Install rules
#install(TARGETS ${NUMERIC_DRIVER} DESTINATION bin)
