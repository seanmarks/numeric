# AUTHOR: Sean M. Marks (https://github.com/seanmarks)

# Based on helper function '_OPENMP_WRITE_SOURCE_FILE' in FindOpenMP.cmake (from CMake 3.15)
# - Input:
#     SRC_FILE_CONTENT_VAR - contains a string with the C++ source to write
#     SRC_FILE_NAME        - name of source file to be generated --> ${SRC_FILE_NAME}.cpp
#     SRC_FILE_DIR         - where the source file should be created under ${CMAKE_FILES_DIRECTORY}
# - Output:
#     SRC_FILE_FULL_PATH   - full path to created source file
function(_write_cxx_source_file  SRC_FILE_CONTENT_VAR  SRC_FILE_NAME  SRC_FILE_DIR  SRC_FILE_FULL_PATH)
	set(WORK_DIR ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${SRC_FILE_DIR})
	set(SRC_FILE "${WORK_DIR}/${SRC_FILE_NAME}.cpp")
	file(WRITE "${SRC_FILE}" "${SRC_FILE_CONTENT_VAR}")

	set(${SRC_FILE_FULL_PATH} "${SRC_FILE}" PARENT_SCOPE)
endfunction()
