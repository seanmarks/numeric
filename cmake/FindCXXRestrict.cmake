# Check for a compiler-dependent equivalent of the C keyword 'restrict'
# that will work with C++ code
#
# AUTHOR: Sean M. Marks (https://github.com/seanmarks)

include(WriteCXXSource)

set(CXX_RESTRICT_TEST_SOURCE
"
double test_add(double* CXX_RESTRICT a, double* CXX_RESTRICT b) {
	return *a + *b;
}
int main(void) {
	double a = 1, b = 2;
	double c = test_add(&a, &b);
}
")

_write_cxx_source_file("${CXX_RESTRICT_TEST_SOURCE}" test_cxx_restrict test_flags CXX_RESTRICT_TEST_SRC)

# Flags to test
set(CXX_RESTRICT_MACRO "CXX_RESTRICT")
set(CXX_RESTRICT_TEST_FLAGS "")
list(APPEND CXX_RESTRICT_TEST_FLAGS "-D${CXX_RESTRICT_MACRO}=__restrict")
list(APPEND CXX_RESTRICT_TEST_FLAGS "-D${CXX_RESTRICT_MACRO}=__restrict__")

foreach(CXX_RESTRICT_TEST_FLAG IN LISTS CXX_RESTRICT_TEST_FLAGS)
	try_compile(HAVE_CXX_RESTRICT ${CMAKE_BINARY_DIR} ${CXX_RESTRICT_TEST_SRC}
							COMPILE_DEFINITIONS "${CXX_RESTRICT_TEST_FLAG}"
							OUTPUT_VARIABLE     CXX_RESTRICT_TEST_OUTPUT)

	if(HAVE_CXX_RESTRICT)
		set(CXX_RESTRICT_FLAG "${CXX_RESTRICT_TEST_FLAG}")
		break()
	endif()
endforeach()

unset(CXX_RESTRICT_TEST_SOURCE)
unset(CXX_RESTRICT_TEST_FLAGS)
