// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#include "AlignedAllocator.h"

namespace numeric {
namespace aligned {

// Aligned malloc
void* malloc(const size_type size, const size_type alignment)
{
	// Require alignment/sizeof(void*) = 2^n for some positive integer 'n'
	size_type n = alignment/(sizeof(void*));
	bool is_power_of_2 = ( (n > 0) and ((n & (n-1)) == 0) );
	if ( not is_power_of_2 ) {
		throw std::runtime_error("Error in aligned::malloc - alignment/sizeof(void*) must be a power of 2");
	}

	if ( n == 0 ) {
		return nullptr;
	}

	// Allocate memory
	int ret;
	void* mem_ptr = nullptr;
#if defined(__GNUC__) || defined(__APPLE__)
	ret = posix_memalign(&mem_ptr, alignment, size);
	if ( ret != 0 ) {
		if ( mem_ptr != nullptr ) { free(mem_ptr); }
		throw std::runtime_error("posix_memalign: ret is nonzero");
	}
#elif defined(__intel__)
	mem_ptr = _mm_alloc(size, alignment);
#elif defined(_WIN32)
	mem_ptr = _aligned_malloc(size, alignment);
#else
	// TODO As a fallback: allocate memory_needed+alignment bytes,
	//  and find an address that matches the criteria
	static_assert(false, "Unsupported aligned malloc");
#endif

	return mem_ptr;
}


// aligned free
void free(void* mem_ptr) noexcept
{
#if defined(__GNUC__) || defined(__APPLE__)
	// Memory allocated by 'posix_memalign()' is deallocated by the usual 'free()'
	::free(mem_ptr);
#elif defined(__intel__)
	_mm_free(mem_ptr);
#elif defined(_WIN32)
	_aligned_free(mem_ptr);
#else
	::free(mem_ptr);
#endif
}


} // end namespace aligned {
} // end namespace numeric {
