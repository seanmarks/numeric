// Aligned dynamic memory allocation
//
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// - Provides general, C-style wrappers around system-specified functions that 
//   allocate/deallocate aligned memory
// - Also provides a std::allocator-style object that can be used with other
//   C++ objects in the STL, such as std::vector
//
// - The code here is based on the folloing prior works
//   - 'Mallocator', by Stephen T. Lavavej
//      - Original: TODO
//      - Updated: slides presented at CppCon 2014 
//        (under https://github.com/CppCon/CppCon2014/blob/master/Presentations)
// 
//   - 'Aligned', by Nick Strupat
//     - GitHub: https://github.com/NickStrupat/Aligned
//
// - Misc. notes
//   - 'construct()' and 'destruct()' are not implemented, as they are redundant in C++11
//

#pragma once
#ifndef ALIGNED_ALLOCATOR_H
#define ALIGNED_ALLOCATOR_H


// TODO C++17 has hardware_destructive_interference_size/hardware_destructive_interference_size
// for avoiding false sharing
// - OpenMP SIMD requires cache-lined
// TODO set from Makefile/CMake to ensure system-dependent
// - For now: a 64-byte cache alignment is a reasonable assumption
#ifndef CACHE_LINE_SIZE
#  define CACHE_LINE_SIZE 64
#endif // ifndef CACHE_LINE_SIZE


// The following headers are required for all allocators.
#include <cstddef>    // size_t and ptrdiff_t
#include <exception>
#include <memory>     // allocator_traits
#include <new>        // placement new, bad_alloc, bad_array_new_length
#include <stdexcept>  // std::length_error

#include <cstdlib>    // malloc() and free()
#include <iostream>   // std::cout
#include <limits>
#include <ostream>    // std::endl

#include "System.h"

namespace numeric {
namespace aligned {

using size_type = std::size_t;

// Aligned malloc
void* malloc(const size_type size, const size_type alignment);

// Aligned free
void free(void* mem_ptr) noexcept;

// STL-compatible allocator for aligned memory
// - If template parameter 'Alignment' is zero or unspecified,
//   data is aligned according to the size of the cache line
template <typename T, size_type Alignment = CACHE_LINE_SIZE>
class Allocator {
 public:
	// The following will be the same for virtually all allocators.
	using value_type      = T;
	using pointer         = T*;
	using const_pointer   = const T*;
	using reference       = T&;
	using const_reference = const T&;
	using difference_type = std::ptrdiff_t;

	// When rebinding, use the same alignment
	template <typename U> 
	struct rebind {
		using other = Allocator<U>;
	};

	// Returns true if and only if storage allocated from *this
	// can be deallocated from other, and vice versa.
	bool operator==(const Allocator& other) const noexcept {
		return true;
	}

	bool operator!=(const Allocator& other) const noexcept {
		// Given the implementation of operator== above, this always returns false
		return (not ( *this == other ));  
	}

	// Default constructor: use template argument as alignment
	Allocator() {
		setAlignment(Alignment);
	}

	// Copy constructor: copy alignment from other allocator
	Allocator(const Allocator& other) {
		setAlignment(other.alignment_);
	}

	// Rebinding constructor, for making allocators of different types
	template <typename U> 
	Allocator(const Allocator<U>& other) {
		setAlignment(other.alignment_);
	}

	~Allocator() { }

	pointer allocate(const size_type n) const {
		void* const pv = aligned::malloc(n * sizeof(value_type), alignment_);

		if ( pv == nullptr ) {
			// Allocation failed
			throw std::bad_alloc();
		}

		return static_cast<pointer>(pv);
	}

	void deallocate(pointer const p, const size_type n) const noexcept {
		aligned::free( static_cast<void*>(p) );
	}

	size_type getAlignment() const noexcept {
		return alignment_;
	}

 private:
	size_type alignment_ = 0;

	void setAlignment(const size_type alignment) {
		if ( alignment != 0 ) {
			alignment_ = alignment;
		}
		else {
			alignment_ = System::get_cache_line_size();
		}
	}

	// Allocators are not required to be assignable, so all allocators should have a private 
	// unimplemented assignment operator
	Allocator& operator=(const Allocator&);
};

// For readability
template<typename T>
using CacheAlignedAllocator = Allocator<T>;

} // end namespace aligned
} // end namespace numeric

#endif // ifndef ALIGNED_ALLOCATOR_H
