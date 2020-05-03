
#pragma once

#ifndef C_ARRAY_REF_H
#define C_ARRAY_REF_H

#include <cstddef>
#include <type_traits>

// Wrapper around a C array with dynamic size
// - Its use assumes that you provide it with pointers to:
//    (1) the pointer that *owns* the data
//    (2) the integral value that *owns* its current length
// - Note that the constructor takes *pointers* to both the data pointer and the size
//   - In the event that the data is reallocated (e.g. the length of the C array
//     increases), the location of the data would change. However, the location of the
//     owning pointer would *not* change.
template<typename T, typename I = int>
class CArrayRef
{
  //typename = std::enable_if_t<std::is_integral<I>::value>
	static_assert(std::is_integral<I>::value, "integral type required");
	// TODO: require sizeof(I) <= sizeof(std::size_t)

 public:
	using value_type = T;
	using size_type  = std::size_t;

	using reference       = value_type&;
	using const_reference = const value_type&;
	using pointer         = value_type*;
	using const_pointer   = const value_type*;

	using iterator               = pointer;
	using const_iterator         = const_pointer;
	using reverse_iterator       = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	CArrayRef() {}

	CArrayRef(pointer* data_ptr, const I* size_ptr) {
		set(data_ptr, size_ptr);
	}

	// Number of elements
	size_type size() const {
		return static_cast<size_type>(*size_ptr_);
	}
	size_type length() const {
		return size();
	}

	void set(pointer* data_ptr, const I* size_ptr)
	{
		data_ptr_ = data_ptr;
		size_ptr_ = size_ptr;
#ifndef NDEBUG
		// TODO Check for null, negative sizes, etc.
#endif
	}

	// Element access
	reference       operator[](const size_type i)       { return (*data_ptr_)[i]; }
	const_reference operator[](const size_type i) const { return (*data_ptr_)[i]; }
	reference       operator()(const size_type i)       { return (*data_ptr_)[i]; }
	const_reference operator()(const size_type i) const { return (*data_ptr_)[i]; }

	// Iteration
	iterator       begin()        { return *data_ptr_; }
	const_iterator cbegin() const { return *data_ptr_; }
	iterator       end()          { return *data_ptr_ + *size_ptr_; }
	const_iterator cend()   const { return *data_ptr_ + *size_ptr_; } 

	// Access the underlying data
	pointer       data()       { return *data_ptr_; }
	const_pointer data() const { return *data_ptr_; }

 private:
	pointer* data_ptr_ = nullptr;  // ptr to data
	const I* size_ptr_ = nullptr;  // ptr to size of block (may change!)
};

// For constant data
template<typename T, typename I = int>
using CArrayConstRef = CArrayRef<const T, I>;

#endif // ifndef C_ARRAY_REF_H
