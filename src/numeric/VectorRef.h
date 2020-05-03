
#pragma once

#ifndef VECTOR_REF_H
#define VECTOR_REF_H

#include <cstddef>
#include <type_traits>

// Wrapper around a vector-like block of memory
// - Owns neither the data nor the size
// TODO: Dealing with const/non-const T
template<
	typename T, typename I = int
>
class VectorRef
{
  //typename = std::enable_if_t<std::is_integral<I>::value>
	static_assert(std::is_integral<I>::value, "integral type required");

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

	VectorRef() {}

	//
	VectorRef(pointer* data_ptr, const I* size_ptr):
		data_ptr_(data_ptr), size_ptr_(size_ptr)
	{}

	bool null() const { return (data_ptr_ == nullptr); }

	size_type size() const {
		return static_cast<size_type>(*size_ptr_);
	}

	reference       operator[](const size_type i)       { return (*data_ptr_)[i]; }
	const_reference operator[](const size_type i) const { return (*data_ptr_)[i]; }
	reference       operator()(const size_type i)       { return (*data_ptr_)[i]; }
	const_reference operator()(const size_type i) const { return (*data_ptr_)[i]; }

	iterator       begin()        { return *data_ptr_; }
	const_iterator cbegin() const { return *data_ptr_; }
	iterator       end()          { return *data_ptr_ + *size_ptr_; }
	const_iterator cend()   const { return *data_ptr_ + *size_ptr_; } 

 private:
	pointer* data_ptr_ = nullptr;  // ptr to ptr to block
	const I* size_ptr_ = nullptr;  // ptr to size of block (may change!)
};

#endif // ifndef VECTOR_REF_H
