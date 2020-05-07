
#pragma once
#ifndef VECTOR_OF_ARRAYS_H
#define VECTOR_OF_ARRAYS_H

#include <cassert>
#include <cstddef>
#include <exception>
#include <stdexcept>
#include <vector>

template<
	typename T,
	template<typename, typename> class VectorType = std::vector,
	template<typename> class AllocType = std::allocator
>
class VectorOfArrays
{
 public:
	// TODO change 'iterator' type to VectorType::iterator, etc.?
	using iterator        = T*;
	using const_iterator  = const T*;
	using size_type       = std::size_t;

	VectorOfArrays() {}

	VectorOfArrays(const size_type size) {
		resize(size);
	}

	// Get/set the length of the vector
	void resize(const size_type new_size);
	size_type size()   const { return size_; }
	size_type length() const { return size_; }

	// Clears the sub-arrays (resets their sizes to zero) without
	// affecting the amount or layout of internal data
	void clearArrays() {
		array_sizes_.assign( array_sizes_.size(), 0 );
	}

	void clear() {
		clearArrays();
		data_.clear();
	}


	// Query information about the 'i'th array
	size_type getArraySize(const size_type i) const {
		return array_sizes_[i];
	}
	size_type getArrayCapacity(const size_type i) const {
		return array_offsets_[i+1] - array_offsets_[i];
	}

	// Access/modify 'j'th entry of 'i'th array
	T& operator()(const size_type i, const size_type j) {
		return data_[array_offsets_[i] + j];
	}
	const T& operator()(const size_type i, const size_type j) const {
		return data_[array_offsets_[i] + j];
	}

	// Access underlying data
	VectorType<T,AllocType<T>>& data() {
		return data_;
	}
	const VectorType<T,AllocType<T>>& data() const {
		return data_;
	}

	// Get/set the default capacity allocated for new sub-arrays
	void set_default_array_capacity(const size_type capacity) {
		default_array_capacity_ = capacity;
	}
	size_type get_default_array_capacity() {
		return default_array_capacity_;
	}

	// TODO generalize to any stream
	void printLayout() const {
		// DEBUG
		std::cout << "size = " << size_ << "\n";
		for ( unsigned i=0; i<size_; ++i ) {
			std::cout << "  i=" << i << ":  "
			          << "size=" << array_sizes_[i] << ", "
			          << "offset=" << array_offsets_[i] << ", "
			          << "capacity=" << getArrayCapacity(i) << "\n";
		}
		std::cout << "array_offsets_[size] = " << array_offsets_[size_] << "\n"
		          << "data().size() = " << data().size() << "\n";
	}

 private:
	VectorType<T,AllocType<T>> data_;       // underlying storage
	size_type  size_ = 0;   // number of arrays

	size_type default_array_capacity_ = 10;
	//double growth_factor_ = 1.2;

	// Start of each array; size_+1 element is data_.size()
	VectorType<size_type,AllocType<size_type>> array_offsets_ = {0};

	// Size of each array
	VectorType<size_type,AllocType<size_type>> array_sizes_ = {};
};


template<
	typename T, template<typename, typename> class VectorType, template<typename> class AllocType
>
void VectorOfArrays<T,VectorType,AllocType>::resize(const size_type new_size) {
	if ( new_size == size_ ) {
		return;  // nothing to do
	}
	// Resize helper arrays
	array_offsets_.resize(new_size+1);
	array_sizes_.resize(new_size);

	size_type default_data_size = default_array_capacity_*new_size;

	// TODO: When reducing vector size, this can lead to the last array having a very large capacity
	// - Would it be better to evenly distribute this "excess" space evenly?
	if ( new_size > size_ ) {
		// Update underlying storage (but don't allow it to get out of control)
		size_type new_data_size = data_.size() + default_array_capacity_*(new_size - size_);
		new_data_size = std::min(new_data_size, 2*default_data_size);
		data_.resize( new_data_size );
		//new_data_size = static_cast<size_type>( static_cast<double>(new_data_size)*growth_factor_ );

		// Initialize new arrays
		for ( size_type i = size_; i < new_size; ++i ) {
			array_sizes_[i] = 0;
			if ( i >= 1 ) {
				array_offsets_[i] = array_offsets_[i-1] + default_array_capacity_;
			}
		}
	}
	else { // new_size < size_
		// For more large/complicated data types, call std::swap() on a default-constructed 
		// object to "free" some of their storage (TODO?)
		//for ( size_type i=; i<
	}
	array_offsets_.back() = data_.size();
	size_ = new_size;
}

#endif // ifndef VECTOR_OF_ARRAYS_H
