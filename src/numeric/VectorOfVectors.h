
#pragma once
#ifndef VECTOR_OF_VECTORS_H
#define VECTOR_OF_VECTORS_H

#include <cstddef>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

#include "Assert.h"

template<
	typename T,
	typename Vector = std::vector<T>
>
class VectorOfVectors
{
 public:
	static_assert(std::is_same<T, typename Vector::value_type>::value, "type mismatch");

  using value_type     = T;
  using size_type      = typename Vector::size_type;
  using allocator_type = typename Vector::allocator_type;

	using reference       = value_type&;
	using const_reference = const value_type&;
	using pointer         = value_type*;
	using const_pointer   = const value_type*;

	/*
	// TODO: These should be iterators over either the outer vector, or over
	// the contents of *all* the subvectors, in linear fashion
	using iterator               = pointer;
	using const_iterator         = const_pointer;
	using reverse_iterator       = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;
	*/

	using SubvectorIterator             = pointer;
	using SubvectorConstIterator        = const_pointer;
	using SubvectorReverseIterator      = pointer;
	using SubvectorConstReverseIterator = const_pointer;

	VectorOfVectors() {}

	// Initialize with specified subvector capacities
	VectorOfVectors(const size_type size, const std::vector<int>& capacities) {
		assignEmptyWithCapacities(size, capacities);
	}

	// Initialize with uniform subvector capacity
	VectorOfVectors(const size_type size, const size_type subvec_capacity): 
		VectorOfVectors(size, std::vector<int>(size, subvec_capacity))
	{}

	// Current sizes
	// - Outer vector
	size_type size() const {
		return subvec_begins_.size();
	}
	// - Inner vector
	size_type size(const size_type i) const {
		return subvec_ends_[i] - subvec_begins_[i];  //return subvec_sizes_[i];
	}

	// Allocated sizes
	// - Outer vector
	size_type capacity() const {
		return totalCapacity();
	}
	size_type totalCapacity() const {
		return data_.size();
	}
	// - Inner vector
	size_type capacity(const size_type i) const {
		return subvec_alloc_ends_[i];
	}


	// TODO
	void push_back(const size_type i, const T& value) {
	}

	// TODO:
	// - More efficient appending to the end
	void resizeWithCapacities(const size_type new_size, const std::vector<int>& capacities) {
		FANCY_ASSERT(new_size == capacities.size(), "inconsistent input");

		size_type old_size = size();
		if ( new_size == old_size ) {
			return;  // TODO: 
		}

		/*
		size_type total_capacity = 0;
		for ( unsigned i=0; i<new_size; ++i ) {
			total_capacity += capacities[i];
		}
		if ( data_.size() < total_capacity ) {
		}
		*/

		if ( new_size > old_size ) {
			// TODO: growth? should probably handle elsewhere to keep this fxn's scope limited

			// Swap data into a new object, which has the new size and subvec capacities
			VectorOfVectors new_obj(new_size, capacities);
			for ( unsigned i=0; i<size; ++i ) {
				unsigned first = subvec_begins_[i];
				unsigned last  = subvec_ends_[i];

				unsigned k_new = new_obj.subvec_begins_[i];
				for ( unsigned k = first; k != last; ++k, ++k_new ) {
					std::swap( new_obj.data[k_new], this->data_[k] );
				}
			}

			// Steal contents
			std::swap(*this, new_obj);
			return;
		}
		else {  // new_size < old_size
			// TODO: give up space
		}
	}

	// Nuke the contents
	void clear() {
		data_.clear();
		subvec_begins_.clear();
		subvec_ends_.clear();
		subvec_alloc_ends_.clear();
	}

 protected:
	// Assigns the outer vector a new size and clears the inner vectors, ensuring that each
	// subvector has the specified minimum capacity
	void assignEmptyWithCapacities(const size_type size, const std::vector<int>& capacities) {
		FANCY_ASSERT(size == capacities.size(), "inconsistent input");

		if ( size == 0 ) {
			clear();
			return;
		}

		// Set total capacity
		size_type total_capacity = std::accumulate(capacities.begin(), capacities.end(), 0);
		total_capacity *= std::max( data_growth_factor_, 1.0 );  // ensure a nonzero total capacity
		data_.resize(total_capacity);

		// Set up offsets
		// - TODO: grow capacities here?
		subvec_begins_.resize(size);       subvec_begins_[0]     = 0;
		subvec_ends_.resize(size);         subvec_ends_[0]       = 0;
		subvec_alloc_ends_.resize(size);   subvec_alloc_ends_[0] = capacities[0];
		for ( unsigned i=1; i<size; ++i ) {
			subvec_begins_[i]     = subvec_alloc_ends_[i-1];
			subvec_ends_[i]       = subvec_begins_[i];
			subvec_alloc_ends_[i] = subvec_begins_[i] + capacities[i];
		}
	}

 private:
	Vector data_ = {};

	// Subvector arrangement
	std::vector<int> subvec_begins_     = {};  // offsets in 'data_' where each subvector begins
	std::vector<int> subvec_ends_       = {};  // one past the last element of each subvector
	std::vector<int> subvec_alloc_ends_ = {};  // one past where the subvector *capacities* end

	// TODO: Worth storing sizes explicitly for fast retrieval?
	//   subvec_sizes_[i] = subvec_ends_[i] = subvec_begins[i]
	//std::vector<int> subvec_sizes_  = {};

	int default_subvec_capacity_ = 10;

	// These determine how quickly the outer and inner vectors grow upon reallocation
	double data_growth_factor_   = 1.1;
	double subvec_growth_factor_ = 1.1;
};




#endif //ifndef VECTOR_OF_VECTORS_H
