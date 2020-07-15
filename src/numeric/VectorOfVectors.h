// VectorOfVectors
// - Vector of subvectors with varying lengths
//   - Implemented in linear storage
//
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// - TODO:
//   - Move implementations out of class body for readability

#pragma once
#ifndef VECTOR_OF_VECTORS_H
#define VECTOR_OF_VECTORS_H

#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

#include "Assert.h"


template<typename T, typename Vector = std::vector<T>>
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

	// For iterating over a particular subvector
	using SubvectorIterator             = pointer;
	using SubvectorConstIterator        = const_pointer;
	using SubvectorReverseIterator      = pointer;
	using SubvectorConstReverseIterator = const_pointer;

	// Initialize empty
	VectorOfVectors() {}

	// Initialize with specified subvector capacities
	// - Optionally, also specify the minimum total capacity the object should have
	// - Will allocate more memory than the minimum
	VectorOfVectors(const std::vector<size_type>& capacities, const size_type total_capacity = 0) {
		assignEmptyWithCapacities(capacities, total_capacity);
	}

	// Initialize as 'size' empty subvectors (with some default uniform capacity)
	VectorOfVectors(const size_type size):
		VectorOfVectors(size, initial_default_subvec_capacity_)
	{}

	// Initialize as 'size' empty subvectors with uniform capacity
	VectorOfVectors(const size_type size, const size_type subvec_capacity): 
		VectorOfVectors(std::vector<size_type>(size, subvec_capacity))
	{}

	// Current sizes
	// - Outer vector
	size_type size() const {
		return subvec_begins_.size();
	}
	// - Inner vector
	size_type size(const size_type i) const {
		FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
		return subvec_ends_[i] - subvec_begins_[i];  //return subvec_sizes_[i];
	}

	// Access jth element of ith subvector
	// - TODO: noexcept iff NDEBUG
	reference       operator()(const size_type i, const size_type j) noexcept;
	const_reference operator()(const size_type i, const size_type j) const noexcept;

	// Iterate over subvector 'i'
	// - TODO: noexcept iff NDEBUG
	SubvectorIterator      begin (const size_type i)       noexcept;
	SubvectorIterator      end   (const size_type i)       noexcept;
	SubvectorConstIterator begin (const size_type i) const noexcept;
	SubvectorConstIterator end   (const size_type i) const noexcept;
	SubvectorConstIterator cbegin(const size_type i) const noexcept;
	SubvectorConstIterator cend  (const size_type i) const noexcept;

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
		FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
		return subvec_alloc_ends_[i] - subvec_begins_[i];
	}

	// Total number of elements in use
	size_type count() const;

	/*
	// TODO: Add a new subvector to the end
	void push_back() {
		FANCY_ASSERT(false, "not yet implemented");
	}
	*/

	// Append an element to the 'j'th subvector
	// - May cause a reallocation
	void push_back(const size_type i, const T& value);

	// Resize the 'i'th subvector to 'new_size'
	// - May cause a reallocation
	void resize(const size_type i, const size_type new_size);

	// Extend the capacity of the 'i'th subvector to be at least 'new_capacity'
	void extendCapacity(const size_type i, const size_type new_capacity);

	/*
	void push_back_at_end(const T& value) {
		size_type excess_capacity = getExcessCapacity();
	}
	*/

	// Resizes the outer vector
	// - If 'new_size' is larger than the current size, new empty subvectors (with
	//   default capacity) will be appended
	// - May cause a reallocation
	void resize(const size_type new_size);

	// Resize the vector and enforce that the subvectors satisfy the specified minimum capacities.
	// - Upon completion:  size() == capacities.size()
	// - May cause a reallocation
	void resizeWithMinimumCapacities(const std::vector<size_type>& capacities);

	// Reset the contents of subvector 'i'
	void reset(const size_type i) {
		subvec_ends_[i] = subvec_begins_[i];
	}
	void clear(const size_type i) {
		reset(i);
	}

	// Reset all subvectors
	void reset() {
		const size_type len = size();
		for ( size_type i=0; i<len; ++i ) {
			reset(i);
		}
	}

	// Nuke the contents
	void clear() {
		data_.clear();
		subvec_begins_.clear();
		subvec_ends_.clear();
		subvec_alloc_ends_.clear();
	}


	//----- Tuning parameters -----//

	void setDefaultSubvectorCapacity(const size_type capacity) {
		default_subvec_capacity_ = capacity;
	}

	void setOverallGrowthFactor(const double fac) {
		FANCY_ASSERT(data_growth_factor_ >= 1.0, "invalid value");
		data_growth_factor_ = fac;
	}

	void setSubvectorGrowthFactor(const double fac) {
		FANCY_ASSERT(subvec_growth_factor_ >= 1.0, "invalid value");
		subvec_growth_factor_ = fac;
	}


	//----- Debugging -----//

	// Dump lots of info to the indicated stream
	std::ostream& dump(std::ostream& os) const;

	// Performs lots of internal santity checks. Throws an exception if any of them fails.
	void checkInternalConsistency() const;


 protected:
	// Assigns the outer vector a new size and clears the inner vectors, ensuring that each
	// subvector has the specified minimum capacity
	// - Outer size is inferred from 'capacities.size()'
	// - Optionally, also enforces a minimum total capacity
	// - All stored data is lost!
	void assignEmptyWithCapacities(
		const std::vector<size_type>& capacities,
		const size_type total_capacity = 0
	);	

	// Performs a full reallocation
	// - Moves all data to a new object that satisfies the specified minimum capacities
	// - The new total capacity must also be at least as large as the old total capacity
	// - TODO: apply subvector capacity growth here?
	void reallocate(std::vector<size_type> capacities);

	// Calculates a new total capcity, including growth
	// - This must be at least as large as the current total capacity
	size_type calculateNewTotalCapacity(const std::vector<size_type>& capacities) const {
		size_type total_capacity = std::accumulate(capacities.begin(), capacities.end(), 0);
		return calculateNewTotalCapacity(total_capacity);
	}
	size_type calculateNewTotalCapacity(const size_type total_capacity) const {
		size_type min_total_capacity = std::max( this->capacity(), total_capacity );  // may not shrink
		return static_cast<size_type>( std::ceil(data_growth_factor_*min_total_capacity) );
	}

	// Calculates a new subvector capcity, including growth
	// - This must be at least as large as the current capacity
	size_type calculateNewSubvectorCapacity(const size_type old_capacity) const {
		return static_cast<size_type>( std::ceil(subvec_growth_factor_*old_capacity) ); 
	}

	// Total excess capacity beyond the amount allocated for the last subvector
	size_type getExcessCapacity() const;

	// Creates a new empty subvector at the end of the outer vector with the given capacity
	// - Takes advantage of excess capacity at the end of 'data_'
	// - Does not do any bounds checking!
	void allocateNewSubvectorAtEnd(const size_type capacity);

 private:
	Vector data_ = {};

	// Subvector arrangement
	std::vector<size_type> subvec_begins_     = {};  // offsets in 'data_' where each subvector begins
	std::vector<size_type> subvec_ends_       = {};  // one past the last element of each subvector
	std::vector<size_type> subvec_alloc_ends_ = {};  // one past where the subvector *capacities* end

	// TODO: Worth storing sizes explicitly for fast retrieval?
	//   subvec_sizes_[i] = subvec_ends_[i] - subvec_begins[i]
	//std::vector<size_type> subvec_sizes_  = {};

	static constexpr size_type initial_default_subvec_capacity_ = 10;
	size_type default_subvec_capacity_ = initial_default_subvec_capacity_;

	// These determine how quickly the outer and inner vectors grow upon reallocation
	// - Must be >1.0
	double data_growth_factor_   = 1.1;
	double subvec_growth_factor_ = 1.5;
};


template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::reference
VectorOfVectors<T,V>::operator()(const size_type i, const size_type j) noexcept
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
	FANCY_DEBUG_ASSERT( subvec_begins_[i]+j < data_.size(),
	                    "index " << subvec_begins_[i]+j << " out of bounds: size = " << data_.size() );
	return data_[subvec_begins_[i] + j];
} 

template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::const_reference
VectorOfVectors<T,V>::operator()(const size_type i, const size_type j) const noexcept
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
	FANCY_DEBUG_ASSERT( subvec_begins_[i]+j < data_.size(),
	                    "index " << subvec_begins_[i]+j << " out of bounds: size = " << data_.size() );
	return data_[subvec_begins_[i] + j];
}

template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::SubvectorIterator
VectorOfVectors<T,V>::begin(const size_type i) noexcept
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
	FANCY_DEBUG_ASSERT( subvec_begins_[i] < data_.size(),
	                    "index " << subvec_begins_[i] << " out of bounds: size = " << data_.size() );
	return &data_[subvec_begins_[i]];
}

template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::SubvectorIterator
VectorOfVectors<T,V>::end(const size_type i) noexcept 
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
	FANCY_DEBUG_ASSERT( subvec_ends_[i] < data_.size(),
	                    "index " << subvec_ends_[i] << " out of bounds: size = " << data_.size() );
	return &data_[subvec_ends_[i]];
}

template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::SubvectorConstIterator
VectorOfVectors<T,V>::begin (const size_type i) const noexcept
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
	FANCY_DEBUG_ASSERT( subvec_begins_[i] < data_.size(),
	                    "index " << subvec_begins_[i] << " out of bounds: size = " << data_.size() );
	return &data_[subvec_begins_[i]];
}

template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::SubvectorConstIterator
VectorOfVectors<T,V>::end(const size_type i) const noexcept
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
	FANCY_DEBUG_ASSERT( subvec_ends_[i] < data_.size(),
	                    "index " << subvec_ends_[i] << " out of bounds: size = " << data_.size() );
	return &data_[subvec_ends_[i]];
}

template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::SubvectorConstIterator
VectorOfVectors<T,V>::cbegin(const size_type i) const noexcept
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
	FANCY_DEBUG_ASSERT( subvec_begins_[i] < data_.size(),
	                    "index " << subvec_begins_[i] << " out of bounds: size = " << data_.size() );
	return &data_[subvec_begins_[i]];
}

template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::SubvectorConstIterator
VectorOfVectors<T,V>::cend(const size_type i) const noexcept
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );
	FANCY_DEBUG_ASSERT( subvec_ends_[i] < data_.size(),
	                    "index " << subvec_ends_[i] << " out of bounds: size = " << data_.size() );
	return &data_[subvec_ends_[i]];
}


template<typename T, typename V>
void VectorOfVectors<T,V>::checkInternalConsistency() const
{
	// Check member variable sizes
	size_type outer_size = this->size();
	FANCY_ASSERT( subvec_begins_.size()     == outer_size, "bad size" );
	FANCY_ASSERT( subvec_ends_.size()       == outer_size, "bad size" );
	FANCY_ASSERT( subvec_alloc_ends_.size() == outer_size, "bad size" );

	if ( outer_size > 0 ) {
		FANCY_ASSERT( subvec_alloc_ends_.back() <= capacity(), "bad capacity" );
	}

	// Check subvector setup
	for ( size_type i=0; i<outer_size; ++i ) {
		FANCY_ASSERT( subvec_begins_[i] <= subvec_ends_[i], "bad range" );
		FANCY_ASSERT( subvec_ends_[i] <= subvec_alloc_ends_[i], "bad range" );
		FANCY_ASSERT( size(i) <= capacity(i), "size is larger than capacity" );
		if ( i > 0 ) {
			FANCY_ASSERT( subvec_begins_[i] == subvec_alloc_ends_[i-1], "gap in storage" );
		}
	}
}


template<typename T, typename V>
std::ostream& VectorOfVectors<T,V>::dump(std::ostream& os) const
{
	os << "VectorOfVectors\n"
	   << "-----\n"
	   << "outer_size      = " << this->size()              << "\n"
	   << "total_capacity  = " << this->capacity()          << "\n"
	   << "num_elements    = " << this->count()             << "\n"
	   << "excess_capacity = " << this->getExcessCapacity() << "\n";

	size_type outer_size = this->size();
	for ( size_type i=0; i<outer_size; ++i ) {
		size_type inner_size = this->size(i);
		os << "[i=" << i << ", "
		   << "size=" << inner_size << ", capacity=" << this->capacity(i) << ", "
		   << "begin=" << subvec_begins_[i] << ", end=" << subvec_ends_[i] << ", alloc_end=" << subvec_alloc_ends_[i] << "]:";
		os << "  (";
		for ( size_type j=0; j<inner_size; ++j ) {
			if ( j > 0 ) {
				os << ", ";
			}
			os << (*this)(i,j);
		}
		os << ")\n";
	}
	return os;
}


template<typename T, typename V>
std::ostream& operator<<(std::ostream& os, const VectorOfVectors<T,V>& vec)
{
	using size_type = typename VectorOfVectors<T,V>::size_type;

	size_type outer_size = vec.size();
	for ( size_type i=0; i<outer_size; ++i ) {
		os << "[i=" << i << "]:  (";
		size_type inner_size = vec.size(i);
		for ( size_type j=0; j<inner_size; ++j ) {
			if ( j > 0 ) {
				os << ", ";
			}
			os << vec(i,j);
		}
		os << ")\n";
	}
	return os;
}


template<typename T, typename V>
typename VectorOfVectors<T,V>::size_type VectorOfVectors<T,V>::count() const
{
	size_type num = 0;
	for ( size_type i=0; i<size(); ++i ) {
		num += size(i);
	}
	return num;
}


template<typename T, typename V>
void VectorOfVectors<T,V>::push_back(const size_type i, const T& value)
{
	FANCY_DEBUG_ASSERT( i < size(), "index " << i << " out of bounds: size = " << size() );

	if ( subvec_ends_[i] < subvec_alloc_ends_[i] ) {
		data_[subvec_ends_[i]] = value;
		++(subvec_ends_[i]);
		return;
	}
	else {
		// Need to allocate more space
		size_type new_capacity = calculateNewSubvectorCapacity( capacity(i) );
		extendCapacity(i, new_capacity);

		data_[subvec_ends_[i]] = value;
		++(subvec_ends_[i]);
		return;
	}
#ifndef NDEBUG
	checkInternalConsistency();
#endif // ifndef NDEBUG
}


template<typename T, typename V>
void VectorOfVectors<T,V>::resize(const size_type i, const size_type new_size)
{
	FANCY_DEBUG_ASSERT(i < size(), "out of bounds (i=" << i << ", size=" << size() << ")" );

	size_type current_capacity = capacity(i);
	if ( new_size <= current_capacity ) {
		subvec_ends_[i] = subvec_begins_[i] + new_size;
		return;
	}
	else {
		size_type new_capacity = calculateNewSubvectorCapacity( new_size );
		extendCapacity(i, new_capacity);
		subvec_ends_[i] = subvec_begins_[i] + new_size;
		return;
	}
#ifndef NDEBUG
	checkInternalConsistency();
#endif // ifndef NDEBUG
}


template<typename T, typename V>
void VectorOfVectors<T,V>::extendCapacity(const size_type i, const size_type new_capacity)
{
	size_type old_capacity = capacity(i);
	if ( new_capacity <= old_capacity ) {
		return;  // nothing to do
	}

	size_type extra_capacity_needed = new_capacity - old_capacity;  // non-negative

	size_type outer_size = size();
	if ( i+1 == outer_size ) {
		// Extending the last subvector
		size_type excess_capacity = getExcessCapacity();
		if ( extra_capacity_needed > excess_capacity ) {
			// Allocate extra space
			data_.resize( data_growth_factor_*(data_.size() + new_capacity) );
		}
		subvec_alloc_ends_.back() += extra_capacity_needed;
		return;
	}
	else {
		// Insert capacity at the end of the subvector
		auto pos = data_.begin() + subvec_alloc_ends_[i];
		data_.insert( pos, extra_capacity_needed, T() );
		subvec_alloc_ends_[i] += extra_capacity_needed;
		for ( size_type k=i+1; k<outer_size; ++k ) {
			subvec_begins_[k]     += extra_capacity_needed;
			subvec_ends_[k]       += extra_capacity_needed;
			subvec_alloc_ends_[k] += extra_capacity_needed;
		}
	}
#ifndef NDEBUG
	checkInternalConsistency();
#endif // ifndef NDEBUG
}


template<typename T, typename V>
void VectorOfVectors<T,V>::resize(const size_type new_size)
{
	size_type old_size = this->size();
	if ( new_size == old_size ) {
		return;  // nothing else to do
	}
	else if ( new_size > old_size ) {
		size_type extra_capacity_needed_at_end = (new_size - old_size)*default_subvec_capacity_;
		if ( extra_capacity_needed_at_end > getExcessCapacity() ) {
			// Not enough space at the end: need to allocate more
			size_type new_capacity = calculateNewTotalCapacity(capacity() + extra_capacity_needed_at_end);
			data_.resize( new_capacity );
		}

		// Set up new subvectors at the end
		for ( size_type i=old_size; i<new_size; ++i ) {
			allocateNewSubvectorAtEnd( default_subvec_capacity_ );
		}
		return;
	}
	else { // new_size < old_size
		// Deallocate space
		subvec_begins_.resize(new_size);
		subvec_ends_.resize(new_size);
		subvec_alloc_ends_.resize(new_size);
	}
#ifndef NDEBUG
	checkInternalConsistency();
#endif // ifndef NDEBUG
}


template<typename T, typename V>
void VectorOfVectors<T,V>::resizeWithMinimumCapacities(const std::vector<size_type>& capacities)
{
	size_type new_size = capacities.size();
	size_type old_size = this->size();

	// Check whether an existing subvector is too small
	bool need_realloc = false;
	size_type i_max = std::min(old_size, new_size);
	for ( size_type i=0; i<i_max; ++i ) {
		if ( this->capacity(i) < capacities[i] ) {
			need_realloc = true;
			break;
		}
	}
	if ( need_realloc ) {
		this->reallocate(capacities);
		return;
	}

	if ( new_size == old_size ) {
		return;  // nothing else to do
	}
	else if ( new_size > old_size ) {
		// Check whether there is enough space at the end
		size_type extra_capacity_needed_at_end = 0;
		for ( size_type i=old_size; i<new_size; ++i ) {
			extra_capacity_needed_at_end += capacities[i];
		}
		if ( extra_capacity_needed_at_end <= getExcessCapacity() ) {
			for ( size_type i=old_size; i<new_size; ++i ) {
				allocateNewSubvectorAtEnd( capacities[i] );
			}
		}
		else {
			this->reallocate(capacities);
			return;
		}
	}
	else { // new_size < old_size
		// Deallocate space
		subvec_begins_.resize(new_size);
		subvec_ends_.resize(new_size);
		subvec_alloc_ends_.resize(new_size);
	}
#ifndef NDEBUG
	checkInternalConsistency();
#endif // ifndef NDEBUG
}


//-----------------------------//
//----- Protected Methods -----//
//-----------------------------//


template<typename T, typename V>
void VectorOfVectors<T,V>::assignEmptyWithCapacities(
	const std::vector<size_type>& capacities, const size_type total_capacity)
{
	const size_type size = capacities.size();
	if ( size == 0 ) {
		clear();
		return;
	}

	// Set total capacity
	size_type new_total_capacity = std::max( calculateNewTotalCapacity(capacities), total_capacity );
	data_.resize(new_total_capacity);

	// Set up offsets
	subvec_begins_.resize(size);       subvec_begins_[0]     = 0;
	subvec_ends_.resize(size);         subvec_ends_[0]       = 0;
	subvec_alloc_ends_.resize(size);   subvec_alloc_ends_[0] = capacities[0];
	for ( size_type i=1; i<size; ++i ) {
		subvec_begins_[i]     = subvec_alloc_ends_[i-1];
		subvec_ends_[i]       = subvec_begins_[i];
		subvec_alloc_ends_[i] = subvec_begins_[i] + capacities[i];
	}

#ifndef NDEBUG
	checkInternalConsistency();
#endif // ifndef NDEBUG
}


template<typename T, typename V>
void VectorOfVectors<T,V>::reallocate(std::vector<size_type> capacities)
{
	// Ensure that no data is lost in existing subvectors that are larger than
	// the minimum requested
	size_type new_size = capacities.size();
	size_type old_size = this->size();
	size_type i_max = std::min(old_size, new_size);
	for ( size_type i=0; i<i_max; ++i ) {
		capacities[i] = std::max( capacities[i], this->capacity(i) );
	}

	// Move data to a new object of the appropriate dimensions
	VectorOfVectors new_obj(capacities, this->capacity());
	for ( size_type i=0; i<i_max; ++i ) {
		size_type first = this->subvec_begins_[i];
		size_type last  = this->subvec_ends_[i];
		size_type k_new = new_obj.subvec_begins_[i];
		for ( size_type k = first; k != last; ++k, ++k_new ) {
			std::swap( new_obj.data[k_new], this->data_[k] );
		}
		new_obj.subvec_ends_[i] = new_obj.subvec_begins_[i] + this->size(i);
	}

#ifndef NDEBUG
	checkInternalConsistency();
#endif // ifndef NDEBUG
}


template<typename T, typename V>
inline
typename VectorOfVectors<T,V>::size_type
VectorOfVectors<T,V>::getExcessCapacity() const
{
	size_type excess_capacity = data_.size();
	if ( this->size() > 0 ) {
		excess_capacity -= subvec_alloc_ends_.back();
	}
	return excess_capacity;
}


template<typename T, typename V>
void VectorOfVectors<T,V>::allocateNewSubvectorAtEnd(const size_type capacity)
{
	FANCY_DEBUG_ASSERT( capacity <= getExcessCapacity(), "improper use (not enough space)" );

	size_type first = 0;
	if ( size() > 0 ) {
		first = subvec_alloc_ends_.back();
	}
	subvec_begins_.push_back( first );
	subvec_ends_.push_back( first );
	subvec_alloc_ends_.push_back( first + capacity );

#ifndef NDEBUG
	checkInternalConsistency();
#endif // ifndef NDEBUG
}

#endif //ifndef VECTOR_OF_VECTORS_H
