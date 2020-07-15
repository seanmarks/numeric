// AUTHOR: Sean M. Marks (https://github.com/seanmarks)

#ifndef NUMERIC_ARRAY_H
#define NUMERIC_ARRAY_H

#include <cstddef>
#include <iostream>
#include <iterator>

namespace numeric {

// Fixed-size array
template<typename T, std::size_t N>
class Array
{
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


	//----- Constructors -----//

	Array() {}

	Array(const_reference value) {
		fill(value);
	}

	template<typename U>
	Array(const Array<U,N>& other) {
		*this = other;
	}


	//----- Assignment -----//

	// Array
	template<typename U>
	Array& operator=(const Array<U,N>& other) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] = other[i];
		}
		return *this;
	}

	// Scalar
	template<typename U>
	Array& operator=(const U& value) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] = value;
		}
		return *this;
	}


	//----- Data access -----//

	reference       operator[](const size_type i)       { return data_[i]; }
	const_reference operator[](const size_type i) const { return data_[i]; }
	reference       operator()(const size_type i)       { return data_[i]; }
	const_reference operator()(const size_type i) const { return data_[i]; }

	iterator       begin()        { return &data_[0]; }
	const_iterator cbegin() const { return &data_[0]; }
	iterator       end()          { return &data_[N]; }
	const_iterator cend()   const { return &data_[N]; } 


	//----- Managing contents -----//

	void fill(const_reference value) {
		for ( reference v : data_ ) {
			v = value;
		}
	}

	// TODO necessary?
	void swap(Array<T,N>& other) noexcept {
		using std::swap;
		for ( unsigned i=0; i<N; ++i ) {
			swap(data_[i], other[i]);
		}
	}
	

	//----- Element-wise, in-place arithmetic -----//

	template<typename U>
	Array<T,N>& operator+=(const Array<U,N>& other) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] += other[i];
		}
		return *this;
	}

	template<typename U>
	Array<T,N>& operator-=(const Array<U,N>& other) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] -= other[i];
		}
		return *this;
	}

	template<typename U>
	Array<T,N>& operator*=(const Array<U,N>& other) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] *= other[i];
		}
		return *this;
	}

	template<typename U>
	Array<T,N>& operator/=(const Array<U,N>& other) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] /= other[i];
		}
		return *this;
	}


	//----- Scalar in-place arithmetic -----//

	template<typename U>
	Array<T,N>& operator+=(const U& value) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] += value;
		}
		return *this;
	}

	template<typename U>
	Array<T,N>& operator-=(const U& value) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] -= value;
		}
		return *this;
	}

	template<typename U>
	Array<T,N>& operator*=(const U& value) {
		for ( unsigned i=0; i<N; ++i ) {
			data_[i] *= value;
		}
		return *this;
	}

	template<typename U>
	Array<T,N>& operator/=(const U& value) {
		U inv_value = 1.0/value;
		*this *= inv_value;
		return *this;
	}


	//----- Output/Debugging -----//

	//template<typename U, std::size_t M>
	//friend std::ostream& operator<<(std::ostream& os, const Array<U,M>& array);

 private:
	T data_[N];
};


template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const Array<T,N>& array)
{
	os << "[";
	for ( unsigned i=0; i<N; ++i ) {
		if ( i > 0 ) { os << ", "; }
		os << array[i];
	}
	os << "]";
	return os;
}


template<typename T, typename U, std::size_t N>
auto operator+(const Array<T,N>& left, const Array<U,N>& right) -> Array<decltype(left[0] + right[0]),N>
{
	Array<decltype(left[0] + right[0]),N> result(left);
	result += right;
	return result;
}

/*
template<typename T, std::size_t N>
Array<T,N> operator+(const Array<T,N>& left, const Array<T,N>& right)
{
	Array<T,N> result(left);
	result += right;
	return result;
}
*/

} // end namespace numeric

#endif // ifndef NUMERIC_ARRAY_H
