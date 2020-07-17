// Matrix: A flexible 2-dimensional array
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// - Implemented using a 1-dimensional array (default: std::vector)

#ifndef NUMERIC_MATRIX_H
#define NUMERIC_MATRIX_H

#include <cmath>
#include <type_traits>
#include <vector>

#include "Assert.h"
#include "OpenMP.h"

namespace numeric {

template<typename T, typename V> class Matrix;

// Forward-declare friend functions
template<typename T, typename V> Matrix<T,V> log(const Matrix<T,V>& other);


template<typename T, class Vector = std::vector<T>>
class Matrix
{
 public:
	static_assert(std::is_same<T, typename Vector::value_type>::value, "type mismatch");

	static constexpr int X_DIM = 0;
	static constexpr int Y_DIM = 1;
	static constexpr int N_DIM = 2;

	using Int2 = std::array<int,N_DIM>;


	//----- Constructors -----//

	Matrix() {
		this->clear();
	}

	Matrix(const Int2& shape): 
		shape_(shape)
	{
		setShape(shape);
	}

	Matrix(const Int2& shape, const T& value): 
		shape_(shape)
	{
		setShape(shape);
		assign(value);
	}


	//----- Size Management -----//

	// Set size
	void setShape(const Int2& shape) {
		// Allocate memory
		int len = shape[X_DIM]*shape[Y_DIM];
		data_.resize(len);

		shape_ = shape;
	}
	void setShape(const int num_rows, const int num_cols) {
		setShape( {{num_rows, num_cols}} );
	}
	void resize(const Int2& shape) {
		setShape(shape);
	}
	void resize(const int num_rows, const int num_cols) {
		setShape( {{num_rows, num_cols}} );
	}

	// Get size(s)
	const Int2& getShape() const {
		return shape_;
	}


	//----- Data Management -----//

	// Assign all entries to a particular value
	void assign(const T& value) {
		data_.assign( data_.size(), value );
	}

	void assign(const Int2& shape, const T& value) {
		this->setShape(shape);  // also resize the array
		this->assign(value);
	}

	void zero() {
		this->assign(0.0);
	}

	// Clear all contents
	void clear() {
		data_.clear();
		shape_.fill(0);
	}


	//----- Copying -----//

	// Copy assignment (from an array of a different type)
	template<typename U>
	Matrix<T>& operator=(const Matrix<U>& other)
	{
		this->setShape( other.getShape() );

		int len = data_.size();
		const auto& other_data = other.data();
		for ( int i=0; i<len; ++i ) {
			this->data_[i] = other_data[i];
		}

		return *this;
	}
	/*
	// Copy from an array of a different type
	template<typename U>
	Matrix(const Matrix<U>& other): 
		shape_(other.getShape())
	{
		setShape(shape);

		int len = data_.size();
		for ( int i=0; i<len; ++i ) {
			this->data_[i] = other.data_[i];
		}
	}
	*/

	//----- Data access -----//

	// Individual elements
	T& operator()(const int i, const int j) {
		FANCY_DEBUG_ASSERT( getLinearIndex(i,j) < static_cast<int>(data_.size()),
		                    "indices (" << i << "," << j << ") are out of bounds "
		                    << "(" shape_[X_DIM] << "," << shape_[Y_DIM] << ")" );
		return data_[ getLinearIndex(i,j) ];
	}
	const T& operator()(const int i, const int j) const {
		// TODO DEBUG MODE: check bounds
		return data_[ getLinearIndex(i,j) ];
	}
	T& operator()(const Int2& indices) {
		// TODO DEBUG MODE: check bounds
		return data_[ getLinearIndex(indices) ];
	}
	const T& operator()(const Int2& indices) const {
		// TODO DEBUG MODE: check bounds
		return data_[ getLinearIndex(indices) ];
	}

	// Underlying 1D array (use with caution!)
	Vector& data() noexcept {
		return data_;
	}
	const Vector& data() const noexcept {
		return data_;
	}


	//----- Array Properties -----//

	// Sum of all elements
	T sum() const {
		int len = this->data_.size();
		double s = 0.0;

		#pragma omp parallel for \
			default(shared) schedule(static,10) reduction(+:s)
		for ( int i=0; i<len; ++i ) {
			s += data_[i];
		}

		return s;
	}


	//----- Arithmetic Operators -----//

	// Between two arrays
	Matrix& operator+=(const Matrix& other) {
		// TODO DEBUG MODE: check shape
		int len = this->data_.size();
		#pragma omp parallel for \
			default(shared) schedule(static,10)
		for ( int i=0; i<len; ++i ) {
			this->data_[i] += other.data_[i];
		}
		return *this;
	}

	// Array and scalar
	template<typename U>
	Matrix& operator+=(const U& value) {
		int len = this->data_.size();
		#pragma omp parallel for \
			default(shared) schedule(static,10)
		for ( int i=0; i<len; ++i ) {
			this->data_[i] += value;
		}
		return *this;
	}
	template<typename U>
	Matrix& operator*=(const U& value) {
		int len = this->data_.size();
		#pragma omp parallel for \
			default(shared) schedule(static,10)
		for ( int i=0; i<len; ++i ) {
			this->data_[i] *= value;
		}
		return *this;
	}


	//----- Friend Functions -----//

	template<typename U, typename V>
	friend Matrix<U,V> log(const Matrix<U,V>& other);


 protected:
	// Map from 2D indices to 1D index of underlying array
	int getLinearIndex(const int i, const int j) const noexcept {
		return i*shape_[Y_DIM] + j;
	}
	int getLinearIndex(const Int2& indices) const noexcept {
		return getLinearIndex( indices[X_DIM], indices[Y_DIM] );
	}

 private:
	Vector data_;    // underlying 1D array
	Int2   shape_;   // array dimensions
};


//----- Friend Functions -----//

template<typename T, typename V>
Matrix<T,V> log(const Matrix<T,V>& other)
{
	Matrix<T,V> arr_out( other.getShape() );
	int len = arr_out.data_.size();

	#pragma omp parallel for \
		default(shared) schedule(static,10)
	for ( int i=0; i<len; ++i ) {
		arr_out.data_[i] = std::log( other.data_[i] );
	}

	return arr_out;
}

} // end namespace numeric

#endif // ifndef NUMERIC_MATRIX_H
