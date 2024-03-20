#ifndef INVERSE_MATRIX_ITERATE_H
#define INVERSE_MATRIX_ITERATE_H

#include "complexmatrix.h"

#include<cstddef>
#include<exception>
#include<limits>

// input A, X_0
// cal A_I
// iterate: X_{n+1} = X_n (2I - A*X_n )
 
class Inverse_Matrix_Iterate
{
	public:
	inline void init(size_t dim_in);
	inline Inverse_Matrix_Iterate &set_X0(const ComplexMatrix &X0_in);
	Inverse_Matrix_Iterate &set_A(ComplexMatrix &A_in) { this->A_ptr=&A_in; return *this;};
	inline const ComplexMatrix &cal_inverse();
	const ComplexMatrix &get_inverse() { return *this->inverse_ptr; }

	private:
	inline void iterate(const ComplexMatrix &old_X, ComplexMatrix &new_X);	
	
	private:
	ComplexMatrix *A_ptr{nullptr};
	ComplexMatrix X_store_1;
	ComplexMatrix X_store_2;
	const ComplexMatrix *X0_ptr{nullptr};
	const ComplexMatrix *inverse_ptr{nullptr};
	ComplexMatrix identity;
	ComplexMatrix matrix_tmp_1;
	ComplexMatrix matrix_tmp_2;
	ComplexMatrix X0_tmp;
	size_t iterate_num{20};
	double tolerate_error2{1E-20};
	double fluctrate_error2{0.0};
	double last_error2{std::numeric_limits<double>::max()};
	struct Matrix_Info{ size_t dim; char symmetry{NULL}; char uplo{NULL};} matrix_info;

	public:
	Inverse_Matrix_Iterate &set_iterate_num(size_t iterate_num_in) { this->iterate_num=iterate_num_in; return *this;}
	Inverse_Matrix_Iterate &set_tolerate_error(double tolerate_error_in) { this->tolerate_error2=tolerate_error_in*tolerate_error_in; return *this;}
	Inverse_Matrix_Iterate &set_fluctrate_error(double fluctrate_error_in) { this->fluctrate_error2=fluctrate_error_in*fluctrate_error_in; return *this;}
	Inverse_Matrix_Iterate &set_symmetric_info(char symmetry_in, char uplo_in) { this->matrix_info.symmetry = symmetry_in; this->matrix_info.uplo = uplo_in; return *this;}
	
	private:
	class Little_Error: public exception {};
	public:
	class Unstable: public exception {};
	private:
};

#endif