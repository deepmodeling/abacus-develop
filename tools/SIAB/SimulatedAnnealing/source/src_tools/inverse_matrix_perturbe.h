#ifndef INVERSE_MATRIX_PERTURBE_H
#define INVERSE_MATRIX_PERTURBE_H

#include "complexmatrix.h"

// know A.inverse(),B
// cal (A+-B).inverse()
// (A-B)^{-1} = A^{-1} + A^{-1}*B*A^{-1} + A^{-1}*B*A^{-1}*B*A^{-1} + ...

class Inverse_Matrix_Perturbe
{
    public:
	void init(const int &dim);
	ComplexMatrix &cal_inverse();
	ComplexMatrix &get_I() { return I; }
	
	private:
    const ComplexMatrix *A_inverse{nullptr};
	const ComplexMatrix *B{nullptr};
	ComplexMatrix I;
	size_t iterate_num;
	int plus_minus_state{0};		// 1 plus, 2 minus
	
	public:
	const ComplexMatrix *get_A_inverse() { return A_inverse; }
	Inverse_Matrix_Perturbe &set_A_inverse(ComplexMatrix *A_inverse_in) { A_inverse = A_inverse_in; return *this;}
	const ComplexMatrix *get_B() { return B; }
	Inverse_Matrix_Perturbe &set_B(ComplexMatrix *B_in) { B = B_in; return *this;}
	const size_t &get_iterate_num() { return iterate_num; }
	Inverse_Matrix_Perturbe &set_iterate_num(size_t iterate_num_in) { iterate_num = iterate_num_in; return *this;}
	const int &get_plus_minus_state() { return plus_minus_state; }
	Inverse_Matrix_Perturbe &set_plus_minus_state(int plus_minus_state_in) { plus_minus_state = plus_minus_state_in; return *this;}
};

#endif