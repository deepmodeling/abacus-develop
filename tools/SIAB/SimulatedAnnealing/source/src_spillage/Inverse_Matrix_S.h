#ifndef INVERSE_MATRIX_S_H
#define INVERSE_MATRIX_S_H

#include "../src_tools/inverse_matrix.h"
#include "../src_tools/inverse_matrix_iterate.h"

class Inverse_Matrix_S
{
public:
	Inverse_Matrix inverse;
	Inverse_Matrix_Iterate inverse_iterate;	

	void init( int dim_in);
	const ComplexMatrix & cal_inverse(ComplexMatrix &S);
	const ComplexMatrix & get_inverse() { return (0==loop_count) ? inverse.A : inverse_iterate.get_inverse(); }
	inline void update();		
	void reset_loop() { loop_count = 0; }

private:
	size_t dim{0};
	size_t loop_count{0};
	size_t exact_loop_num{1};
	const ComplexMatrix *inverse_old{nullptr};
	
	inline void copy_matrix_U2L(ComplexMatrix &matrix);
	inline void zero_matrix_L(ComplexMatrix &matrix);
	
public:
	size_t get_exact_loop_num() { return this->exact_loop_num; }
	size_t get_loop_count() { return this->loop_count; }
	Inverse_Matrix_S & set_exact_loop_num(size_t exact_loop_num_in) { this->exact_loop_num = exact_loop_num_in; return *this;}
	size_t get_dim() { return this->dim; }
};

#endif
