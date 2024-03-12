#include "Inverse_Matrix_S.h"
#include "Inverse_Matrix_S_inline.h"
#include "../src_tools/inverse_matrix_iterate_inline.h"

void Inverse_Matrix_S::init( int dim_in)
{
	this->dim = dim_in;
	this->inverse.init(this->dim);
	this->inverse_iterate.init(this->dim);
	this->inverse_iterate.set_iterate_num(20);		// test, undetermined
	this->inverse_iterate.set_tolerate_error(1e-10);		// test, undetermined
	this->inverse_iterate.set_symmetric_info('H', 'U');
	this->set_exact_loop_num(1);					// test, undetermined
	return;	
}

const ComplexMatrix & Inverse_Matrix_S::cal_inverse(ComplexMatrix &S)
{
	if(0==loop_count)
	{
		inverse.using_zpotrf(S);				
	}
	else
	{
		inverse_iterate.set_X0(*this->inverse_old)
		               .set_A(S);				
		try
		{
			if (1==loop_count) { copy_matrix_U2L(const_cast<ComplexMatrix &>(*this->inverse_old)); }			
//			copy_matrix_U2L(const_cast<ComplexMatrix &>(*this->inverse_old));
//			copy_matrix_U2L(S);

			inverse_iterate.cal_inverse();
		
//			zero_matrix_L(const_cast<ComplexMatrix &>(*this->inverse_old));
//			zero_matrix_L(S);
		}
		catch(const Inverse_Matrix_Iterate::Unstable &e)
		{			
//			zero_matrix_L(const_cast<ComplexMatrix &>(*this->inverse_old));
//			zero_matrix_L(S);
			
			loop_count = 0;
			inverse.using_zpotrf(S);
		}
	}
	return this->get_inverse();
}