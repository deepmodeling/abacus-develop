#include "Inverse_Matrix_S.h"

void Inverse_Matrix_S::copy_matrix_U2L(ComplexMatrix &matrix)
{
	for( size_t iw1(0); iw1<this->dim; ++iw1)
	{
		for( size_t iw2(0); iw2<iw1; ++iw2)
		{
			matrix(iw1,iw2) = conj(matrix(iw2,iw1));
		}
	}	
}

void Inverse_Matrix_S::zero_matrix_L(ComplexMatrix &matrix)
{
	for( size_t iw1(0); iw1<this->dim; ++iw1)
	{
		for( size_t iw2(0); iw2<iw1; ++iw2)
		{
			matrix(iw1,iw2) = 0;
		}
	}	
}

void Inverse_Matrix_S::update()
{
	inverse_old = &get_inverse();
	loop_count = (loop_count+1)%exact_loop_num;
}