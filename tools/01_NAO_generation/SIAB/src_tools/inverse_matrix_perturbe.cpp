#include "inverse_matrix_perturbe.h"
#include<iostream>
ComplexMatrix &Inverse_Matrix_Perturbe::cal_inverse()
{
//std::cout<<"@";
	this->I = *this->A_inverse;
	const ComplexMatrix b_a_inverse(*this->B * *this->A_inverse);
	ComplexMatrix add_item(*this->A_inverse);
	for(size_t i(1); i<=iterate_num; ++i)
	{
		add_item = add_item * b_a_inverse;
		if(plus_minus_state==2) { this->I += add_item; }
		else if(plus_minus_state==1)
		{
			if((i&1) == 0) { this->I += add_item; }
			else { this->I -= add_item; }
		}	
		else { throw invalid_argument("set plus_minus_state"); }
	}
	return this->I;
}

void Inverse_Matrix_Perturbe::init(const int &dim)
{
    this->I.create(dim, dim);
    return;
}