#ifndef INVERSE_MATRIX_ITERATE_UNITTEST_H
#define INVERSE_MATRIX_ITERATE_UNITTEST_H

#include "../../src_tools/inverse_matrix_iterate.h"
#include "../../src_tools/inverse_matrix_iterate_inline.h"
#include "common_test.h"

void Inverse_Matrix_Iterate_Unittest()
{
	const size_t dim(2);
	
	ComplexMatrix A(dim,dim), B(dim,dim), A_I(dim,dim), A_B(dim,dim);
	A(0,0)=1; A(0,1)=3; A(1,0)=5; A(1,1)=6;
	B(0,0)=0.2; B(0,1)=0.1; B(1,0)=0.4; B(1,1)=0.3;
	A_I(0,0)=-0.66666667; A_I(0,1)=0.33333333; A_I(1,0)=0.55555556; A_I(1,1)=-0.11111111;
	A_B = A-B;

	Inverse_Matrix_Iterate inverse;	
	inverse.init(dim);
	inverse.set_X0(A_I).set_A(A_B).set_tolerate_error(1E-10).set_fluctrate_error(1E-4);
	
	for(size_t iterate_num(0); iterate_num<=10; ++iterate_num)
	{
		cout<<iterate_num<<endl;
		const ComplexMatrix ii(inverse.set_iterate_num(iterate_num).cal_inverse());
		cout_matrix(ii);
	}
}

#endif

/* output
0
(-0.666667,0)   (0.333333,0)
(0.555556,0)    (-0.111111,0)

1
(-0.648148,0)   (0.32963,0)
(0.523457,0)    (-0.091358,0)

2
(-0.649201,0)   (0.330294,0)
(0.523918,0)    (-0.0911178,0)

3
(-0.649203,0)   (0.330296,0)
(0.523918,0)    (-0.0911162,0)

4
(-0.649203,0)   (0.330296,0)
(0.523918,0)    (-0.0911162,0)

5
(-0.649203,0)   (0.330296,0)
(0.523918,0)    (-0.0911162,0)

6
(-0.649203,0)   (0.330296,0)
(0.523918,0)    (-0.0911162,0)

7
(-0.649203,0)   (0.330296,0)
(0.523918,0)    (-0.0911162,0)

8
(-0.649203,0)   (0.330296,0)
(0.523918,0)    (-0.0911162,0)

9
(-0.649203,0)   (0.330296,0)
(0.523918,0)    (-0.0911162,0)

10
(-0.649203,0)   (0.330296,0)
(0.523918,0)    (-0.0911162,0)
*/