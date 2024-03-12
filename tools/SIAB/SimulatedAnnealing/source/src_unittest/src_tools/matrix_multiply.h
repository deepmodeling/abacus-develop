#include "../../src_tools/complexmatrix.h"
#include "../../src_tools/complexmatrix_inline.h"

#include "common_test.h"

inline void matrix_multiply_test_1()
{
	const size_t dim1(2),dim2(3),dim3(4);
	
	// C = A * B
	ComplexMatrix A(dim1,dim2), B(dim2,dim3), C(dim1,dim3);

	A(0,0)=1; A(0,1)=3;A(0,2)=4; 
	A(1,0)=5; A(1,1)=6;A(1,2)=8;
	
	B(0,0)=2; B(0,1)=3;B(0,2)=5;B(0,3)=4; 
	B(1,0)=5; B(1,1)=6;B(1,2)=8;B(1,3)=6;
	B(2,0)=7; B(2,1)=6;B(2,2)=5;B(2,3)=4;

	C = A * B;

	cout_matrix(C);
}

inline void matrix_multiply_test_2()
{
	int M(2),N(2);
	ComplexMatrix A(M,M), B(M,N), C(M,N);
	
	A(0,0)=1; A(0,1)=2;
	A(1,0)=3; A(1,1)=3;
	
	B(0,0)=4; B(0,1)=1;
	B(1,0)=4; B(1,1)=2;

	C = A*B;
	cout_matrix(C);
	
	C = multiply_special('H','R','U',A,B);
	cout_matrix(C);
	
	C = multiply_special('H','L','L',A,B);
	cout_matrix(C);	
}

/* output:
(12.0000000000,0.0000000000)    (5.0000000000,0.0000000000)
(24.0000000000,0.0000000000)    (9.0000000000,0.0000000000)

(6.0000000000,0.0000000000)     (11.0000000000,0.0000000000)
(8.0000000000,0.0000000000)     (14.0000000000,0.0000000000)

(16.0000000000,0.0000000000)    (7.0000000000,0.0000000000)
(24.0000000000,0.0000000000)    (9.0000000000,0.0000000000)
*/

