#ifndef MATRIX_UNITTEST_H
#define MATRIX_UNITTEST_H

#include "module_base/matrix.h"

static void matrix_mutiply_test()
{
	const size_t dim1=2,dim2=3,dim3=4;
	
	// C = A * B
	matrix A(dim1,dim2), B(dim2,dim3), C(dim1,dim3);

	A(0,0)=1; A(0,1)=3;A(0,2)=4; 
	A(1,0)=5; A(1,1)=6;A(1,2)=8;
	
	B(0,0)=2; B(0,1)=3;B(0,2)=5;B(0,3)=4; 
	B(1,0)=5; B(1,1)=6;B(1,2)=8;B(1,3)=6;
	B(2,0)=7; B(2,1)=6;B(2,2)=5;B(2,3)=4;

	C = A * B;

	cout<<C<<endl;
}

/* output:
	45	45	49	38
	96	99	113	88
*/

#endif