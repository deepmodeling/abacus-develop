#ifndef INVERSE_MBTRIX_UNITTEST_H
#define INVERSE_MBTRIX_UNITTEST_H

#include "../../src_tools/inverse_matrix.h"
#include "../../src_tools/complexmatrix_inline.h"

#include "common_test.h"

void Inverse_Matrix_Unittest()
{
	const int dim(3);
	ComplexMatrix B(dim,dim);
	B(0,0)=1;	B(0,1)=2; B(0,2)=4;
	B(1,1)=13;	B(1,2)=23;
	B(2,2)=77;
		
	Inverse_Matrix inverse;
	inverse.init(dim);
	
	inverse.using_zpotrf(B);
	cout_matrix(inverse.A);
}

#endif

/* output:
(1.4567901235,0.0000000000)     (-0.1913580247,-0.0000000000)   (-0.0185185185,-0.0000000000)
(0.0000000000,0.0000000000)     (0.1882716049,0.0000000000)     (-0.0462962963,-0.0000000000)
(0.0000000000,0.0000000000)     (0.0000000000,0.0000000000)     (0.0277777778,0.0000000000)
*/