#include "../inverse_matrix.h"
#include "gtest/gtest.h"
#include <iostream>

/************************************************
 *  unit test of class Inverse_Matrix_Complex
 ***********************************************/

/**
 * - Tested Functions:
 *   - Inverse
 *     - use Inverse_Matrix_Complex to inverse a Hermite matrix
 *     - functions: init and using_zheev
 *
 */

//a mock function of WARNING_QUIT, to avoid the uncorrected call by matrix.cpp at line 37.
namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {return ;}
}

TEST(InverseMatrixComplexTest,Inverse)
{
	int dim = 10;
	ModuleBase::ComplexMatrix B(dim,dim);
	ModuleBase::ComplexMatrix C(dim,dim);
	ModuleBase::ComplexMatrix D(dim,dim);
	double a;
	double b;
	double c;
	// construct a Hermite matrix
	for(int j=0;j<dim;j++)
	{
		for(int i=0;i<=j;i++)
		{
			if (i==j)
			{
				c = std::rand();
				B(i,j) = std::complex<double> (c,0.0);
			}
			else
			{
				a = std::rand();
				b = std::rand();
				B(i,j) = std::complex<double> (a,b);
				B(j,i) = conj(B(i,j));
			}
		}
	}
	ModuleBase::Inverse_Matrix_Complex IMC;
	IMC.init(dim);
	IMC.using_zheev(B,C);
	D = B * C;
	for(int i=0;i<dim;i++)
	{
		EXPECT_NEAR(D(i,i).real(),1.0,1e-14);
		EXPECT_NEAR(D(i,i).imag(),0.0,1e-14);
		//std::cout << D(i,i).real() << " " << D(i,i).imag()  << std::endl;
	}
}
