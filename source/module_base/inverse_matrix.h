#ifndef INVERSE_MATRIX_H
#define INVERSE_MATRIX_H

#include "global_function.h"
#include "global_variable.h"
#include "complexmatrix.h"
namespace ModuleBase
{

class Inverse_Matrix_Complex
{
	public:

	Inverse_Matrix_Complex();
	~Inverse_Matrix_Complex();
	
	ModuleBase::ComplexMatrix A;

	void using_zheev(const ModuleBase::ComplexMatrix &in, ModuleBase::ComplexMatrix &out);
	void init( const int &dim_in);

	private:
	int dim;
	double *e;
	int lwork;
	std::complex<double> *work2;
	double* rwork;
	int info;
	bool allocate; //mohan add 2012-04-02

	ModuleBase::ComplexMatrix EA;
};


}
#endif
