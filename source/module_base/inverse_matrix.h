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
	int dim=0;
	double *e=nullptr;
	int lwork=0;
	std::complex<double> *work2=nullptr;
	double* rwork=nullptr;
	int info=0;
	bool allocate=false; //mohan add 2012-04-02

	ModuleBase::ComplexMatrix EA;
};

/**
 * @brief computes the inverse of a dim*dim real matrix "in" using LAPACK routines
 * "out" contains the inverse on output; "in" is unchanged
 *
 * @param dim [in] dimension of the matrix
 * @param in [in] input matrix
 * @param out [out] output matrix
 */
void Inverse_Matrix_Real(const int dim, const double* in, double* out);
}
#endif
