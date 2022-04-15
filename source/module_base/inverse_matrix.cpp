#include "inverse_matrix.h"
#include "tool_quit.h"
#include "lapack_connector.h"
#include "timer.h"
#include "complexmatrix.h"

namespace ModuleBase
{

Inverse_Matrix_Complex::Inverse_Matrix_Complex()
{
	allocate=false;	
}
Inverse_Matrix_Complex::~Inverse_Matrix_Complex()
{
	if(allocate)
	{
		delete[] e; //mohan fix bug 2012-04-02
		delete[] work2;
		delete[] rwork;
		allocate=false;
	}	
}

void Inverse_Matrix_Complex::init(const int &dim_in)
{
//	GlobalV::ofs_running << " allocate=" << allocate << std::endl;
	if(allocate)
	{
		delete[] e; //mohan fix bug 2012-04-02
		delete[] work2;
		delete[] rwork;
		allocate=false;
	}	

	this->dim = dim_in;

	assert(dim>0);
	this->e = new double[dim];
	this->lwork = 2*dim;

	assert(lwork>0);
	this->work2 = new std::complex<double>[lwork];

	assert(3*dim-2>0);
	this->rwork = new double[3*dim-2];
	this->info = 0;
	this->A.create(dim, dim);
	this->EA.create(dim, dim);

	this->allocate = true;

	return;
}


void Inverse_Matrix_Complex::using_zheev( const ModuleBase::ComplexMatrix &Sin, ModuleBase::ComplexMatrix &Sout)
{
	ModuleBase::timer::tick("Inverse","using_zheev");
	this->A = Sin;

    LapackConnector::zheev('V', 'U', dim, this->A, dim, e, work2, lwork, rwork, &info);
	
	for(int i=0; i<dim; i++)
	{
		for(int j=0; j<dim; j++)
		{
			EA(i,j)= conj( this->A(j,i) ) / e[i] ;
		}
	}

    Sout = this->A * this->EA;
	ModuleBase::timer::tick("Inverse","using_zheev");
    return;
}

}