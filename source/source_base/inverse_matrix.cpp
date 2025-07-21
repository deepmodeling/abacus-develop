#include "inverse_matrix.h"
#include "tool_quit.h"
#include "module_external/lapack_connector.h"
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
		allocate=false;
	}	
}

void Inverse_Matrix_Complex::init(const int &dim_in)
{
//	GlobalV::ofs_running << " allocate=" << allocate << std::endl;
	if(allocate)
	{
		delete[] e; //mohan fix bug 2012-04-02
		allocate=false;
	}	

	this->dim = dim_in;

	assert(dim>0);
	this->e = new double[dim];

	assert(lwork>0);

	assert(3*dim-2>0);
	this->A.create(dim, dim);
	this->EA.create(dim, dim);

	this->allocate = true;

	return;
}


void Inverse_Matrix_Complex::using_zheev( const ModuleBase::ComplexMatrix &Sin, ModuleBase::ComplexMatrix &Sout)
{
	ModuleBase::timer::tick("Inverse","using_zheev");
	this->A = Sin;

    LapackConnector::heev(LapackConnector::RowMajor, 'V', 'U', dim, this->A.c, dim, e);
	
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

void Inverse_Matrix_Real(const int dim, const double* in, double* out)
{
    int lda = dim;
    std::vector<int> ipiv(dim);

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            out[i * dim + j] = in[i * dim + j];
        }
    }

    LapackConnector::getrf(LapackConnector::ColMajor, dim, dim, out, lda, ipiv.data());
    LapackConnector::getri(LapackConnector::ColMajor, dim, out, lda, ipiv.data());
}
}