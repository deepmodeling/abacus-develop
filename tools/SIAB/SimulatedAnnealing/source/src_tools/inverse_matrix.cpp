#include "inverse_matrix.h"
#include "../src_spillage/Inverse_Matrix_S.h"
#include "complexmatrix_inline.h"

Inverse_Matrix::Inverse_Matrix(){}
Inverse_Matrix::~Inverse_Matrix(){}
void Inverse_Matrix::init(const int &dim_in)
{
    this->dim = dim_in;
    this->e = new double[dim];
    this->lwork = 2*dim;
    this->work2 = new complex<double>[lwork];
    this->rwork = new double[3*dim-2];
    this->info = 0;
    this->A.create(dim, dim);
    this->EA.create(dim, dim);
    return;
}

void Inverse_Matrix::using_zheev( const ComplexMatrix &Sin, ComplexMatrix &Sout)
{
    timer::tick("Inverse","using_zheev");
    this->A = Sin;

    //LapackConnector::zhegv( 1, 'V', 'U', nwan , A ,  nwan , B , nwan , e, work2 , 80 , rwork , &info );
    LapackConnector::zheev('V', 'U', dim, this->A, dim, e, work2, lwork, rwork, &info);

    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            EA(i,j)= conj( this->A(j,i) ) / e[i] ;
        }
    }

    Sout = this->A * this->EA;
    timer::tick("Inverse","using_zheev");
    return;
}

void Inverse_Matrix::using_zpotrf( const ComplexMatrix &Sin)
{
//  timer::tick("Inverse","using_zpotrf");

    for(int i=0; i<dim; i++)
    {
        for(int j=i; j<dim; j++)
        {
            A(i,j) = Sin(i,j);
        }
    }

    LapackConnector::zpotrf('U',dim,A,dim,info);	
//  cout << "\n info=" << info << endl;
    LapackConnector::zpotri('U',dim,A,dim,info);
//  cout << "\n info=" << info << endl;
    //PRINTCM("A_zpotrf",A);

//  timer::tick("Inverse","using_zpotrf");
    return;
}

