#ifndef INVERSE_MATRIX_H
#define INVERSE_MATRIX_H

#include "../src_spillage/common.h"

class Inverse_Matrix
{
    public:
    ComplexMatrix A;
    void using_zpotrf( const ComplexMatrix &Sin);

    void using_zheev(const ComplexMatrix &in, ComplexMatrix &out);
    void init( const int &dim_in);

    Inverse_Matrix();
    ~Inverse_Matrix();

    private:
    int dim;
    double *e;
    int lwork;
    complex<double> *work2;
    double* rwork;
    int info;

    ComplexMatrix EA;
};

#endif
