//==========================================================
// AUTHOR : wangjp
// Data :2009-04
// Last Update:
//
// 09-05-10 modify SchmitOrth() diag_zhegvx() as static
// member function
//==========================================================

#ifndef Diago_David_H
#define Diago_David_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"

#if ((defined __CUDA) || (defined __ROCM))

#ifdef __CUDA
#include "hamilt_pw.cuh"
#else
#include "hamilt_pw_hip.h"
#endif

#else
#include "hamilt_pw.h"
#endif

class Diago_David
{
public:

    Diago_David(Hamilt_PW* phamilt);
    ~Diago_David();

    void SchmitOrth(
        const int& npw,
        const int n_band,
        const int m,
        const ModuleBase::ComplexMatrix &psi,
        std::complex<double>* psi_m,
        std::complex<double>* spsi);

    static void diag_zhegvx(
        const int& n,
        const int& m,
        const ModuleBase::ComplexMatrix &hc,
        const ModuleBase::ComplexMatrix &sc,
        const int& ldh,
        double* e,
        ModuleBase::ComplexMatrix &vc);

    void diag(
        ModuleBase::ComplexMatrix &psi,
        double *en,
        const int &npw,
        const int &nband,
        const double *precondition,
        const int order,
        const double &eps,
        const int &maxiter,
        int &notconv,
        double &avg_iter);

private:

    Hamilt_PW* hpw;

    int test_david;

    void cal_grad(
        const int& npw,
        const int& nbase,
        const int& notconv,
        ModuleBase::ComplexMatrix &basis,
        ModuleBase::ComplexMatrix &hp,
        ModuleBase::ComplexMatrix &sp,
        const ModuleBase::ComplexMatrix &vc,
        const int* unconv,
        const double* precondition,
        const double* en,
        std::complex<double>* hpsi,
        std::complex<double>* spsi,
        std::complex<double>* ppsi,
        std::complex<double>* respsi);

    void cal_elem(
        const int& npw,
        int& nbase,
        const int& notconv,
        const ModuleBase::ComplexMatrix &basis,
        const ModuleBase::ComplexMatrix &hp,
        const ModuleBase::ComplexMatrix &sp,
        ModuleBase::ComplexMatrix &hc,
        ModuleBase::ComplexMatrix &sc);


    void refresh(
        const int& npw,
        const int& nband,
        int& nbase,
        const double* en,
        const ModuleBase::ComplexMatrix &psi,
        ModuleBase::ComplexMatrix &basis,
        ModuleBase::ComplexMatrix &hp,
        ModuleBase::ComplexMatrix &sp,
        ModuleBase::ComplexMatrix &hc,
        ModuleBase::ComplexMatrix &sc,
        ModuleBase::ComplexMatrix &vc);

    void cal_err(
        const int& npw,
        const int& nband,
        const int& nbase,
        const ModuleBase::ComplexMatrix &vc,
        const ModuleBase::ComplexMatrix &hp,
        const ModuleBase::ComplexMatrix &basis,
        const double* en,
        std::complex<double>* respsi);

};

#endif
