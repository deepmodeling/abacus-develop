//==========================================================
// AUTHOR : wangjp
// Data :2009-04
// Last Update:
//
// 09-05-10 modify SchmitOrth() diag_zhegvx() as static
// member function
//==========================================================

#ifndef DIAGODAVID_H
#define DIAGODAVID_H

#include "diagh.h"
#include "module_base/complexmatrix.h"

#if ((defined __CUDA) || (defined __ROCM))

#ifdef __CUDA
#include "src_pw/hamilt_pw.cuh"
#else
#include "src_pw/hamilt_pw_hip.h"
#endif

#else
#include "src_pw/hamilt_pw.h"
#endif
#include "src_pw/pw_basis.h"

namespace ModuleHSolver
{

class DiagoDavid : public DiagH
{
public:

    DiagoDavid(Hamilt_PW* hpw_in, const double *precondition_in);
    ~DiagoDavid();

    //this is the override function diag() for CG method
    void diag(
        ModuleHamilt::Hamilt* phm_in,
        ModulePsi::Psi<std::complex<double>> &phi,
        double *eigenvalue_in) override;

    static int PW_DIAG_NDIM;

private:

    int test_david;

    /// record for how many bands not have convergence eigenvalues
    int notconv;

    void cal_grad(
        const int& npw,
        const int& nbase,
        const int& notconv,
        ModuleBase::ComplexMatrix &basis,
        ModuleBase::ComplexMatrix &hp,
        ModuleBase::ComplexMatrix &sp,
        const ModuleBase::ComplexMatrix &vc,
        const int* unconv,
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
        const ModulePsi::Psi<std::complex<double>> &psi,
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
    
    void SchmitOrth(
        const int& npw,
        const int n_band,
        const int m,
        const ModuleBase::ComplexMatrix &psi,
        std::complex<double>* psi_m,
        std::complex<double>* spsi);

    void diag_zhegvx(
        const int& n,
        const int& m,
        const ModuleBase::ComplexMatrix &hc,
        const ModuleBase::ComplexMatrix &sc,
        const int& ldh,
        double* e,
        ModuleBase::ComplexMatrix &vc);

    void diag_mock(ModulePsi::Psi<std::complex<double>> &psi,
        double *eigenvalue_in);

    Hamilt_PW* hpw;
    const double* precondition;


};

}

#endif
