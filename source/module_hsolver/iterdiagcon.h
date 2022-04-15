#ifndef ITERDIAGCONTROL_H
#define ITERDIAGCONTROL_H

#include "module_psi/psi.h"
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

namespace ModuleHSolver
{

class IterDiagControl
{
    public: 
    static double PW_DIAG_THR;
    static int PW_DIAG_NMAX;
    
    /// average steps of last cg diagonalization for each band.
    static double avg_iter;

    static void diagH_subspace(
        Hamilt_PW* phm,
        const ModulePsi::Psi<std::complex<double>> &psi,
        ModulePsi::Psi<std::complex<double>> &evc,
        double *en,
        int n_band=0);

    static void diagH_LAPACK(
        const int nstart,
        const int nbands,
        const ModuleBase::ComplexMatrix &hc,
        const ModuleBase::ComplexMatrix &sc,
        const int ldh, // nstart
        double *e,
        ModuleBase::ComplexMatrix &hvec);

    static bool test_exit_cond(const int &ntry, const int &notconv);
};

}

#endif