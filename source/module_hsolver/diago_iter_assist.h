#ifndef DIAGOITERASSIST_H
#define DIAGOITERASSIST_H

#include "module_base/complexmatrix.h"
#include "module_psi/psi.h"
#include "module_hamilt/hamilt.h"

namespace hsolver
{

class DiagoIterAssist
{
  public:
    static double PW_DIAG_THR;
    static int PW_DIAG_NMAX;

    /// average steps of last cg diagonalization for each band.
    static double avg_iter;
    static bool need_subspace;

    // for CG diagonalization only
    static void diagH_subspace(hamilt::Hamilt* pHamilt,
                               const psi::Psi<std::complex<double>> &psi,
                               psi::Psi<std::complex<double>> &evc,
                               double *en,
                               int n_band = 0);
    // for initializing wave function , this is a template function
    static void diagH_subspace_init(hamilt::Hamilt* pHamilt,
                               const ModuleBase::ComplexMatrix &psi,
                               psi::Psi<std::complex<double>> &evc,
                               double *en);

    static void diagH_LAPACK(const int nstart,
                             const int nbands,
                             const ModuleBase::ComplexMatrix &hc,
                             const ModuleBase::ComplexMatrix &sc,
                             const int ldh, // nstart
                             double *e,
                             ModuleBase::ComplexMatrix &hvec);

    static bool test_exit_cond(const int &ntry, const int &notconv);
};

} // namespace hsolver

#endif