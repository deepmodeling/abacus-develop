#include "hsolver_lcao.h"

#include "diago_blas.h"
#include "diago_elpa.h"

namespace hsolver
{

template <typename T>
void HSolverLCAO::solveTemplate(hamilt::Hamilt* pHamilt, psi::Psi<T>& psi, elecstate::ElecState* pes)
{
    // select the method of diagonalization
    if (this->method == "genelpa")
        pdiagh = new DiagoElpa();
    else if (this->method == "scalapack_gvx")
        pdiagh = new DiagoBlas();
    else
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This method of DiagH is not supported!");

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

        psi.fix_k(ik);

        /// solve eigenvector and eigenvalue for H(k)
        double* p_eigenvalues = &(pes->ekb(ik, 0));
        this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);
    }
    pes->psiToRho(psi);
}

int HSolverLCAO::out_mat_hs = 0;
int HSolverLCAO::out_mat_hsR = 0;

void HSolverLCAO::solve(hamilt::Hamilt* pHamilt, psi::Psi<std::complex<double>>& psi, elecstate::ElecState* pes)
{
    this->solveTemplate(pHamilt, psi, pes);
}
void HSolverLCAO::solve(hamilt::Hamilt* pHamilt, psi::Psi<double>& psi, elecstate::ElecState* pes)
{
    this->solveTemplate(pHamilt, psi, pes);
}

void HSolverLCAO::hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue)
{
    pdiagh->diag(hm, psi, eigenvalue);
}

void HSolverLCAO::hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<double>& psi, double* eigenvalue)
{
    pdiagh->diag(hm, psi, eigenvalue);
}

} // namespace hsolver