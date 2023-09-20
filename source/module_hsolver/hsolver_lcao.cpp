#include "hsolver_lcao.h"

#include "diago_blas.h"
#include "module_base/timer.h"
#include "module_io/write_HS.h"

#ifdef __ELPA
#include "diago_elpa.h"
#endif

namespace hsolver
{

template <typename T>
void HSolverLCAO::solveTemplate(hamilt::Hamilt<std::complex<double>>* pHamilt,
                                psi::Psi<T>& psi,
                                elecstate::ElecState* pes,
                                const std::string method_in,
                                const bool skip_charge)
{
    ModuleBase::TITLE("HSolverLCAO", "solve");
    ModuleBase::timer::tick("HSolverLCAO", "solve");
    // select the method of diagonalization
    this->method = method_in;
    if (this->method == "scalapack_gvx")
    {
        if (pdiagh != nullptr)
        {
            if (pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = nullptr;
            }
        }
        if (pdiagh == nullptr)
        {
            pdiagh = new DiagoBlas();
            pdiagh->method = this->method;
        }
    }
#ifdef __ELPA
    else if (this->method == "genelpa")
    {
        if (pdiagh != nullptr)
        {
            if (pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = nullptr;
            }
        }
        if (pdiagh == nullptr)
        {
            pdiagh = new DiagoElpa();
            pdiagh->method = this->method;
        }
    }
#endif
    else if (this->method == "lapack")
    {
        ModuleBase::WARNING_QUIT("hsolver_lcao", "please fix lapack solver!!!");
        // We are not supporting diagonalization with lapack
        // until the obsolete globalc::hm is removed from
        // diago_lapack.cpp
        /*
        if (pdiagh != nullptr)
        {
            if (pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = nullptr;
            }
        }
        if (pdiagh == nullptr)
        {
            pdiagh = new DiagoLapack();
            pdiagh->method = this->method;
        }
        */
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This method of DiagH is not supported!");
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This method of DiagH is not supported!");
    }

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

        hamilt::MatrixBlock<T> h_mat, s_mat;

        psi.fix_k(ik);

        /// solve eigenvector and eigenvalue for H(k)
        double* p_eigenvalues = &(pes->ekb(ik, 0));
        this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);

        if (skip_charge)
        {
            pes->print_psi(psi);
        }
    }

    if (this->method != "genelpa" && this->method != "scalapack_gvx" && this->method != "lapack")
    {
        delete pdiagh;
        pdiagh = nullptr;
    }

    // used in nscf calculation
    if (skip_charge)
    {
        ModuleBase::timer::tick("HSolverLCAO", "solve");
        return;
    }

    // calculate charge by psi
    // called in scf calculation
    pes->psiToRho(psi);
    ModuleBase::timer::tick("HSolverLCAO", "solve");
}

int HSolverLCAO::out_mat_hs = 0;
int HSolverLCAO::out_mat_hsR = 0;
int HSolverLCAO::out_mat_t = 0;
int HSolverLCAO::out_mat_dh = 0;

void HSolverLCAO::solve(hamilt::Hamilt<std::complex<double>>* pHamilt,
                        psi::Psi<std::complex<double>>& psi,
                        elecstate::ElecState* pes,
                        const std::string method_in,
                        const bool skip_charge)
{
    this->solveTemplate(pHamilt, psi, pes, method, skip_charge);
}
void HSolverLCAO::solve(hamilt::Hamilt<std::complex<double>>* pHamilt,
                        psi::Psi<double>& psi,
                        elecstate::ElecState* pes,
                        const std::string method_in,
                        const bool skip_charge)
{
    this->solveTemplate(pHamilt, psi, pes, method, skip_charge);
}

void HSolverLCAO::hamiltSolvePsiK(hamilt::Hamilt<std::complex<double>>* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue)
{
    ModuleBase::TITLE("HSolverLCAO", "hamiltSolvePsiK");
    ModuleBase::timer::tick("HSolverLCAO", "hamiltSolvePsiK");
    pdiagh->diag(hm, psi, eigenvalue);
    ModuleBase::timer::tick("HSolverLCAO", "hamiltSolvePsiK");
}

void HSolverLCAO::hamiltSolvePsiK(hamilt::Hamilt<std::complex<double>>* hm, psi::Psi<double>& psi, double* eigenvalue)
{
    ModuleBase::TITLE("HSolverLCAO", "hamiltSolvePsiK");
    ModuleBase::timer::tick("HSolverLCAO", "hamiltSolvePsiK");
    pdiagh->diag(hm, psi, eigenvalue);
    ModuleBase::timer::tick("HSolverLCAO", "hamiltSolvePsiK");
}

} // namespace hsolver