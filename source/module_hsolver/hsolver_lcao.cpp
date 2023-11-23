#include "hsolver_lcao.h"

#include "diago_blas.h"
#include "module_base/timer.h"
#include "module_io/write_HS.h"
#ifdef __CUSOLVER_LCAO
#include "diago_cusolver.h"
#endif
#ifdef __ELPA
#include "diago_elpa.h"
#endif

namespace hsolver
{

template <typename T>
void HSolverLCAO<T>::solveTemplate(hamilt::Hamilt<T>* pHamilt,
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
        if (this->pdiagh != nullptr)
        {
            if (this->pdiagh->method != this->method)
            {
                delete[] this->pdiagh;
                this->pdiagh = nullptr;
            }
        }
        if (this->pdiagh == nullptr)
        {
            this->pdiagh = new DiagoBlas<T>();
            this->pdiagh->method = this->method;
        }
    }
#ifdef __ELPA
    else if (this->method == "genelpa")
    {
        if (this->pdiagh != nullptr)
        {
            if (this->pdiagh->method != this->method)
            {
                delete[] this->pdiagh;
                this->pdiagh = nullptr;
            }
        }
        if (this->pdiagh == nullptr)
        {
            this->pdiagh = new DiagoElpa<T>();
            this->pdiagh->method = this->method;
        }
    }
#endif
#ifdef __CUSOLVER_LCAO
    else if (this->method == "cusolver")
    {
        if (this->pdiagh != nullptr)
        {
            if (this->pdiagh->method != this->method)
            {
                delete[] this->pdiagh;
                this->pdiagh = nullptr;
            }
        }
        if (this->pdiagh == nullptr)
        {
            this->pdiagh = new DiagoCusolver<T>();
            this->pdiagh->method = this->method;
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
        if (this->pdiagh != nullptr)
        {
            if (this->pdiagh->method != this->method)
            {
                delete[] this->pdiagh;
                this->pdiagh = nullptr;
            }
        }
        if (this->pdiagh == nullptr)
        {
            this->pdiagh = new DiagoLapack();
            this->pdiagh->method = this->method;
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

    if (this->method != "genelpa" && this->method != "scalapack_gvx" && this->method != "lapack"
#ifdef __CUSOLVER_LCAO
        && this->method != "cusolver"
#endif
    )
    {
        delete this->pdiagh;
        this->pdiagh = nullptr;
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
template <typename T>
int HSolverLCAO<T>::out_mat_hs = 0;
template <typename T>
int HSolverLCAO<T>::out_mat_hsR = 0;
template <typename T>
int HSolverLCAO<T>::out_mat_t = 0;
template <typename T>
int HSolverLCAO<T>::out_mat_dh = 0;

template <typename T>
void HSolverLCAO<T>::solve(hamilt::Hamilt<T>* pHamilt,
    psi::Psi<T>& psi,
    elecstate::ElecState* pes,
    const std::string method_in,
    const bool skip_charge)
{
    this->solveTemplate(pHamilt, psi, pes, this->method, skip_charge);
}

template <typename T>
void HSolverLCAO<T>::hamiltSolvePsiK(hamilt::Hamilt<T>* hm, psi::Psi<T>& psi, double* eigenvalue)
{
    ModuleBase::TITLE("HSolverLCAO", "hamiltSolvePsiK");
    ModuleBase::timer::tick("HSolverLCAO", "hamiltSolvePsiK");
    this->pdiagh->diag(hm, psi, eigenvalue);
    ModuleBase::timer::tick("HSolverLCAO", "hamiltSolvePsiK");
}

template class HSolverLCAO<double>;
template class HSolverLCAO<std::complex<double>>;

} // namespace hsolver