#include "hsolver_lcao.h"

#include "diago_blas.h"
#include "diago_elpa.h"
#include "diago_lapack.h"
#include "module_base/timer.h"
#include "src_io/write_HS.h"
#include "src_pw/global.h"

namespace hsolver
{

template <typename T>
void HSolverLCAO::solveTemplate(hamilt::Hamilt* pHamilt,
                                psi::Psi<T>& psi,
                                elecstate::ElecState* pes,
                                const std::string method_in,
                                const bool skip_charge)
{
    ModuleBase::TITLE("HSolverLCAO", "solve");
    ModuleBase::timer::tick("HSolverLCAO", "solve");
    // select the method of diagonalization
    this->method = method_in;
    if (this->method == "genelpa")
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
    else if (this->method == "scalapack_gvx")
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
    else if (this->method == "lapack")
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
            pdiagh = new DiagoLapack();
            pdiagh->method = this->method;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This method of DiagH is not supported!");
    }

    pHamilt->constructHamilt();

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {  
        /// update H(k) for each k point
        pHamilt->updateHk(ik);
        hamilt::MatrixBlock<T> h_mat, s_mat;
        pHamilt->matrix(h_mat, s_mat);

        if(this->out_mat_hs_k)
        {
            bool write_binary = false; // LiuXh, 2017-03-21
            // if set write_binary = true, there would be error in soc-multi-core calculation, noted by zhengdy-soc
            int ind_k;//which k point
            if(GlobalV::NSPIN==1 || GlobalV::NSPIN==4)
            {
                ind_k = ik;
            }
            else
            {
                ind_k = ik%(GlobalC::kv.nkstot/2);
            }
            HS_Matrix::save_HS(
                ik,
                h_mat.p,
                s_mat.p,
                write_binary,
                "data_k" + std::to_string(ind_k+1) + "_s" + std::to_string(GlobalC::kv.isk[ik]+1),
                this->ParaV[0]); // LiuXh, 2017-03-21
        }

        psi.fix_k(ik);

        /// solve eigenvector and eigenvalue for H(k)
        double* p_eigenvalues = &(pes->ekb(ik, 0));
        this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);

        if(skip_charge) 
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

int HSolverLCAO::out_mat_hs_k = 0;
int HSolverLCAO::out_mat_hs_r = 0;

void HSolverLCAO::solve(hamilt::Hamilt* pHamilt,
                        psi::Psi<std::complex<double>>& psi,
                        elecstate::ElecState* pes,
                        const std::string method_in,
                        const bool skip_charge)
{
    this->solveTemplate(pHamilt, psi, pes, method, skip_charge);
}
void HSolverLCAO::solve(hamilt::Hamilt* pHamilt,
                        psi::Psi<double>& psi,
                        elecstate::ElecState* pes,
                        const std::string method_in,
                        const bool skip_charge)
{
    this->solveTemplate(pHamilt, psi, pes, method, skip_charge);
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