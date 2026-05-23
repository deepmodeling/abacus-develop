#include <mpi.h>
#include <complex>
#include "source_io/module_parameter/parameter.h"
#include <memory>
#ifdef __PEXSI
#include "diago_pexsi.h"
#include "source_base/tool_title.h"
#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "module_pexsi/pexsi_solver.h"

#include <utility>

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{
template <typename T>
std::vector<double> DiagoPexsi<T>::mu_buffer;

template <typename T>
DiagoPexsi<T>::DiagoPexsi(const Parallel_Orbitals* ParaV_in, std::unique_ptr<pexsi::IPexsiSolver> solver_in)
{
    this->ParaV = ParaV_in;
    this->ps = std::move(solver_in);
    if (this->ps == nullptr)
    {
        this->ps = std::make_unique<pexsi::PEXSI_Solver>();
    }

    int nspin = PARAM.inp.nspin;
    if (PARAM.inp.nspin == 4)
    {
        nspin = 1;
    }
    mu_buffer.assign(nspin, pexsi::PEXSI_Solver::pexsi_mu);

    const int local_size = ParaV->nrow * ParaV->ncol;

    this->DM.resize(nspin);
    this->EDM.resize(nspin);
    this->dm_buffer_.resize(nspin);
    this->edm_buffer_.resize(nspin);
    for (int i = 0; i < nspin; i++)
    {
        this->dm_buffer_[i].assign(local_size, T{});
        this->edm_buffer_[i].assign(local_size, T{});
        this->DM[i] = this->dm_buffer_[i].data();
        this->EDM[i] = this->edm_buffer_[i].data();
    }

}

template <typename T>
DiagoPexsi<T>::~DiagoPexsi()
{
}

template <>
void DiagoPexsi<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);
    int ik = psi.get_current_k();
    if (ik < 0 || ik >= static_cast<int>(this->DM.size()))
    {
        ModuleBase::WARNING_QUIT("DiagoPEXSI",
                                 "PEXSI real path only has density buffers for Gamma/spin channels; multi-k requires "
                                 "the complex expert-routine path");
    }
    this->ps->prepare(this->ParaV->blacs_ctxt,
                      this->ParaV->nb,
                      this->ParaV->nrow,
                      this->ParaV->ncol,
                      h_mat.p,
                      s_mat.p,
                      DM[ik],
                      EDM[ik]);
    this->ps->solve(mu_buffer[ik]);
    this->totalFreeEnergy = this->ps->get_totalFreeEnergy();
    this->totalEnergyH = this->ps->get_totalEnergyH();
    this->totalEnergyS = this->ps->get_totalEnergyS();
    mu_buffer[ik] = this->ps->get_mu();
}

template <>
void DiagoPexsi<std::complex<double>>::diag(hamilt::Hamilt<std::complex<double>>* phm_in,
                                            psi::Psi<std::complex<double>>& psi,
                                            double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    ModuleBase::WARNING_QUIT("DiagoPEXSI", "PEXSI is not completed for multi-k case");
}

template class DiagoPexsi<double>;
template class DiagoPexsi<std::complex<double> >;

} // namespace hsolver
#endif
