#include <mpi.h>
#include <complex>
#include "source_io/module_parameter/parameter.h"
#include <memory>
#ifdef __PEXSI
#include "diago_pexsi.h"
#include "source_base/tool_title.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_global.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "module_pexsi/simple_pexsi.h"
#include "module_pexsi/pexsi_solver.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <utility>

extern MPI_Comm DIAG_WORLD;

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
    if (static_cast<int>(mu_buffer.size()) < nspin)
    {
        mu_buffer.resize(nspin, pexsi::PEXSI_Solver::pexsi_mu);
    }
    this->resize_density_buffers(nspin);

}

template <typename T>
DiagoPexsi<T>::~DiagoPexsi()
{
}

template <typename T>
void DiagoPexsi<T>::resize_density_buffers(const int count)
{
    const int local_size = this->ParaV->nrow * this->ParaV->ncol;
    this->DM.resize(count);
    this->EDM.resize(count);
    this->dm_buffer_.resize(count);
    this->edm_buffer_.resize(count);
    for (int i = 0; i < count; ++i)
    {
        if (static_cast<int>(this->dm_buffer_[i].size()) != local_size)
        {
            this->dm_buffer_[i].assign(local_size, T{});
        }
        if (static_cast<int>(this->edm_buffer_[i].size()) != local_size)
        {
            this->edm_buffer_[i].assign(local_size, T{});
        }
        this->DM[i] = this->dm_buffer_[i].data();
        this->EDM[i] = this->edm_buffer_[i].data();
    }
    if (static_cast<int>(mu_buffer.size()) < count)
    {
        mu_buffer.resize(count, pexsi::PEXSI_Solver::pexsi_mu);
    }
}

template <typename T>
void DiagoPexsi<T>::begin_mu_search()
{
    this->has_mu_lower_ = false;
    this->has_mu_upper_ = false;
    this->mu_lower_ = 0.0;
    this->mu_upper_ = 0.0;
}

template <typename T>
void DiagoPexsi<T>::begin_k_loop()
{
    this->num_electron_sum_ = 0.0;
    this->num_electron_derivative_sum_ = 0.0;
    this->totalEnergyH = 0.0;
    this->totalEnergyS = 0.0;
    this->totalFreeEnergy = 0.0;
}

template <typename T>
void DiagoPexsi<T>::set_k_weight(const int ik, const double weight)
{
    if (ik >= static_cast<int>(this->k_weights_.size()))
    {
        this->k_weights_.resize(ik + 1, 1.0);
    }
    this->k_weights_[ik] = weight;
}

template <typename T>
bool DiagoPexsi<T>::finish_k_loop(const double target_nelec)
{
    return true;
}

template <>
bool DiagoPexsi<std::complex<double>>::finish_k_loop(const double target_nelec)
{
    const double residual = target_nelec - this->num_electron_sum_;
    if (std::abs(residual) <= pexsi::PEXSI_Solver::pexsi_elec_thr)
    {
        return true;
    }
    if (residual > 0.0)
    {
        this->has_mu_lower_ = true;
        this->mu_lower_ = mu_buffer[0];
    }
    else
    {
        this->has_mu_upper_ = true;
        this->mu_upper_ = mu_buffer[0];
    }

    double delta_mu = 0.0;
    if (std::abs(this->num_electron_derivative_sum_) <= DBL_MIN)
    {
        if (this->has_mu_lower_ && this->has_mu_upper_)
        {
            mu_buffer[0] = 0.5 * (this->mu_lower_ + this->mu_upper_);
            return false;
        }
        const double fallback_step = std::max(pexsi::PEXSI_Solver::pexsi_mu_guard, 0.5);
        delta_mu = residual > 0.0 ? fallback_step : -fallback_step;
    }
    else
    {
        delta_mu = residual / this->num_electron_derivative_sum_;
    }
    if (pexsi::PEXSI_Solver::pexsi_mu_guard > 0.0 && std::abs(this->num_electron_derivative_sum_) > DBL_MIN)
    {
        delta_mu = std::max(-pexsi::PEXSI_Solver::pexsi_mu_guard,
                            std::min(pexsi::PEXSI_Solver::pexsi_mu_guard, delta_mu));
    }
    if (mu_buffer.empty())
    {
        mu_buffer.push_back(pexsi::PEXSI_Solver::pexsi_mu);
    }
    mu_buffer[0] += delta_mu;
    return false;
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
    matcd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    const int ik = psi.get_current_k();
    if (ik < 0)
    {
        ModuleBase::WARNING_QUIT("DiagoPEXSI", "invalid k-point index for complex PEXSI path");
    }
    if (ik >= static_cast<int>(this->DM.size()))
    {
        this->resize_density_buffers(ik + 1);
    }
    if (mu_buffer.empty())
    {
        mu_buffer.push_back(pexsi::PEXSI_Solver::pexsi_mu);
    }

    MPI_Group grid_group;
    MPI_Group world_group;
    int grid_np = 0;
    MPI_Comm_size(DIAG_WORLD, &grid_np);
    MPI_Comm_group(DIAG_WORLD, &world_group);
    int grid_proc_range[3] = {0, (GlobalV::NPROC / grid_np) * grid_np - 1, GlobalV::NPROC / grid_np};
    MPI_Group_range_incl(world_group, 1, &grid_proc_range, &grid_group);

    double num_electron = 0.0;
    double num_electron_derivative = 0.0;
    double total_energy_h = 0.0;
    double total_energy_s = 0.0;
    double total_free_energy = 0.0;
    pexsi::simplePEXSIComplex(DIAG_WORLD,
                              DIAG_WORLD,
                              grid_group,
                              this->ParaV->blacs_ctxt,
                              PARAM.globalv.nlocal,
                              this->ParaV->nb,
                              this->ParaV->nrow,
                              this->ParaV->ncol,
                              'c',
                              h_mat.p,
                              s_mat.p,
                              PARAM.inp.nelec,
                              "PEXSIOPTION",
                              this->DM[ik],
                              this->EDM[ik],
                              total_energy_h,
                              total_energy_s,
                              total_free_energy,
                              mu_buffer[0],
                              mu_buffer[0],
                              &num_electron,
                              &num_electron_derivative);

    const double k_weight = ik < static_cast<int>(this->k_weights_.size()) ? this->k_weights_[ik] : 1.0;
    this->num_electron_sum_ += k_weight * num_electron;
    this->num_electron_derivative_sum_ += k_weight * num_electron_derivative;
    this->totalEnergyH += k_weight * total_energy_h;
    this->totalEnergyS += k_weight * total_energy_s;
    this->totalFreeEnergy += k_weight * total_free_energy;

    MPI_Group_free(&grid_group);
    MPI_Group_free(&world_group);
}

template class DiagoPexsi<double>;
template class DiagoPexsi<std::complex<double> >;

} // namespace hsolver
#endif
