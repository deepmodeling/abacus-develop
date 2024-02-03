#include <mpi.h>
#include <complex>
#ifdef __PEXSI
#include "c_pexsi_interface.h"
#include "diago_pexsi.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_pexsi/pexsi_solver.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{
template <>
int DiagoPexsi<double>::pexsi_npole = 0;
template <>
bool DiagoPexsi<double>::pexsi_inertia = 0;
template <>
int DiagoPexsi<double>::pexsi_nmax = 0;
// template <>
// int DiagoPexsi<double>::pexsi_symbolic = 0;
template <>
bool DiagoPexsi<double>::pexsi_comm = 0;
template <>
bool DiagoPexsi<double>::pexsi_storage = 0;
template <>
int DiagoPexsi<double>::pexsi_ordering = 0;
template <>
int DiagoPexsi<double>::pexsi_row_ordering = 0;
template <>
int DiagoPexsi<double>::pexsi_nproc = 0;
template <>
bool DiagoPexsi<double>::pexsi_symm = 0;
template <>
bool DiagoPexsi<double>::pexsi_trans = 0;
template <>
int DiagoPexsi<double>::pexsi_method = 0;
template <>
int DiagoPexsi<double>::pexsi_nproc_pole = 0;
// template <>
// double DiagoPexsi<double>::pexsi_spin = 2;
template <>
double DiagoPexsi<double>::pexsi_temp = 0.0;
template <>
double DiagoPexsi<double>::pexsi_gap = 0.0;
template <>
double DiagoPexsi<double>::pexsi_delta_e = 0.0;
template <>
double DiagoPexsi<double>::pexsi_mu_lower = 0.0;
template <>
double DiagoPexsi<double>::pexsi_mu_upper = 0.0;
template <>
double DiagoPexsi<double>::pexsi_mu = 0.0;
template <>
double DiagoPexsi<double>::pexsi_mu_thr = 0.0;
template <>
double DiagoPexsi<double>::pexsi_mu_expand = 0.0;
template <>
double DiagoPexsi<double>::pexsi_mu_guard = 0.0;
template <>
double DiagoPexsi<double>::pexsi_elec_thr = 0.0;
template <>
double DiagoPexsi<double>::pexsi_zero_thr = 0.0;

template <>
int DiagoPexsi<std::complex<double>>::pexsi_npole = 0;
template <>
bool DiagoPexsi<std::complex<double>>::pexsi_inertia = 0;
template <>
int DiagoPexsi<std::complex<double>>::pexsi_nmax = 0;
// template <>
// int DiagoPexsi<std::complex<double>>::pexsi_symbolic = 0;
template <>
bool DiagoPexsi<std::complex<double>>::pexsi_comm = 0;
template <>
bool DiagoPexsi<std::complex<double>>::pexsi_storage = 0;
template <>
int DiagoPexsi<std::complex<double>>::pexsi_ordering = 0;
template <>
int DiagoPexsi<std::complex<double>>::pexsi_row_ordering = 0;
template <>
int DiagoPexsi<std::complex<double>>::pexsi_nproc = 0;
template <>
bool DiagoPexsi<std::complex<double>>::pexsi_symm = 0;
template <>
bool DiagoPexsi<std::complex<double>>::pexsi_trans = 0;
template <>
int DiagoPexsi<std::complex<double>>::pexsi_method = 0;
template <>
int DiagoPexsi<std::complex<double>>::pexsi_nproc_pole = 0;
// template <>
// double DiagoPexsi<std::complex<double>>::pexsi_spin = 2;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_temp = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_gap = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_delta_e = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_mu_lower = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_mu_upper = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_mu = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_mu_thr = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_mu_expand = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_mu_guard = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_elec_thr = 0.0;
template <>
double DiagoPexsi<std::complex<double>>::pexsi_zero_thr = 0.0;

template <>
void DiagoPexsi<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);
    MPI_Comm COMM_DIAG = MPI_COMM_WORLD;
    int ik = psi.get_current_k();
    this->ps = new pexsi::PEXSI_Solver(this->ParaV->blacs_ctxt,
                                       this->ParaV->nb,
                                       this->ParaV->nrow,
                                       this->ParaV->ncol,
                                       h_mat.p,
                                       s_mat.p,
                                       this->totalEnergyH,
                                       this->totalEnergyS,
                                       this->totalFreeEnergy);
    this->ps->solve(mu_buffer[ik]);
    this->EDM.push_back(this->ps->get_EDM());
    this->DM.push_back(this->ps->get_DM());
    this->totalFreeEnergy = this->ps->get_totalFreeEnergy();
    this->totalEnergyH = this->ps->get_totalEnergyH();
    this->totalEnergyS = this->ps->get_totalEnergyS();
    this->mu_buffer[ik] = this->ps->get_mu();
}

template <>
void DiagoPexsi<std::complex<double>>::diag(hamilt::Hamilt<std::complex<double>>* phm_in,
                                            psi::Psi<std::complex<double>>& psi,
                                            double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    ModuleBase::WARNING_QUIT("DiagoPEXSI", "PEXSI is not completed for multi-k case");
}

} // namespace hsolver
#endif