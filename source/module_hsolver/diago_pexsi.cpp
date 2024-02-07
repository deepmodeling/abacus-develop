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