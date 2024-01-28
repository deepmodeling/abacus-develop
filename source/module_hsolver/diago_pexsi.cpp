#include <complex>
#ifdef __PEXSI
#include "diago_pexsi.h"

#include "c_pexsi_interface.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "pexsi/pexsi_solver.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{

template<>
void DiagoPexsi<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);
    MPI_Comm COMM_DIAG = MPI_COMM_WORLD;
    this->ps = new pexsi::PEXSI_Solver(this->ParaV->blacs_ctxt,
                                this->ParaV->nb,
                                this->ParaV->nrow,
                                this->ParaV->ncol,
                                h_mat.p,
                                s_mat.p,
                                this->DM,
                                this->EDM,
                                this->totalEnergyH,
                                this->totalEnergyS,
                                this->totalFreeEnergy);
    this->ps->solve();
    this->EDM = this->ps->EDM;
    this->DM = this->ps->DM; // loc.dm_gamma[ik] loc.dm_gamma[0]?
    this->totalFreeEnergy = this->ps->totalFreeEnergy;
    this->totalEnergyH = this->ps->totalEnergyH;
    this->totalEnergyS = this->ps->totalEnergyS;
}

template<>
void DiagoPexsi<std::complex<double>>::diag(hamilt::Hamilt<std::complex<double>>* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    ModuleBase::WARNING_QUIT("DiagoPEXSI", "PEXSI is not completed for multi-k case");
    
}

} // namespace hsolver
#endif