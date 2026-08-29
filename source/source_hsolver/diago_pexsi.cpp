#include <mpi.h>
#include <complex>
#include <memory>
#ifdef __PEXSI
#include "diago_pexsi.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "module_pexsi/pexsi_solver.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{
template <typename T>
std::vector<double> DiagoPexsi<T>::mu_buffer;

template <typename T>
DiagoPexsi<T>::DiagoPexsi(const Parallel_Orbitals* ParaV_in,
                          const int nspin_in,
                          const int nlocal_in,
                          const double nelec_in)
{
    this->nspin_dm = (nspin_in == 4) ? 1 : nspin_in;
    this->nlocal = nlocal_in;
    this->nelec = nelec_in;

    mu_buffer.resize(this->nspin_dm);
    for (int i = 0; i < this->nspin_dm; i++)
    {
        mu_buffer[i] = pexsi::PEXSI_Solver::pexsi_mu;
    }

    this->ParaV = ParaV_in;
    this->ps.reset(new pexsi::PEXSI_Solver());

    this->DM.resize(this->nspin_dm);
    this->EDM.resize(this->nspin_dm);
    for (int i = 0; i < this->nspin_dm; i++)
    {
        this->DM[i] = new T[ParaV->nrow * ParaV->ncol];
        this->EDM[i] = new T[ParaV->nrow * ParaV->ncol];
    }

}

template <typename T>
DiagoPexsi<T>::~DiagoPexsi()
{
    for (int i = 0; i < this->nspin_dm; i++)
    {
        delete[] this->DM[i];
        delete[] this->EDM[i];
    }

}

template <>
void DiagoPexsi<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoPEXSI", "diag");
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    int ik = psi.get_current_k();
    this->ps->prepare(this->ParaV->blacs_ctxt,
                      this->ParaV->nb,
                      this->ParaV->nrow,
                      this->ParaV->ncol,
                      this->nlocal,
                      this->nelec,
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