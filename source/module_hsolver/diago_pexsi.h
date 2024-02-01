#ifndef DIGAOPEXSI_H
#define DIGAOPEXSI_H

#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_pexsi/pexsi_solver.h"

namespace hsolver
{

template <typename T>
class DiagoPexsi : public DiagH<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    DiagoPexsi(const Parallel_Orbitals* ParaV_in)
    {
        this->ParaV = ParaV_in;
    }
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in) override;
    const Parallel_Orbitals* ParaV;
    T* DM;
    double* EDM;
    double totalEnergyH;
    double totalEnergyS;
    double totalFreeEnergy;
    pexsi::PEXSI_Solver* ps;

    //==========================================================
    // PEXSI related variables
    //==========================================================
    static int pexsi_npole;
    static int pexsi_inertia;
    static int pexsi_nmax;
    // static int pexsi_symbolic;
    static int pexsi_comm;
    static int pexsi_storage;
    static int pexsi_ordering;
    static int pexsi_row_ordering;
    static int pexsi_nproc;
    static int pexsi_symm;
    static int pexsi_trans;
    static int pexsi_method;
    static int pexsi_nproc_pole;
    // static double pexsi_spin = 2;
    static double pexsi_temp;
    static double pexsi_gap;
    static double pexsi_delta_e;
    static double pexsi_mu_lower;
    static double pexsi_mu_upper;
    static double pexsi_mu;
    static double pexsi_mu_thr;
    static double pexsi_mu_expand;
    static double pexsi_mu_guard;
    static double pexsi_elec_thr;
    static double pexsi_zero_thr;

    static MPI_Group grid_group;
};
} // namespace hsolver

#endif
