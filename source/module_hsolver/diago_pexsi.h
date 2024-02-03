#ifndef DIGAOPEXSI_H
#define DIGAOPEXSI_H

#include <vector>
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
    std::vector<T*> DM;
    std::vector<T*> EDM;
    double totalEnergyH;
    double totalEnergyS;
    double totalFreeEnergy;
    pexsi::PEXSI_Solver* ps;

    //==========================================================
    // PEXSI related variables
    //==========================================================
    /** 
     * @brief  Number of terms in the pole expansion.
     */ 
    static int pexsi_npole;
    /** 
     * @brief  Whether inertia counting is used at the very beginning.
     */ 
    static bool pexsi_inertia;
    /** 
     * @brief  Maximum number of PEXSI iterations after each inertia counting procedure.
     */ 
    static int pexsi_nmax;
    /** 
     * @brief  Whether to construct PSelInv communication pattern.
     */ 
    static bool pexsi_comm;
    /** 
     * @brief  Whether to use symmetric storage space used by the Selected Inversion algorithm for symmetric matrices.  
     */ 
    static bool pexsi_storage;
    /** 
     * @brief  Ordering strategy for factorization and selected inversion. 
     */ 
    static int pexsi_ordering;
    /** 
     * @brief  row permutation strategy for factorization and selected inversion.  
     */ 
    static int pexsi_row_ordering;
    /** 
     * @brief  Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering == 0.
     */ 
    static int pexsi_nproc;
    /** 
     * @brief  Matrix structure.
     * - = 0   : Unsymmetric matrix
     * - = 1   : Symmetric matrix (default).
     */ 
    static bool pexsi_symm;
    /** 
     * @brief  Transpose.
     * - = 0   : Factor non transposed matrix (default).
     * - = 1   : Factor transposed matrix.
     */ 
    static bool pexsi_trans;
    /** 
     * @brief  The pole expansion method to be used.
     * - = 1   : Cauchy Contour Integral method used.
     * - = 2   : Moussa optimized method.
     */ 
    static int pexsi_method;
    /** 
     * @brief  The point parallelizaion of PEXSI.
     * - = 2  : Recommend two points parallelization
     */ 
    static int pexsi_nproc_pole;
    /** 
     * @brief  Temperature, in the same unit as H 
     */ 
    static double pexsi_temp;
    /** 
     * @brief  Spectral gap. **Note** This can be set to be 0 in most cases.
     */ 
    static double pexsi_gap;
    /** 
     * @brief  An upper bound for the spectral radius of \f$S^{-1} H\f$.
     */ 
    static double pexsi_delta_e;
    /** 
     * @brief  Initial guess of lower bound for mu.
     */ 
    static double pexsi_mu_lower;
    /** 
     * @brief  Initial guess of upper bound for mu.
     */ 
    static double pexsi_mu_upper;
    /** 
     * @brief  Initial guess for mu (for the solver) (AG)
     */ 
    static double pexsi_mu;
    /** 
     * @brief  Stopping criterion in terms of the chemical potential for the inertia counting procedure.
     */ 
    static double pexsi_mu_thr;
    /** 
     * @brief  If the chemical potential is not in the initial interval, the interval is expanded by muInertiaExpansion.
     */ 
    static double pexsi_mu_expand;
    /** 
     * @brief  Safe guard criterion in terms of the chemical potential to reinvoke the inertia counting procedure.
     */ 
    static double pexsi_mu_guard;
    /** 
     * @brief  Stopping criterion of the %PEXSI iteration in terms of the number of electrons compared to numElectronExact.
     */ 
    static double pexsi_elec_thr;
    /** 
     * @brief  Stopping criterion for the zero threshold.
     */ 
    static double pexsi_zero_thr;
};
} // namespace hsolver

#endif
