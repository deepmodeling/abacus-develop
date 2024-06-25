#ifndef LOCAL_ORBITAL_WFC
#define LOCAL_ORBITAL_WFC

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_control.h" // mohan add 2021-05-24
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_psi/psi.h"

class Local_Orbital_wfc
{
  public:
    Local_Orbital_wfc();
    ~Local_Orbital_wfc();
    ///=========================================
    /// grid wfc
    /// used to generate density matrix: LOC.DM_R,
    /// which is used to calculate the charge density.
    /// which is got after the diagonalization of
    /// std::complex Hamiltonian matrix.
    ///=========================================
    std::complex<double>*** wfc_k_grid; // [NK, GlobalV::NBANDS, GlobalV::NLOCAL]
    std::complex<double>* wfc_k_grid2;  // [NK*GlobalV::NBANDS*GlobalV::NLOCAL]

    // pointer to const Parallel_Orbitals object
    const Parallel_Orbitals* ParaV;
    // pointer to const Grid_Technique object, although have no idea about what it is...
    // the name of Grid_Technique should be changed to be more informative
    const Grid_Technique* gridt;

    /// read wavefunction coefficients: WFC_NAO_K/GAMMA*.txt
    void gamma_file(psi::Psi<double>* psid, elecstate::ElecState* pelec);
    void allocate_k(const int& lgd,
                    psi::Psi<std::complex<double>>* psi,
                    elecstate::ElecState* pelec,
                    const int& nks,
                    const int& nkstot,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                    const int& istep);

    //=========================================
    // Init Cij, make it satisfy 2 conditions:
    // (1) Unit
    // (2) Orthogonal <i|S|j>= \delta{ij}
    //=========================================
    // void init_Cij(const bool change_c = 1);

    ///=========================================
    /// Parallel: map of index in 2D distribution: global<->local
    ///=========================================
    static int globalIndex(int localindex, int nblk, int nprocs, int myproc);
    static int localIndex(int globalindex, int nblk, int nprocs, int& myproc);

    int error = 0;

  private:
    bool wfck_flag;
    bool complex_flag;
    bool allocate_flag;
    int nks;
};
#endif
