#ifndef STOCHASTIC_WF_H
#define STOCHASTIC_WF_H

#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_psi/psi.h"

//----------------------------------------------
// Generate stochastic wave functions
//----------------------------------------------
template <typename T, typename Device = base_device::DEVICE_CPU>
class Stochastic_WF
{
  public:
    Stochastic_WF();

    ~Stochastic_WF();

    void init(K_Vectors* p_kv, const int npwx_in);

    // origin stochastic wavefunctions in CPU
    psi::Psi<T, base_device::DEVICE_CPU>* chi0_cpu = nullptr;
    // origin stochastic wavefunctions in GPU or CPU
    psi::Psi<T, Device>* chi0 = nullptr;
    // stochastic wavefunctions after in reciprocal space orthogonalized with KS wavefunctions
    psi::Psi<T, Device>* chiortho = nullptr;
    // sqrt(f(H))|chi>
    psi::Psi<T, Device>* shchi = nullptr;
    int nchi = 0;         ///< Total number of stochatic obitals
    int* nchip = nullptr; ///< The number of stochatic orbitals in current process of each k point.
    int nchip_max = 0;    ///< Max number of stochastic orbitals among all k points.
    int nks = 0;          ///< number of k-points
    int npwx = 0;         ///< max ngk[ik] in all processors
    int nbands_diag = 0;  ///< number of bands obtained from diagonalization
    int nbands_total = 0; ///< number of bands in total, nbands_total=nchi+nbands_diag;
    std::vector<int> ngk;   ///< ngk in klist
  public:
    // Tn(H)|chi>
    psi::Psi<T, Device>* chiallorder = nullptr;
    // allocate chiallorder
    void allocate_chiallorder(const int& norder);
    // chiallorder cost too much memories and should be cleaned after scf.
    void clean_chiallorder();

  public:
    // init stochastic orbitals
    void init_sto_orbitals(const int seed_in);
    // init stochastic orbitals from a large Ecut
    // It can test the convergence of SDFT with respect to Ecut
    void init_sto_orbitals_Ecut(const int seed_in,
                                const K_Vectors& kv,
                                const ModulePW::PW_Basis_K& wfcpw,
                                const int max_ecut);
    // allocate chi0
    void allocate_chi0();
    // update stochastic orbitals
    void update_sto_orbitals(const int seed_in);
    // init complete orbitals
    void init_com_orbitals();
    // sync chi0 from CPU to GPU
    void sync_chi0();

  protected:
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
};
#endif
