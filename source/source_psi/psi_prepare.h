#ifndef PSI_PREPARE_H
#define PSI_PREPARE_H
#include "source_hamilt/hamilt.h"
#include "source_psi/psi_base.h"
#include "source_psi/psi_prepare_base.h"

namespace psi
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class PSIPrepare : public PSIPrepareBase
{
  public:
    PSIPrepare(const std::string& init_wfc_in,
            const std::string& ks_solver_in,
            const std::string& basis_type_in,
            const int& rank,
            const UnitCell& ucell,
            const Structure_Factor& sf,
            const std::vector<int>& ik2iktot,
            const int& nkstot,
            const int& lmaxkb,
            const ModulePW::PW_Basis_K& pw_wfc);

    ~PSIPrepare(){};

    ///@brief prepare the wavefunction initialization
    ///@param random_seed seed for random initialization
    ///@param istep current ion/relax step; informational warnings are only
    ///       printed on the first step to avoid spamming relax output
    void prepare_init(const int& random_seed, const int istep);

    /**
     * @brief initialize the wavefunction
     *
     * @param psi store the wavefunction
     * @param kspw_psi Kohn-Sham wavefunction in plane-wave basis
     * @param p_hamilt Hamiltonian operator
     * @param ofs_running output stream for running information
     */
    void initialize_psi(Psi<std::complex<double>>* psi,
                        psi::Psi<T, Device>* kspw_psi,
                        hamilt::Hamilt<T, Device>* p_hamilt,
                        std::ofstream& ofs_running);

    /**
     * @brief initialize NAOs in plane wave basis, only for LCAO_IN_PW
     *
     */
    void initialize_lcao_in_pw(Psi<T>* psi_local, std::ofstream& ofs_running);

    std::unique_ptr<psi_base<T>> psi_initer;

  private:

    // wavefunction initialization type
    std::string init_wfc = "none";

    // Kohn-Sham solver type
    std::string ks_solver = "none";

    // basis type
    std::string basis_type = "none";

    // pw basis
    const ModulePW::PW_Basis_K& pw_wfc;

    // local->global k-point mapping
    const std::vector<int>& ik2iktot_;
    // total number of k-points
    const int nkstot_;

    // unit cell
    const UnitCell& ucell;

    // structure factor
    const Structure_Factor& sf;

    // max angular momentum for non-local projectors
    const int lmaxkb;

    Device* ctx = {};                      ///< device
    base_device::DEVICE_CPU* cpu_ctx = {}; ///< CPU device
    const int rank;                        ///< MPI rank

    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
};

void allocate_psi(Psi<std::complex<double>>*& psi, const int& nks, const std::vector<int>& ngk, const int& nbands, const int& npwx);

} // namespace psi
#endif
