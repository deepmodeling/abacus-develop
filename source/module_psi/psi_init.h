#ifndef PSI_INIT_H
#define PSI_INIT_H
#include "module_hamilt_general/hamilt.h"
#include "module_psi/psi_initializer.h"

namespace psi
{

// This class is used to initialize the wavefunction
template <typename T, typename Device = base_device::DEVICE_CPU>
class PSIInit
{
  public:
    PSIInit(const std::string& init_wfc_in,
            const std::string& ks_solver_in,
            const std::string& basis_type_in,
            const int& rank,
            const UnitCell& ucell,
            const Structure_Factor& sf,
            const K_Vectors& kv_in,
            const pseudopot_cell_vnl& nlpp,
            const ModulePW::PW_Basis_K& pw_wfc);
    ~PSIInit(){};

    ///@brief prepare the wavefunction initialization
    void prepare_init(const int& random_seed);

    //------------------------ only for psi_initializer --------------------
    /**
     * @brief initialize the wavefunction
     *
     * @param psi store the wavefunction
     * @param p_hamilt Hamiltonian operator
     * @param ofs_running output stream for running information
     * @param is_already_initpsi whether psi has been initialized
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

    // psi_initializer<T, Device>* psi_initer = nullptr;
    // change to use smart pointer to manage the memory, and avoid memory leak
    // while the std::make_unique() is not supported till C++14,
    // so use the new and std::unique_ptr to manage the memory, but this makes new-delete not symmetric
    std::unique_ptr<psi_initializer<T>> psi_initer;

  private:
    // wavefunction initialization type
    std::string init_wfc = "none";

    // Kohn-Sham solver type
    std::string ks_solver = "none";

    // basis type
    std::string basis_type = "none";

    // pw basis
    const ModulePW::PW_Basis_K& pw_wfc;

    // parallel kpoints
    const K_Vectors& kv;

    // unit cell
    const UnitCell& ucell;

    // structure factor
    const Structure_Factor& sf;

    // nonlocal pseudopotential
    const pseudopot_cell_vnl& nlpp;

    Device* ctx = {};                      ///< device
    base_device::DEVICE_CPU* cpu_ctx = {}; ///< CPU device
    const int rank;                        ///< MPI rank

    //-------------------------OP--------------------------------------------
    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
};

///@brief allocate the wavefunction
void allocate_psi(Psi<std::complex<double>>*& psi, const int& nks, const std::vector<int>& ngk, const int& nbands, const int& npwx);

} // namespace psi
#endif