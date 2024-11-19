#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_WFINIT_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_WFINIT_H
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_psi/psi_initializer.h"

namespace psi
{

// This class is used to initialize the wavefunction
template <typename T, typename Device = base_device::DEVICE_CPU>
class WFInit
{
  public:
    WFInit(const std::string& init_wfc_in,
           const std::string& ks_solver_in,
           const std::string& basis_type_in,
           const bool& use_psiinitializer_in,
           wavefunc* p_wf_in,
           ModulePW::PW_Basis_K* pw_wfc_in);
    ~WFInit(){};

    // prepare the wavefunction initialization
    void prepare_init(Structure_Factor* p_sf, //< structure factor
                      UnitCell* p_ucell,      //< unit cell
                      const int& random_seed, //< random seed
#ifdef __MPI
                      Parallel_Kpoints* = nullptr, //< parallel kpoints
                      const int& rank = 0,         //< rank
#endif
                      pseudopot_cell_vnl* = nullptr); //< nonlocal pseudopotential

    // allocate the wavefunction
    void allocate_psi(Psi<std::complex<double>>*& psi, //< psi: wavefunction
                      const int nkstot,                //< total number of k-points for all pools
                      const int nks,                   //< number of k-points in the current pool
                      const int* ngk,                  //< number of G-vectors in the current pool
                      const int npwx,                  //< max number of plane waves of all pools
                      Structure_Factor* p_sf);         //< structure factor

    // make interpolate table
    void make_table(const int nks, Structure_Factor* p_sf);

    //------------------------ only for psi_initializer --------------------
    /**
     * @brief initialize the wavefunction
     *
     * @param psi store the wavefunction
     * @param p_hamilt Hamiltonian operator
     * @param ofs_running output stream for running information
     */
    void initialize_psi(Psi<std::complex<double>>* psi,
                        psi::Psi<T, Device>* kspw_psi,
                        hamilt::Hamilt<T, Device>* p_hamilt,
                        std::ofstream& ofs_running);

    /**
     * @brief get the psi_initializer
     *
     * @return psi_initializer<T, Device>*
     */
    std::weak_ptr<psi::Psi<T, Device>> get_psig() const
    {
        return this->psi_init->share_psig();
    }
    //----------------------------------------------------------------------

  private:
    // psi_initializer<T, Device>* psi_init = nullptr;
    // change to use smart pointer to manage the memory, and avoid memory leak
    // while the std::make_unique() is not supported till C++14,
    // so use the new and std::unique_ptr to manage the memory, but this makes new-delete not symmetric
    std::unique_ptr<psi_initializer<T, Device>> psi_init;
    // temporary
    // wavefunc pointer
    wavefunc* p_wf = nullptr;
    // whether to use psi_initializer
    bool use_psiinitializer = false;
    // wavefunction initialization type
    std::string init_wfc = "none";
    // Kohn-Sham solver type
    std::string ks_solver = "none";
    // basis type
    std::string basis_type = "none";
    // pw basis
    ModulePW::PW_Basis_K* pw_wfc = nullptr;

    Device* ctx = {};
};

} // namespace psi
#endif