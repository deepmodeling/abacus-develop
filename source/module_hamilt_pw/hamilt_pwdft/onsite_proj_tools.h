#ifndef MODULEHAMILTPW_ONSITEPROJTOOLS_H
#define MODULEHAMILTPW_ONSITEPROJTOOLS_H

#include "module_base/module_device/device.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"
#include "module_base/kernels/math_kernel_op.h"
#include "module_psi/psi.h"

#include <complex>

namespace hamilt
{

/**
 * @brief Nonlocal pseudopotential tools in plane wave basis set.
 * used for calculating force and stress for different algorithm
 * the main functions are:
 * 1. cal_becp: calculate the becp = <psi|beta> for all beta functions
 * 2. cal_dbecp_s: calculate the dbecp_{ij} = <psi|\partial beta/\partial varepsilon_{ij}> for all beta functions
 *                 stress_{ij} = -1/omega \sum_{n,k}f_{nk} \sum_I \sum_{lm,l'm'}D_{l,l'}^{I} becp * dbecp_{ij} also
 * calculated
 * 3. cal_dbecp_f: calculate the dbecp_i = <psi|\partial beta/\partial \tau^I_i> for all beta functions
 * 4. cal_force: calculate the force^I_i = - \sum_{n,k}f_{nk} \sum_{lm,l'm'}D_{l,l'}^{I} becp * dbecp_i
 */
template <typename FPTYPE, typename Device>
class Onsite_Proj_tools
{
  public:
    Onsite_Proj_tools(const pseudopot_cell_vnl* nlpp_in,
                      const UnitCell* ucell_in,
                      const psi::Psi<std::complex<FPTYPE>, Device>* psi_in,
                      const K_Vectors* kv_in,
                      const ModulePW::PW_Basis_K* wfc_basis_in,
                      const Structure_Factor* sf_in,
                      const ModuleBase::matrix& wg,
                      const ModuleBase::matrix& ekb);

    // a more general constructor is in the following
    Onsite_Proj_tools(const std::vector<int>& nproj,     // number of projectors for each atom type
                      const std::vector<int>& lproj,
                      const ModuleBase::realArray& tab,  // radials' spherical bessel transform
                      const ModuleBase::matrix& nhtol,
                      std::complex<FPTYPE>* vkb_buf,
                      const UnitCell* ucell_in,
                      const psi::Psi<std::complex<FPTYPE>, Device>* psi_in,
                      const K_Vectors* kv_in,
                      const ModulePW::PW_Basis_K* wfc_basis_in,
                      const Structure_Factor* sf_in,
                      const ModuleBase::matrix& wg,
                      const ModuleBase::matrix& ekb);

    ~Onsite_Proj_tools();

    /**
     * @brief calculate the becp = <psi|beta> for all beta functions
     */
    void cal_becp(int ik, int npm, std::complex<FPTYPE>* becp_in = nullptr, const std::complex<FPTYPE>* ppsi_in = nullptr);
    /**
     * @brief calculate the dbecp_{ij} = <psi|\partial beta/\partial varepsilon_{ij}> for all beta functions
     *       stress_{ij} = -1/omega \sum_{n,k}f_{nk} \sum_I \sum_{lm,l'm'}D_{l,l'}^{I} becp * dbecp_{ij} also calculated
     */
    void cal_dbecp_s(int ik, int npm, int ipol, int jpol);
    /**
     * @brief calculate the dbecp_i = <psi|\partial beta/\partial \tau^I_i> for all beta functions
     */
    void cal_dbecp_f(int ik, int npm, int ipol);

    void cal_force_dftu(int ik, int npm, FPTYPE* force, const int* orbital_corr, const std::complex<FPTYPE>* vu, const int size_vu, const FPTYPE* h_wg);
    void cal_force_dspin(int ik, int npm, FPTYPE* force, const ModuleBase::Vector3<double>* lambda, const FPTYPE* h_wg);
    void cal_stress_dftu(int ik, int npm, FPTYPE* stress, const int* orbital_corr, const std::complex<FPTYPE>* vu, const int size_vu, const FPTYPE* h_wg);
    void cal_stress_dspin(int ik, int npm, FPTYPE* stress, const ModuleBase::Vector3<double>* lambda, const FPTYPE* h_wg);


    std::complex<FPTYPE>* get_becp() { return becp; }
    std::complex<FPTYPE>* get_dbecp() { return dbecp; }

  private:
    /**
     * @brief allocate the memory for the variables
     */
    void allocate_memory(const ModuleBase::matrix& wg, 
                         const ModuleBase::matrix& ekb,
                         const std::vector<int>& nproj,
                         const std::vector<int>& nch);
    /**
     * @brief delete the memory for the variables
     */
    void delete_memory();

  private:
    /// pointers to access the data without memory arrangement
    const Structure_Factor* sf_;
    const pseudopot_cell_vnl* nlpp_;
    const UnitCell* ucell_;
    const psi::Psi<std::complex<FPTYPE>, Device>* psi_;
    const K_Vectors* kv_;
    const ModulePW::PW_Basis_K* wfc_basis_;

    /// the following variables are used for the calculation
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};
    int nkb;
    int nbands;
    int deeq_dims[4] = {0, 0, 0, 0};    // deeq can be something other than that in pseudopotentials
    int deeq_nc_dims[4] = {0, 0, 0, 0};

    int current_ik = -1;

    int max_nh = 0;
    int max_npw = 0;
    int ntype;
    bool nondiagonal;
    int pre_ik_s = -1;
    int pre_ik_f = -1;

    int* atom_nh = nullptr;
    int* atom_na = nullptr;
    std::vector<int> h_atom_nh;
    std::vector<int> h_atom_na;
    std::vector<int> nproj;

    /// ------------------------- Key optimization -------------------------
    /// @brief the following variables are used for transfer gcar and reuse the values of vkb for force calculation
    int* gcar_zero_indexes = nullptr;
    int gcar_zero_counts[3] = {0, 0, 0};
    std::complex<FPTYPE>* vkb_save = nullptr;
    /// @brief count zero gcar indexes and prepare zero_indexes, do gcar_y /= gcar_x, gcar_z /= gcar_y
    void transfer_gcar(int npw, int npw_max, const FPTYPE* gcar_in);
    /// @brief save the 0-value dvkbs for calculating the dbecp_i in the force calculation
    void save_vkb(int npw, int ipol);
    /// @brief revert the 0-value dvkbs for calculating the dbecp_i in the force calculation
    void revert_vkb(int npw, int ipol);
    /// ---------------------------------------------------------------------

    /// pointers to access the data without memory arrangement
    const ModuleBase::realArray* tabtpr = nullptr;
    const ModuleBase::matrix* nhtol = nullptr;
    int lprojmax = -1;

    FPTYPE* d_wg = nullptr;
    FPTYPE* d_ekb = nullptr;
    FPTYPE* gcar = nullptr;

    FPTYPE* deeq = nullptr;
    std::complex<FPTYPE>* deeq_nc = nullptr;

    FPTYPE* kvec_c = nullptr;
    FPTYPE* qq_nt = nullptr;
    /// --------------------- Key variable ---------------------
    /// borrow the memory from the vkb in VNL_in_pw to calculate vkb and dvkb
    std::complex<FPTYPE>* ppcell_vkb = nullptr;
    /// ---------------------------------------------------------
    /// the following variables are used for the calculation
    /// allocate memory on CPU device only
    std::vector<FPTYPE> g_plus_k;
    /// allocate memory on CPU/GPU device
    FPTYPE* hd_ylm = nullptr;              // (lmax + 1) * (lmax + 1) * npw
    FPTYPE* hd_ylm_deri = nullptr;         // 3 * (lmax + 1) * (lmax + 1) * npw
    FPTYPE* hd_vq = nullptr;               // this->ucell->atoms[it].ncpp.nbeta * npw
    FPTYPE* hd_vq_deri = nullptr;          // this->ucell->atoms[it].ncpp.nbeta * npw
    std::complex<FPTYPE>* hd_sk = nullptr; // this->ucell->nat * npw
    /// allocate global memory on GPU device only
    FPTYPE* d_g_plus_k = nullptr;              // npw * 5
    FPTYPE* d_pref = nullptr;                  // this->ucell->atoms[it].ncpp.nh
    FPTYPE* d_gk = nullptr;                    // this->ucell->atoms[it].ncpp.nh * npw
    FPTYPE* d_vq_tab = nullptr;                // this->ucell->atoms[it].ncpp.nbeta * npw
    std::vector<int> dvkb_indexes;             // this->ucell->atoms[it].ncpp.nh * 4
    int* d_dvkb_indexes = nullptr;             // this->ucell->atoms[it].ncpp.nh * 4
    std::complex<FPTYPE>* d_pref_in = nullptr; // this->ucell->atoms[it].ncpp.nh

    /// becp and dbecp:
    std::complex<FPTYPE>* dbecp = nullptr; // nbands * nkb (for stress) or nbands * nkb * 3 (for force)
    std::complex<FPTYPE>* becp = nullptr;  // nbands * nkb

    /// @brief rename the operators for CPU/GPU device
    using gemm_op = ModuleBase::gemm_op<std::complex<FPTYPE>, Device>;
    using cal_stress_nl_op = hamilt::cal_stress_nl_op<FPTYPE, Device>;
    using cal_dbecp_noevc_nl_op = hamilt::cal_dbecp_noevc_nl_op<FPTYPE, Device>;

    using resmem_complex_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using resmem_complex_h_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>;
    using setmem_complex_op = base_device::memory::set_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_h_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU>;
    using syncmem_complex_h2d_op
        = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2h_op
        = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, base_device::DEVICE_CPU, Device>;

    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using resmem_var_h_op = base_device::memory::resize_memory_op<FPTYPE, base_device::DEVICE_CPU>;
    using setmem_var_op = base_device::memory::set_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using delmem_var_h_op = base_device::memory::delete_memory_op<FPTYPE, base_device::DEVICE_CPU>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<FPTYPE, base_device::DEVICE_CPU, Device>;

    using resmem_int_op = base_device::memory::resize_memory_op<int, Device>;
    using delmem_int_op = base_device::memory::delete_memory_op<int, Device>;
    using syncmem_int_h2d_op = base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>;

    using cal_vq_op = hamilt::cal_vq_op<FPTYPE, Device>;
    using cal_vq_deri_op = hamilt::cal_vq_deri_op<FPTYPE, Device>;

    using cal_vkb_op = hamilt::cal_vkb_op<FPTYPE, Device>;
    using cal_vkb_deri_op = hamilt::cal_vkb_deri_op<FPTYPE, Device>;
};

} // namespace hamilt

#endif