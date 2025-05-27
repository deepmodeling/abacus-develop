#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_KERNELS_FORCE_OP_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_KERNELS_FORCE_OP_H
#include "module_parameter/parameter.h"

#include "module_psi/psi.h"

#include <complex>

namespace hamilt
{

template <typename FPTYPE, typename Device>
struct cal_vkb1_nl_op
{
    /// @brief The prestep to calculate the final forces
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nkb - number of k point
    /// @param npwx - number of planewaves
    /// @param vkb_nc - the second dimension of vkb matrix
    /// @param nbasis - number of planewaves of current k point
    /// @param ipol - 0,1,2
    /// @param NEG_IMAG_UNIT - ModuleBase::NEG_IMAG_UNIT
    /// @param vkb - result of getvnl
    /// @param gcar - GlobalC::wfcpw->gcar
    ///
    /// Output Parameters
    /// @param vkb1 - output vkb matrix
    void operator()(const Device* ctx,
                    const int& nkb,
                    const int& npwx,
                    const int& vkb_nc,
                    const int& nbasis,
                    const int& ipol,
                    const std::complex<FPTYPE>& NEG_IMAG_UNIT,
                    const std::complex<FPTYPE>* vkb,
                    const FPTYPE* gcar,
                    std::complex<FPTYPE>* vkb1);
};

template <typename FPTYPE, typename Device>
struct cal_force_nl_op
{
    /// @brief Calculate the final forces for multi-device
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nondiagonal - control flag
    /// @param nbands_occ - number of occupied bands
    /// @param ntype - total atomic type
    /// @param spin - current spin
    /// @param deeq_2 - the second dimension of deeq
    /// @param deeq_3 - the third dimension of deeq
    /// @param deeq_4 - the forth dimension of deeq
    /// @param forcenl_nc - the second dimension of matrix forcenl
    /// @param nbands - NBANDS
    /// @param nkb - number of k point
    /// @param atom_nh - ucell.atoms[ii].ncpp.nh
    /// @param atom_na - ucell.atoms[ii].na
    /// @param tpiba - ucell.tpiba
    /// @param d_wg - input parameter wg
    /// @param occ - if use the occupation of the bands
    /// @param d_ekb - input parameter ekb
    /// @param qq_nt - ppcell.qq_nt
    /// @param deeq - ppcell.deeq
    /// @param becp - intermediate matrix with PARAM.inp.nbands * nkb
    /// @param dbecp - intermediate matrix with 3 * PARAM.inp.nbands * nkb
    ///
    /// Output Parameters
    /// @param force - output forces
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const bool& nondiagonal,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
    // interface for nspin=4 only
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const std::complex<FPTYPE>* deeq_nc,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
    /// kernel for DFT+U
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const std::complex<FPTYPE>* vu,
                    const int* orbital_corr,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
    /// kernel for DeltaSpin
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const FPTYPE* lambda,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                                          FPTYPE* force);
};

template <typename FPTYPE, typename Device>
struct cal_force_loc_sincos_op
{
    /// @brief Calculate local pseudopotential forces (sincos loop only)
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nat - number of atoms
    /// @param npw - number of plane waves
    /// @param ntype - number of atom types
    /// @param gcar - G-vector Cartesian coordinates [npw * 3]
    /// @param tau - atomic positions [nat * 3]
    /// @param iat2it - atom to type mapping [nat]
    /// @param vloc_factors - precomputed vloc factors [ntype * npw]
    /// @param aux - charge density in G-space [npw]
    ///
    /// Output Parameters
    /// @param force - output forces [nat * 3]
    void operator()(const Device* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ntype,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* vloc_factors,
                    const std::complex<FPTYPE>* aux,
                    FPTYPE* force);
};

template <typename FPTYPE, typename Device>
struct cal_force_ew_parallel_op
{
    /// @brief Calculate Ewald forces (parallel region only)
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nat - number of atoms
    /// @param npw - number of plane waves
    /// @param ntype - number of atom types
    /// @param ig_gge0 - index of G=0 vector (-1 if not present)
    /// @param gcar - G-vector Cartesian coordinates [npw * 3]
    /// @param tau - atomic positions [nat * 3]
    /// @param iat2it - atom to type mapping [nat]
    /// @param it_facts - precomputed it_fact for each atom [nat]
    /// @param atom_na - number of atoms of each type [ntype]
    /// @param zv_values - valence charge for each type [ntype]
    /// @param aux - structure factor related array [npw]
    /// @param alpha - Ewald parameter
    /// @param fact - scaling factor
    /// @param lat0 - lattice constant
    /// @param latvec - lattice vectors [3 * 3]
    /// @param G - reciprocal lattice vectors [3 * 3]
    /// @param use_pwa - whether PAW is used
    /// Output Parameters
    /// @param force - output forces [nat * 3]
    void operator()(const Device* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ntype,
                    const int& ig_gge0,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* it_facts,
                    const int* atom_na,
                    const FPTYPE* zv_values,
                    const std::complex<FPTYPE>* aux,
                    const FPTYPE& alpha,
                    const FPTYPE& fact,
                    const FPTYPE& lat0,
                    const FPTYPE* latvec,
                    const FPTYPE* G,
                    const bool& use_paw,
                    FPTYPE* force);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_vkb1_nl_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nkb,
                    const int& npwx,
                    const int& vkb_nc,
                    const int& nbasis,
                    const int& ipol,
                    const std::complex<FPTYPE>& NEG_IMAG_UNIT,
                    const std::complex<FPTYPE>* vkb,
                    const FPTYPE* gcar,
                    std::complex<FPTYPE>* vkb1);
};

template <typename FPTYPE>
struct cal_force_nl_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const bool& nondiagonal,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
    // interface for nspin=4 only
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const std::complex<FPTYPE>* deeq_nc,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
    /// kernel for DFT+U
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const std::complex<FPTYPE>* vu,
                    const int* orbital_corr,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
    /// kernel for DeltaSpin
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const FPTYPE* lambda,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force);
};

template <typename FPTYPE>
struct cal_force_loc_sincos_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ntype,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* vloc_factors,
                    const std::complex<FPTYPE>* aux,
                    FPTYPE* force);
};

template <typename FPTYPE>
struct cal_force_ew_parallel_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ig_gge0,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* it_facts,
                    const std::complex<FPTYPE>* aux,
                    FPTYPE* force);
};

/**
 * @brief revert the vkb values for force_nl calculation
 */
template <typename FPTYPE>
void revertVkbValues(const int* gcar_zero_ptrs,
                     std::complex<FPTYPE>* vkb_ptr,
                     const std::complex<FPTYPE>* vkb_save_ptr,
                     int nkb,
                     int gcar_zero_counts,
                     int npw,
                     int ipol,
                     int npwx,
                     const std::complex<FPTYPE> coeff);

/**
 * @brief save the vkb values for force_nl calculation
 */
template <typename FPTYPE>
void saveVkbValues(const int* gcar_zero_ptrs,
                   const std::complex<FPTYPE>* vkb_ptr,
                   std::complex<FPTYPE>* vkb_save_ptr,
                   int nkb,
                   int gcar_zero_counts,
                   int npw,
                   int ipol,
                   int npwx);

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
} // namespace hamilt
#endif // W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_PW_HAMILT_PWDFT_KERNELS_FORCE_OP_H