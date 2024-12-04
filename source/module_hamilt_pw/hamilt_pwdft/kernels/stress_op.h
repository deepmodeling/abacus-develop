#ifndef SRC_PW_STRESS_MULTI_DEVICE_H
#define SRC_PW_STRESS_MULTI_DEVICE_H
#include "module_parameter/parameter.h"

#include "module_psi/psi.h"

#include <complex>
#include <module_base/macros.h>

namespace hamilt
{

template <typename FPTYPE, typename Device>
struct cal_dbecp_noevc_nl_op
{
    /// @brief The prestep to calculate the final stresses
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param ipol - loop of 0, 1, 2
    /// @param jpol - loop of 0, 1, 2
    /// @param nkb - number of k points
    /// @param npw - number of planewaves
    /// @param npwx - max number of planewaves
    /// @param ik - the current k point
    /// @param tpiba - ucell.tpiba
    /// @param gcar - GlobalC::wfcpw->gcar
    /// @param kvec_c - GlobalC::wfcpw->kvec_c
    /// @param vkbi - _vkb0[ipol]
    /// @param vkbj - _vkb0[jpol]
    /// @param vkb - result of getvnl
    /// @param vkb1 - intermeidate matrix with size nkb * GlobalC::wf.npwx
    /// @param vkb2 - intermeidate matrix with size nkb * GlobalC::wf.npwx
    ///
    /// Output Parameters
    /// @param dbecp_noevc - result matrix with size nkb * GlobalC::wf.npwx
    void operator()(const Device* ctx,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& npw,
                    const int& npwx,
                    const int& ik,
                    const FPTYPE& tpiba,
                    const FPTYPE* gcar,
                    const FPTYPE* kvec_c,
                    std::complex<FPTYPE>* vkbi,
                    std::complex<FPTYPE>* vkbj,
                    std::complex<FPTYPE>* vkb,
                    std::complex<FPTYPE>* vkb1,
                    std::complex<FPTYPE>* vkb2,
                    std::complex<FPTYPE>* dbecp_noevc);
};

template <typename FPTYPE, typename Device>
struct cal_stress_nl_op
{
    /// @brief Calculate the final stresses for multi-device
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nondiagonal - control flag
    /// @param ipol - loop of 0, 1, 2
    /// @param jpol - loop of 0, 1, 2
    /// @param nkb - number of k point
    /// @param nbands_occ - number of occupied bands
    /// @param ntype - total atomic type
    /// @param spin - current spin
    /// @param deeq_2 - the second dimension of deeq
    /// @param deeq_3 - the third dimension of deeq
    /// @param deeq_4 - the forth dimension of deeq
    /// @param atom_nh - ucell.atoms[ii].ncpp.nh
    /// @param atom_na - ucell.atoms[ii].na
    /// @param d_wg - input parameter wg
    /// @param occ - if use the occupation of the bands
    /// @param d_ekb - input parameter ekb
    /// @param qq_nt - ppcell.qq_nt
    /// @param deeq - ppcell.deeq
    /// @param becp - intermediate matrix with PARAM.inp.nbands * nkb
    /// @param dbecp - intermediate matrix with 3 * PARAM.inp.nbands * nkb
    ///
    /// Output Parameters
    /// @param stress - output stresses
    void operator()(const Device* ctx,
                    const bool& nondiagonal,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* stress);
    // interface for nspin=4 only
    void operator()(const Device* ctx,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const std::complex<FPTYPE>* deeq_nc,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* stress);
};

template <typename T, typename Device>
struct cal_stress_mgga_op
{
    using Real = typename GetTypeReal<T>::type;
    void operator()(const int& spin, const int& nrxx, const Real& w1, const T* gradwfc, Real* crosstaus);
};

// cpu version first, gpu version later
template <typename FPTYPE, typename Device>
struct cal_vkb_op
{
    void operator()(const Device* ctx,
                    const int nh,
                    const int npw,
                    const int* indexes,
                    const FPTYPE* vqs_in,
                    const FPTYPE* ylms_in,
                    const std::complex<FPTYPE>* sk_in,
                    const std::complex<FPTYPE>* pref_in,
                    std::complex<FPTYPE>* vkbs_out);
};

// cpu version first, gpu version later
template <typename FPTYPE, typename Device>
struct cal_vkb_deri_op
{
    void operator()(const Device* ctx,
                    const int nh,
                    const int npw,
                    const int ipol,
                    const int jpol,
                    const int* indexes,
                    const FPTYPE* vqs_in,
                    const FPTYPE* vqs_deri_in,
                    const FPTYPE* ylms_in,
                    const FPTYPE* ylms_deri_in,
                    const std::complex<FPTYPE>* sk_in,
                    const std::complex<FPTYPE>* pref_in,
                    const FPTYPE* gk_in,
                    std::complex<FPTYPE>* vkbs_out);
};

// cpu version first, gpu version later
template <typename FPTYPE, typename Device>
struct cal_vq_op
{
    void operator()(const Device* ctx,
                    const FPTYPE* tab,
                    int it,
                    const FPTYPE* gk,
                    int npw,
                    const int tab_2,
                    const int tab_3,
                    const FPTYPE table_interval,
                    const int nbeta,
                    FPTYPE* vq);
};

// cpu version first, gpu version later
template <typename FPTYPE, typename Device>
struct cal_vq_deri_op
{
    void operator()(const Device* ctx,
                    const FPTYPE* tab,
                    int it,
                    const FPTYPE* gk,
                    int npw,
                    const int tab_2,
                    const int tab_3,
                    const FPTYPE table_interval,
                    const int nbeta,
                    FPTYPE* vq);
};


template <typename FPTYPE, typename Device>
struct cal_stress_drhoc_aux_op{
    void operator()(
        const FPTYPE* r, const FPTYPE* rhoc, 
        const FPTYPE *gx_arr, const FPTYPE *rab, FPTYPE *drhocg, 
        const int mesh, const int igl0, const int ngg, const double omega,
        int type
    );
};

template <typename FPTYPE, typename Device>
struct cal_force_npw_op{
    void operator()(const std::complex<FPTYPE> *psiv,
                    const FPTYPE* gv_x, const FPTYPE* gv_y, const FPTYPE* gv_z,
                    const FPTYPE* rhocgigg_vec,
                    FPTYPE* force,
                    const FPTYPE pos_x, const FPTYPE pos_y, const FPTYPE pos_xz,
                    const int npw,
                    const FPTYPE omega, const FPTYPE tpiba
    );
};

template <typename FPTYPE, typename Device>
struct cal_multi_dot_op{
    FPTYPE operator()(const int& npw,
                    const FPTYPE& fac,
                    const FPTYPE* gk1,
                    const FPTYPE* gk2,
                    const FPTYPE* d_kfac,
                    const std::complex<FPTYPE>* psi);
};


#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
template <typename FPTYPE>
struct cal_dbecp_noevc_nl_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& npw,
                    const int& npwx,
                    const int& ik,
                    const FPTYPE& tpiba,
                    const FPTYPE* gcar,
                    const FPTYPE* kvec_c,
                    std::complex<FPTYPE>* vkbi,
                    std::complex<FPTYPE>* vkbj,
                    std::complex<FPTYPE>* vkb,
                    std::complex<FPTYPE>* vkb1,
                    std::complex<FPTYPE>* vkb2,
                    std::complex<FPTYPE>* dbecp_noevc);
};

template <typename FPTYPE>
struct cal_stress_nl_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const bool& nondiagonal,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* stress);
    // interface for nspin=4 only
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const std::complex<FPTYPE>* deeq_nc,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* stress);
};

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vkb_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int nh,
                    const int npw,
                    const int* indexes,
                    const FPTYPE* vqs_in,
                    const FPTYPE* ylms_in,
                    const std::complex<FPTYPE>* sk_in,
                    const std::complex<FPTYPE>* pref_in,
                    std::complex<FPTYPE>* vkbs_out);
};

template <typename FPTYPE>
struct cal_vkb_deri_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int nh,
                    const int npw,
                    const int ipol,
                    const int jpol,
                    const int* indexes,
                    const FPTYPE* vqs_in,
                    const FPTYPE* vqs_deri_in,
                    const FPTYPE* ylms_in,
                    const FPTYPE* ylms_deri_in,
                    const std::complex<FPTYPE>* sk_in,
                    const std::complex<FPTYPE>* pref_in,
                    const FPTYPE* gk_in,
                    std::complex<FPTYPE>* vkbs_out);
};

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vq_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const FPTYPE* tab,
                    int it,
                    const FPTYPE* gk,
                    int npw,
                    const int tab_2,
                    const int tab_3,
                    const FPTYPE table_interval,
                    const int nbeta,
                    FPTYPE* vq);
};

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vq_deri_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const FPTYPE* tab,
                    int it,
                    const FPTYPE* gk,
                    int npw,
                    const int tab_2,
                    const int tab_3,
                    const FPTYPE table_interval,
                    const int nbeta,
                    FPTYPE* vq);
};

template <typename FPTYPE>
struct cal_multi_dot_op<FPTYPE, base_device::DEVICE_GPU>{
    FPTYPE operator()(const int& npw,
                    const FPTYPE& fac,
                    const FPTYPE* gk1,
                    const FPTYPE* gk2,
                    const FPTYPE* d_kfac,
                    const std::complex<FPTYPE>* psi);
};

/**
 * The operator is used to compute the auxiliary amount of stress /force 
 * in parallel on the GPU. They identify type with the type provided and 
 * select different calculation methods,
 *
 * The function is called by the module as follows
 *      Type = 0 -> stress_cc
 *      Type = 1 -> stress_cc, force_cc
 *      Type = 2 -> force_scc
 *      Type = 3 -> stress_loc
 *
 *  Int the function aux is obtained by traversing the `ngg` and `mesh` firstly,
 *  and then aux is processed by Simpson integral method to obtain auxiliary 
 *  quantities drhocg.
 *
 * In the GPU operator, temporary array space of mesh size is required in order 
 * not to apply Simpson interpolation (which causes GPU memory overflow). 
 * The Simpson integral is then reconstructed in the loop body of the mesh, 
 * using the Simpson integral computed in the loop, rather than executed once 
 * after the loop. After that, in order to reduce the if condition judgment brought 
 * by Simpson interpolation in the loop body, lambda expression is used to shift the 
 * boundary condition out.
 */
template <typename FPTYPE>
struct cal_stress_drhoc_aux_op<FPTYPE, base_device::DEVICE_GPU>{
    void operator()(
        const FPTYPE* r, const FPTYPE* rhoc, 
        const FPTYPE *gx_arr, const FPTYPE *rab, FPTYPE *drhocg, 
        const int mesh, const int igl0, const int ngg, const double omega,
        int type
    );
};


/**
 * This operator is used to compute the force force in three directions for each atom in force_cc 
 * in parallel on GPU [0~3], which is:
 * Force_p =    (2* pi * tpiba * omega * rhocg[ig] * gv_p[ig] 
 *              * (gv_x[ig] * pos_x + gv_y[ig] * pos_y + gv_z[ig] * pos_z)
 *              * complex(sinp, cosp) * psiv[ig]).real()
 *
 * The operator splits NPW into blocks on the GPU in parallel, and the block size is t_size = 1024.
 */
template <typename FPTYPE>
struct cal_force_npw_op<FPTYPE, base_device::DEVICE_GPU>{
    void operator()(const std::complex<FPTYPE> *psiv,
                    const FPTYPE* gv_x, const FPTYPE* gv_y, const FPTYPE* gv_z,
                    const FPTYPE* rhocgigg_vec,
                    FPTYPE* force,
                    const FPTYPE pos_x, const FPTYPE pos_y, const FPTYPE pos_xz,
                    const int npw,
                    const FPTYPE omega, const FPTYPE tpiba
    );
};



#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

template <typename Device>
struct pointer_array_malloc
{
    void operator()(void** ptr, const int n);
};

template <typename Device>
struct synchronize_ptrs
{
    void operator()(void** ptr_out, const void** ptr_in, const int size);
};

} // namespace hamilt
#endif // SRC_PW_STRESS_MULTI_DEVICE_H