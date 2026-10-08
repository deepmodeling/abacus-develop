#ifndef HAMILTPW_NONLOCAL_MATHS_H
#define HAMILTPW_NONLOCAL_MATHS_H

#include "source_base/math_ylmreal.h"
#include "source_base/module_device/device.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/vnl_pw.h"
#include "source_pw/module_pwdft/kernels/stress_op.h"
#include "source_base/kernels/math_kernel_op.h"

namespace hamilt
{

template <typename FPTYPE, typename Device>
class Nonlocal_maths
{
  public:
    Nonlocal_maths(const pseudopot_cell_vnl* nlpp_in, const UnitCell* ucell_in);
    Nonlocal_maths(const ModuleBase::matrix& nhtol, const int lmax, const UnitCell* ucell_in);

  private:
    ModuleBase::matrix nhtol_;
    int lmax_ = 0;
    const UnitCell* ucell_ = nullptr;

    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};

  public:
    // functions
    /**
     * @brief this function prepares all the q (G+k) information in one contiguous memory block
     * including the x, y and z components, its norm and the reciprocal of its norm
     *
     * @param ik index of k point
     * @param pw_basis the plane wave basis
     * @return std::vector<FPTYPE> 1d contiguous memory block containing all the q information. The
     * first 3*npw are data of x, y and z components, the next 2*npw are data of norm and 1/norm.
     * This is beneficial for GPU memory access.
     */
    std::vector<FPTYPE> cal_gk(int ik, const ModulePW::PW_Basis_K* pw_basis);
    /**
     * @brief calculate the real spherical harmonic functions on cpu (and optionally send to gpu,
     * if gpu is available)
     *
     * @param lmax [in] maximum angular momentum to calculate
     * @param npw [in] number of G+k vectors
     * @param gk_in [in] the G+k vectors
     * @param ylm [out] the spherical harmonic functions
     */
    void cal_ylm(int lmax, int npw, const FPTYPE* gk_in, FPTYPE* ylm);
    /// calculate the derivate of the sperical bessel function for projections
    void cal_ylm_deri(int lmax, int npw, const FPTYPE* gk_in, FPTYPE* ylm_deri);
    /// calculate the (-i)^l factors
    std::vector <std::complex<FPTYPE>> cal_pref(int it, const int nh);
    /// calculate the vkb matrix for this atom
    /// vkb = sum_lm (-i)^l * ylm(g^) * vq(g^) * sk(g^)
    void cal_vkb(int it,
                 int ia,
                 int npw,
                 const FPTYPE* vq_in,
                 const FPTYPE* ylm_in,
                 const std::complex<FPTYPE>* sk_in,
                 const std::complex<FPTYPE>* pref_in,
                 std::complex<FPTYPE>* vkb_out);
    /// calculate the dvkb matrix for this atom
    void cal_vkb_deri(int it,
                      int ia,
                      int npw,
                      int ipol,
                      int jpol,
                      const FPTYPE* vq_in,
                      const FPTYPE* vq_deri_in,
                      const FPTYPE* ylm_in,
                      const FPTYPE* ylm_deri_in,
                      const std::complex<FPTYPE>* sk_in,
                      const std::complex<FPTYPE>* pref_in,
                      const FPTYPE* gk_in,
                      std::complex<FPTYPE>* vkb_out);

    /// calculate the ptr used in vkb_op
    void prepare_vkb_ptr(int nbeta,
                         double* nhtol,
                         int nhtol_nc,
                         int npw,
                         int it,
                         std::complex<FPTYPE>* vkb_out,
                         std::complex<FPTYPE>** vkb_ptrs,
                         FPTYPE* ylm_in,
                         FPTYPE** ylm_ptrs,
                         FPTYPE* vq_in,
                         FPTYPE** vq_ptrs);

    /// calculate the indexes used in vkb_deri_op
    /// indexes save (lm, nb, dylm_lm_ipol, dylm_lm_jpol) for nh
    void cal_dvkb_index(const int nbeta,
                        const double* nhtol,
                        const int nhtol_nc,
                        const int npw,
                        const int it,
                        const int ipol,
                        const int jpol,
                        int* indexes);

    static void dylmr2(const int nylm, const int ngy, const FPTYPE* gk, FPTYPE* dylm, const int ipol);
    /// polynomial interpolation tool for calculate derivate of vq
    static FPTYPE Polynomial_Interpolation_nl(const ModuleBase::realArray& table,
                                              const int& dim1,
                                              const int& dim2,
                                              const FPTYPE& table_interval,
                                              const FPTYPE& x);
};

} // namespace hamilt

#endif
