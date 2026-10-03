#ifndef VELOCITY_PW_H
#define VELOCITY_PW_H
#include "op_pw.h"
#include "source_base/module_device/types.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/projector_gradient.h"
#include "source_pw/module_pwdft/velocity_workspace.h"
#include "source_pw/module_pwdft/vnl_pw.h"
namespace hamilt
{

// velocity operator mv = im/\hbar * [H,r] =  p + im/\hbar [V_NL, r]
template <typename FPTYPE, typename Device = base_device::DEVICE_CPU>
class Velocity
{
  public:
    Velocity(const ModulePW::PW_Basis_K* wfcpw_in,
             const int* isk_in,
             pseudopot_cell_vnl* ppcell_in,
             const UnitCell* ucell_in,
             const bool nonlocal_in = true,
             const typename GetTypeReal<FPTYPE>::type* vtau_in = nullptr,
             const int vtau_col_in = 0,
             const int vtau_row_in = 0);

    ~Velocity();

    void init(const int ik_in);
    /** @brief Refresh momentum and projectors with a Hartree-unit vector potential. */
    void init(const int ik_in, const ModuleBase::Vector3<double>& vector_potential);
    /** @brief Refresh borrowed spin and meta-GGA potential views before reuse. */
    void set_state(const int* isk, const FPTYPE* vtau, const int cols, const int rows)
    {
        this->isk = isk;
        vtau_ = vtau;
        vtau_col_ = cols;
        vtau_row_ = rows;
    }

    /**
     * @brief calculate \hat{v}|\psi>
     *
     * @param psi_in Psi class which contains some information
     * @param n_npwx nbands * NPOL
     * @param tmpsi_in |\psi_i>    size: n_npwx*npwx
     * @param tmvpsi \hat{v}|\psi> size: 3*n_npwx*npwx
     * @param add true : tmvpsi = tmvpsi + v|\psi>  false: tmvpsi = v|\psi>
     *
     */
    void act(const psi::Psi<std::complex<FPTYPE>, Device>* psi_in,
             const int n_npwx,
             const std::complex<FPTYPE>* tmpsi_in,
             std::complex<FPTYPE>* tmvpsi,
             const bool add = false) const;

    bool nonlocal = true;

  private:
    const ModulePW::PW_Basis_K* wfcpw = nullptr;

    const int* isk = nullptr;

    pseudopot_cell_vnl* ppcell = nullptr;

    const UnitCell* ucell = nullptr;

    int ik = 0;

    double tpiba = 0.0;
    const typename GetTypeReal<FPTYPE>::type* vtau_ = nullptr; ///< [CPU] meta-GGA vtau on real grid (nspin x nrxx_smooth)
    int vtau_col_ = 0;                                         ///< number of grid points per spin for vtau
    int vtau_row_ = 0;                                         ///< number of spin channels stored in vtau_
    mutable std::complex<FPTYPE>* porter1_ = nullptr;          ///< workspace on real grid / recip grid
    mutable std::complex<FPTYPE>* porter2_ = nullptr;          ///< workspace on real grid / recip grid
    int momentum_capacity_ = 0;
    int projector_capacity_ = 0;
    mutable int porter_capacity_ = 0;
    ProjectorGradient<FPTYPE, Device> gradient_;
    mutable VelocityWorkspace<FPTYPE, Device> contraction_;
    Device* ctx = {};

  private:
    FPTYPE* gx_ = nullptr;                    ///<[Device, npwx] x component of G+K
    FPTYPE* gy_ = nullptr;                    ///<[Device, npwx] y component of G+K
    FPTYPE* gz_ = nullptr;                    ///<[Device, npwx] z component of G+K
    std::complex<FPTYPE>* vkb_ = nullptr;     ///<[Device, nkb * npwk_max] nonlocal pseudopotential vkb
    std::complex<FPTYPE>* gradvkb_ = nullptr; ///<[Device, 3*nkb * npwk_max] gradient of nonlocal pseudopotential gradvkb

    using Complex = std::complex<FPTYPE>;
    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<FPTYPE, Device, base_device::DEVICE_CPU>;
    using resmem_complex_op = base_device::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
};
} // namespace hamilt
#endif
