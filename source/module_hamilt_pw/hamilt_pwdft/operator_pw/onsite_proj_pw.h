#ifndef MODULEHAMILTPW_ONSITE_PROJ_PW_H
#define MODULEHAMILTPW_ONSITE_PROJ_PW_H

#include "operator_pw.h"

#include "module_cell/unitcell.h"
#include "module_hsolver/kernels/math_kernel_op.h"

namespace hamilt {

#ifndef ONSITETEMPLATE_H
#define ONSITETEMPLATE_H

template<class T> class OnsiteProj : public T {};
// template<typename Real, typename Device = base_device::DEVICE_CPU>
// class OnsiteProj : public OperatorPW<T, Device> {};

#endif

template<typename T, typename Device>
class OnsiteProj<OperatorPW<T, Device>> : public OperatorPW<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;
  public:
    OnsiteProj(const int* isk_in,
             const UnitCell* ucell_in,
             const bool cal_delta_spin,
             const bool cal_dftu);

    template<typename T_in, typename Device_in = Device>
    explicit OnsiteProj(const OnsiteProj<OperatorPW<T_in, Device_in>>* onsite_proj);

    virtual ~OnsiteProj();

    virtual void init(const int ik_in)override;

    virtual void act(const int nbands,
        const int nbasis,
        const int npol,
        const T* tmpsi_in,
        T* tmhpsi,
        const int ngk = 0,
        const bool is_first_node = false)const override;

    const int *get_isk() const {return this->isk;}
    const UnitCell *get_ucell() const {return this->ucell;}

  private:
    void cal_ps_delta_spin(const int npol, const int m) const;
    void cal_ps_dftu(const int npol, const int m) const;
    void update_becp(const T* psi_in, const int npol, const int m) const;
    void add_onsite_proj(T *hpsi_in, const int npol, const int m) const;

    const int* isk = nullptr;

    const UnitCell* ucell = nullptr;

    mutable int* ip_iat = nullptr;
    mutable T* lambda_coeff = nullptr;
    mutable int* orb_l_iat = nullptr;
    mutable int* ip_m = nullptr;
    mutable int* vu_begin_iat = nullptr;
    mutable T* vu_device = nullptr;

    mutable int nkb_m = 0;

    bool has_delta_spin = false;
    bool has_dftu = false;

    mutable bool init_dftu = false;
    mutable bool init_delta_spin = false;

    mutable T *ps = nullptr;
    int tnp = 0;
    Device* ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};

    using gemv_op = hsolver::gemv_op<T, Device>;
    using gemm_op = hsolver::gemm_op<T, Device>;
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using syncmem_complex_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    using resmem_int_op = base_device::memory::resize_memory_op<int, Device>;
    using resmem_real_op = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_int_op = base_device::memory::delete_memory_op<int, Device>;
    using delmem_real_op = base_device::memory::delete_memory_op<Real, Device>;
    using syncmem_int_h2d_op = base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>;
    using syncmem_real_h2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;

    T one{1, 0};
    T zero{0, 0};
};

} // namespace hamilt

#endif