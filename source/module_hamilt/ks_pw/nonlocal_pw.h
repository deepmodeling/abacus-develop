#ifndef NONLOCALPW_H
#define NONLOCALPW_H

#include "operator_pw.h"

#include "module_cell/unitcell.h"
#include "module_hamilt/include/nonlocal.h"
#include "module_hsolver/include/math_kernel.h"

#include "src_pw/VNL_in_pw.h"

namespace hamilt {

#ifndef NONLOCALTEMPLATE_H
#define NONLOCALTEMPLATE_H

template<class T> class Nonlocal : public T
{};
// template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
// class Nonlocal : public OperatorPW<FPTYPE, Device> {};

#endif

template<typename FPTYPE, typename Device>
class Nonlocal<OperatorPW<FPTYPE, Device>> : public OperatorPW<FPTYPE, Device>
{
  public:
    Nonlocal(const int* isk_in,const pseudopot_cell_vnl* ppcell_in,const UnitCell* ucell_in);

    template<typename T_in, typename Device_in = Device>
    explicit Nonlocal(const Nonlocal<OperatorPW<T_in, Device_in>>* nonlocal);

    virtual ~Nonlocal();

    virtual void init(const int ik_in)override;

    virtual void act(
        const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
        const int n_npwx, 
        const std::complex<FPTYPE>* tmpsi_in, 
        std::complex<FPTYPE>* tmhpsi
    )const override;

    const int *get_isk() const {return this->isk;}
    const pseudopot_cell_vnl *get_ppcell() const {return this->ppcell;}
    const UnitCell *get_ucell() const {return this->ucell;}

  private:
    void add_nonlocal_pp(std::complex<FPTYPE> *hpsi_in, const std::complex<FPTYPE> *becp, const int m) const;

    mutable int max_npw = 0;

    mutable int npw = 0;

    mutable int npol = 0;

    const int* isk = nullptr;

    const pseudopot_cell_vnl* ppcell = nullptr;

    const UnitCell* ucell = nullptr;

    mutable std::complex<FPTYPE> *ps = nullptr;
    mutable std::complex<FPTYPE> *vkb = nullptr;
    mutable std::complex<FPTYPE> *becp = nullptr;
    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};
    FPTYPE * deeq = nullptr;
    std::complex<FPTYPE> * deeq_nc = nullptr;
    // using nonlocal_op = nonlocal_pw_op<FPTYPE, Device>;
    using gemv_op = hsolver::gemv_op<FPTYPE, Device>;
    using gemm_op = hsolver::gemm_op<FPTYPE, Device>;
    using nonlocal_op = nonlocal_pw_op<FPTYPE, Device>;
    using setmem_complex_op = psi::memory::set_memory_op<std::complex<FPTYPE>, Device>;
    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using syncmem_complex_h2d_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, psi::DEVICE_CPU>;
};

} // namespace hamilt

#endif