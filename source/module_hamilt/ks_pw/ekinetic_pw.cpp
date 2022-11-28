#include "ekinetic_pw.h"

#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_psi/include/device.h"

using hamilt::Ekinetic;
using hamilt::OperatorPW;

template<typename FPTYPE, typename Device>
Ekinetic<OperatorPW<FPTYPE, Device>>::Ekinetic(
    FPTYPE tpiba2_in,
    const FPTYPE* gk2_in,
    const int gk2_row,
    const int gk2_col)
{
  this->classname = "Ekinetic";
  this->cal_type = pw_ekinetic;
  this->tpiba2 = tpiba2_in;
  this->gk2 = gk2_in;
  this->gk2_in = gk2_in;
  this->gk2_row = gk2_row;
  this->gk2_col = gk2_col;
  this->device = psi::device::get_device_type<Device>(this->ctx);
#if ((defined __CUDA) || (defined __ROCM))
  if (this->device == psi::GpuDevice) {
    resmem_var_op()(this->ctx, this->_gk2, this->gk2_row * this->gk2_col);
    syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, this->_gk2, gk2_in, this->gk2_row * this->gk2_col);
    this->gk2 = this->_gk2;
  }
#endif
  if( this->tpiba2 < 1e-10 || this->gk2 == nullptr) {
      ModuleBase::WARNING_QUIT("EkineticPW", "Constuctor of Operator::EkineticPW is failed, please check your code!");
  }
}

template<typename FPTYPE, typename Device>
Ekinetic<OperatorPW<FPTYPE, Device>>::~Ekinetic() {
#if ((defined __CUDA) || (defined __ROCM))
  if (this->device == psi::GpuDevice) {
    delmem_var_op()(this->ctx, this->_gk2);
  }
#endif // __CUDA || __ROCM
}

template<typename FPTYPE, typename Device>
void Ekinetic<OperatorPW<FPTYPE, Device>>::act(
    const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
    const int n_npwx, 
    const std::complex<FPTYPE>* tmpsi_in, 
    std::complex<FPTYPE>* tmhpsi)const
{
  ModuleBase::timer::tick("Operator", "EkineticPW");
  const int npw = psi_in->get_ngk(this->ik);
  this->max_npw = psi_in->get_nbasis() / psi_in->npol;

  const FPTYPE *gk2_ik = &(this->gk2[this->ik * this->gk2_col]);
  // denghui added 20221019
  ekinetic_op()(this->ctx, n_npwx, npw, this->max_npw, tpiba2, gk2_ik, tmhpsi, tmpsi_in);
  // for (int ib = 0; ib < n_npwx; ++ib)
  // {
  //     for (int ig = 0; ig < npw; ++ig)
  //     {
  //         tmhpsi[ig] += gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
  //     }
  //     tmhpsi += this->max_npw;
  //     tmpsi_in += this->max_npw;
  // }
  ModuleBase::timer::tick("Operator", "EkineticPW");
}

// copy construct added by denghui at 20221105
template<typename FPTYPE, typename Device>
template<typename T_in, typename Device_in>
hamilt::Ekinetic<OperatorPW<FPTYPE, Device>>::Ekinetic(const Ekinetic<OperatorPW<T_in, Device_in>> *ekinetic) {
    this->classname = "Ekinetic";
    this->cal_type = pw_ekinetic;
    this->ik = ekinetic->get_ik();
    this->tpiba2 = ekinetic->get_tpiba2();
    this->gk2 = ekinetic->get_gk2();
    this->gk2_row = ekinetic->get_gk2_row();
    this->gk2_col = ekinetic->get_gk2_col();
    this->device = psi::device::get_device_type<Device>(this->ctx);
#if ((defined __CUDA) || (defined __ROCM))
    if (this->device == psi::GpuDevice) {
      resmem_var_op()(this->ctx, this->_gk2, this->gk2_row * this->gk2_col);
      psi::memory::synchronize_memory_op<FPTYPE, Device, Device_in>()(
          this->ctx, ekinetic->get_ctx(),
          this->_gk2, ekinetic->get_gk2(),
          this->gk2_row * this->gk2_col); 
      this->gk2 = this->_gk2;
    }
#endif
    if( this->tpiba2 < 1e-10 || this->gk2 == nullptr) {
        ModuleBase::WARNING_QUIT("EkineticPW", "Copy Constuctor of Operator::EkineticPW is failed, please check your code!");
    }
}

namespace hamilt{
template class Ekinetic<OperatorPW<double, psi::DEVICE_CPU>>;
template Ekinetic<OperatorPW<double, psi::DEVICE_CPU>>::Ekinetic(const Ekinetic<OperatorPW<double, psi::DEVICE_CPU>> *ekinetic);
#if ((defined __CUDA) || (defined __ROCM))
template class Ekinetic<OperatorPW<double, psi::DEVICE_GPU>>;
template Ekinetic<OperatorPW<double, psi::DEVICE_CPU>>::Ekinetic(const Ekinetic<OperatorPW<double, psi::DEVICE_GPU>> *ekinetic);
template Ekinetic<OperatorPW<double, psi::DEVICE_GPU>>::Ekinetic(const Ekinetic<OperatorPW<double, psi::DEVICE_CPU>> *ekinetic);
template Ekinetic<OperatorPW<double, psi::DEVICE_GPU>>::Ekinetic(const Ekinetic<OperatorPW<double, psi::DEVICE_GPU>> *ekinetic);
#endif
} // namespace hamilt