#include "ekinetic_pw.h"

#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_psi/kernels/device.h"


namespace hamilt {

template<typename T, typename Device>
Ekinetic<OperatorPW<T, Device>>::Ekinetic(
    Real tpiba2_in,
    const Real* gk2_in,
    const int gk2_row,
    const int gk2_col)
{
  this->classname = "Ekinetic";
  this->cal_type = pw_ekinetic;
  this->tpiba2 = tpiba2_in;
  this->gk2 = gk2_in;
  this->gk2_row = gk2_row;
  this->gk2_col = gk2_col;
  this->device = psi::device::get_device_type<Device>(this->ctx);
  if( this->tpiba2 < 1e-10 || this->gk2 == nullptr) {
      ModuleBase::WARNING_QUIT("EkineticPW", "Constuctor of Operator::EkineticPW is failed, please check your code!");
  }
}

template<typename T, typename Device>
Ekinetic<OperatorPW<T, Device>>::~Ekinetic() {}

template<typename T, typename Device>
void Ekinetic<OperatorPW<T, Device>>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const T* tmpsi_in,
    T* tmhpsi,
    const int ngk_ik)const
{
    ModuleBase::timer::tick("Operator", "EkineticPW");
    int max_npw = nbasis / npol;

  const Real *gk2_ik = &(this->gk2[this->ik * this->gk2_col]);
  // denghui added 20221019
  ekinetic_op()(this->ctx, nbands, ngk_ik, max_npw, tpiba2, gk2_ik, tmhpsi, tmpsi_in);
  // for (int ib = 0; ib < nbands; ++ib)
  // {
  //     for (int ig = 0; ig < ngk_ik; ++ig)
  //     {
  //         tmhpsi[ig] += gk2_ik[ig] * tpiba2 * tmpsi_in[ig];
  //     }
  //     tmhpsi += max_npw;
  //     tmpsi_in += max_npw;
  // }
  ModuleBase::timer::tick("Operator", "EkineticPW");
}

// copy construct added by denghui at 20221105
template<typename T, typename Device>
template<typename T_in, typename Device_in>
hamilt::Ekinetic<OperatorPW<T, Device>>::Ekinetic(const Ekinetic<OperatorPW<T_in, Device_in>> *ekinetic) {
    this->classname = "Ekinetic";
    this->cal_type = pw_ekinetic;
    this->ik = ekinetic->get_ik();
    this->tpiba2 = ekinetic->get_tpiba2();
    this->gk2 = ekinetic->get_gk2();
    this->gk2_row = ekinetic->get_gk2_row();
    this->gk2_col = ekinetic->get_gk2_col();
    this->device = psi::device::get_device_type<Device>(this->ctx);
    if( this->tpiba2 < 1e-10 || this->gk2 == nullptr) {
        ModuleBase::WARNING_QUIT("EkineticPW", "Copy Constuctor of Operator::EkineticPW is failed, please check your code!");
    }
}

template class Ekinetic<OperatorPW<std::complex<float>, psi::DEVICE_CPU>>;
template class Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_CPU>>;
// template Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_CPU>>::Ekinetic(const Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_CPU>> *ekinetic);
#if ((defined __CUDA) || (defined __ROCM))
template class Ekinetic<OperatorPW<std::complex<float>, psi::DEVICE_GPU>>;
template class Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_GPU>>;
// template Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_CPU>>::Ekinetic(const Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_GPU>> *ekinetic);
// template Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_GPU>>::Ekinetic(const Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_CPU>> *ekinetic);
// template Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_GPU>>::Ekinetic(const Ekinetic<OperatorPW<std::complex<double>, psi::DEVICE_GPU>> *ekinetic);
#endif
} // namespace hamilt