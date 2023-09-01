#include "module_base/timer.h"
#include "module_hamilt_general/operator.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.h"

using namespace hamilt;

template<typename FPTYPE, typename Device>
OperatorPW<FPTYPE, Device>::~OperatorPW(){};

template<typename FPTYPE, typename Device>
typename OperatorPW<FPTYPE, Device>::hpsi_info OperatorPW<FPTYPE, Device>::hPsi(
    hpsi_info& input) const 
{
  ModuleBase::timer::tick("OperatorPW", "hPsi");
  auto psi_input = std::get<0>(input);
  std::tuple<const std::complex<FPTYPE>*, int> psi_info = psi_input->to_range(std::get<1>(input));
  int nbands = std::get<1>(psi_info);

  std::complex<FPTYPE> *tmhpsi = this->get_hpsi(input);
  const std::complex<FPTYPE> *tmpsi_in = std::get<0>(psi_info);
  //if range in hpsi_info is illegal, the first return of to_range() would be nullptr
  if(tmpsi_in == nullptr)
  {
      ModuleBase::WARNING_QUIT("OperatorPW", "please choose correct range of psi for hPsi()!");
  }

  this->act(nbands, psi_input->get_nbasis(), psi_input->npol, tmpsi_in, tmhpsi, psi_input->get_ngk(this->ik));
  OperatorPW* node((OperatorPW*)this->next_op);
  while(node != nullptr)
  {
      node->act(nbands, psi_input->get_nbasis(), psi_input->npol, tmpsi_in, tmhpsi, psi_input->get_ngk(node->ik));
      node = (OperatorPW*)(node->next_op);
  }

  ModuleBase::timer::tick("OperatorPW", "hPsi");
  
  //if in_place, copy temporary hpsi to target hpsi_pointer, then delete hpsi and new a wrapper for return
  std::complex<FPTYPE>* hpsi_pointer = std::get<2>(input);
  if(this->in_place)
  {
      // ModuleBase::GlobalFunc::COPYARRAY(this->hpsi->get_pointer(), hpsi_pointer, this->hpsi->size());
      syncmem_complex_op()(this->ctx, this->ctx, hpsi_pointer, this->hpsi->get_pointer(), this->hpsi->size());
      delete this->hpsi;
      this->hpsi = new psi::Psi<std::complex<FPTYPE>, Device>(hpsi_pointer, *psi_input, 1, nbands / psi_input->npol);
  }      
  return hpsi_info(this->hpsi, psi::Range(1, 0, 0, nbands / psi_input->npol), hpsi_pointer);
}

namespace hamilt {
template class OperatorPW<float, psi::DEVICE_CPU>;
template class OperatorPW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class OperatorPW<float, psi::DEVICE_GPU>;
template class OperatorPW<double, psi::DEVICE_GPU>;
#endif
}