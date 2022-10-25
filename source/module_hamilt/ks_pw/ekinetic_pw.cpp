#include "ekinetic_pw.h"

#include "module_base/timer.h"
#include "module_base/tool_quit.h"

using hamilt::Ekinetic;
using hamilt::OperatorPW;

template<typename FPTYPE, typename Device>
Ekinetic<OperatorPW<FPTYPE, Device>>::Ekinetic(
    FPTYPE tpiba2_in,
    const FPTYPE* gk2_in,
    const int gk2_row,
    const int gk2_col)
{
  this->cal_type = pw_ekinetic;
  this->tpiba2 = tpiba2_in;
  this->gk2_row = gk2_row;
  this->gk2_col = gk2_col;
#if ((defined __CUDA) || (defined __ROCM))
  resize_memory_op()(this->ctx, this->gk2, this->gk2_row * this->gk2_col);
  synchronize_memory_op()(this->ctx, this->cpu_ctx, this->gk2, gk2_in, this->gk2_row * this->gk2_col);  
#else
  this->gk2 = gk2_in;
#endif
  if( this->tpiba2 < 1e-10 || this->gk2 == nullptr) {
      ModuleBase::WARNING_QUIT("EkineticPW", "Constuctor of Operator::EkineticPW is failed, please check your code!");
  }
}

template<typename FPTYPE, typename Device>
Ekinetic<OperatorPW<FPTYPE, Device>>::~Ekinetic() {
#if ((defined __CUDA) || (defined __ROCM))
  delete_memory_op()(this->ctx, this->gk2);
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

namespace hamilt{
template class Ekinetic<OperatorPW<double, psi::DEVICE_CPU>>;
#if ((defined __CUDA) || (defined __ROCM))
template class Ekinetic<OperatorPW<double, psi::DEVICE_GPU>>;
#endif
} // namespace hamilt