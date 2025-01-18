#include "elecstate_pw.h"
#include "elecstate_getters.h"

namespace elecstate {

template<typename T, typename Device>
void ElecStatePW<T, Device>::cal_tau(const psi::Psi<T, Device>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "cal_tau");
    for(int is=0; is<PARAM.inp.nspin; is++)
	{
        setmem_var_op()(this->kin_r[is], 0,  this->charge->nrxx);
	}

    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        int npw = psi.get_current_ngk();
        int current_spin = 0;
        if (PARAM.inp.nspin == 2)
        {
            current_spin = this->klist->isk[ik];
        }
        int nbands = psi.get_nbands();
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            this->basis->recip_to_real(this->ctx, &psi(ibnd,0), this->wfcr, ik);

            const auto w1 = static_cast<Real>(this->wg(ik, ibnd) / ucell->omega);

            // kinetic energy density
            for (int j = 0; j < 3; j++)
            {
                setmem_complex_op()(this->wfcr, 0,  this->charge->nrxx);

                meta_op()(this->ctx,
                            ik,
                            j,
                            npw,
                            this->basis->npwk_max,
                            static_cast<Real>(ucell->tpiba),
                            this->basis->template get_gcar_data<Real>(),
                            this->basis->template get_kvec_c_data<Real>(),
                            &psi(ibnd, 0),
                            this->wfcr);

                this->basis->recip_to_real(this->ctx, this->wfcr, this->wfcr, ik);

                elecstate_pw_op()(this->ctx, current_spin, this->charge->nrxx, w1, this->kin_r, this->wfcr);
            }
        }
    }
    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single") {
        for (int ii = 0; ii < PARAM.inp.nspin; ii++) {
            castmem_var_d2h_op()(this->charge->kin_r[ii], this->kin_r[ii], this->charge->nrxx);
        }
    }
    this->parallelK();
    ModuleBase::TITLE("ElecStatePW", "cal_tau");
}

template class ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ElecStatePW<std::complex<float>, base_device::DEVICE_GPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>;
#endif 
} // namespace elecstate
