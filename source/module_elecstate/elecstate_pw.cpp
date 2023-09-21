#include "elecstate_pw.h"

#include "elecstate_getters.h"
#include "module_base/constants.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_psi/kernels/device.h"

namespace elecstate {

template<typename T, typename Device>
ElecStatePW<T, Device>::ElecStatePW(ModulePW::PW_Basis_K *wfc_basis_in, Charge* chg_in, K_Vectors *pkv_in, ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis_Big* bigpw_in) : basis(wfc_basis_in)  
{
    this->classname = "ElecStatePW";
    this->init_ks(chg_in, pkv_in, pkv_in->nks, rhopw_in, bigpw_in);
}

template<typename T, typename Device>
ElecStatePW<T, Device>::~ElecStatePW() 
{
    if (psi::device::get_device_type<Device>(this->ctx) == psi::GpuDevice) {
        delmem_var_op()(this->ctx, this->rho_data);
        if (get_xc_func_type() == 3)
        {
            delmem_var_op()(this->ctx, this->kin_r_data);
        }
    }
    delmem_complex_op()(this->ctx, this->wfcr);
    delmem_complex_op()(this->ctx, this->wfcr_another_spin);
}

template<typename T, typename Device>
void ElecStatePW<T, Device>::init_rho_data() 
{
    if (GlobalV::device_flag == "gpu" || GlobalV::precision_flag == "single") {
        this->rho = new Real*[this->charge->nspin];
        resmem_var_op()(this->ctx, this->rho_data, this->charge->nspin * this->charge->nrxx);
        for (int ii = 0; ii < this->charge->nspin; ii++) {
            this->rho[ii] = this->rho_data + ii * this->charge->nrxx;
        }
        if (get_xc_func_type() == 3)
        {
            this->kin_r = new Real*[this->charge->nspin];
            resmem_var_op()(this->ctx, this->kin_r_data, this->charge->nspin * this->charge->nrxx);
            for (int ii = 0; ii < this->charge->nspin; ii++) {
                this->kin_r[ii] = this->kin_r_data + ii * this->charge->nrxx;
            }
        }
    }
    else {
        this->rho = reinterpret_cast<Real **>(this->charge->rho);
        if (get_xc_func_type() == 3)
        {
            this->kin_r = reinterpret_cast<Real **>(this->charge->kin_r);
        }
    }
    resmem_complex_op()(this->ctx, this->wfcr, this->basis->nmaxgr, "ElecSPW::wfcr");
    resmem_complex_op()(this->ctx, this->wfcr_another_spin, this->charge->nrxx, "ElecSPW::wfcr_a");
    this->init_rho = true;
}

template<typename T, typename Device>
void ElecStatePW<T, Device>::psiToRho(const psi::Psi<T, Device>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "psiToRho");
    ModuleBase::timer::tick("ElecStatePW", "psiToRho");

    if (!this->init_rho) {
        this->init_rho_data();
    }
    this->calculate_weights();

    this->calEBand();

    for(int is=0; is<GlobalV::NSPIN; is++)
	{
        // denghui replaced at 20221110
		// ModuleBase::GlobalFunc::ZEROS(this->rho[is], this->charge->nrxx);
        setmem_var_op()(this->ctx, this->rho[is], 0,  this->charge->nrxx);
        if (get_xc_func_type() == 3)
        {
            // ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
            setmem_var_op()(this->ctx, this->kin_r[is], 0,  this->charge->nrxx);
        }
	}

    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        this->updateRhoK(psi);
    }
    if (GlobalV::device_flag == "gpu" || GlobalV::precision_flag == "single") {
        for (int ii = 0; ii < GlobalV::NSPIN; ii++) {
            castmem_var_d2h_op()(cpu_ctx, this->ctx, this->charge->rho[ii], this->rho[ii], this->charge->nrxx);
            if (get_xc_func_type() == 3)
            {
                castmem_var_d2h_op()(cpu_ctx, this->ctx, this->charge->kin_r[ii], this->kin_r[ii], this->charge->nrxx);
            }
        }
    }
    this->parallelK();
    ModuleBase::timer::tick("ElecStatePW", "psiToRho");
}

template<typename T, typename Device>
void ElecStatePW<T, Device>::updateRhoK(const psi::Psi<T, Device>& psi)
{
    this->rhoBandK(psi);
}

template<typename T, typename Device>
void ElecStatePW<T, Device>::parallelK()
{
#ifdef __MPI
    this->charge->rho_mpi();
    if(GlobalV::ESOLVER_TYPE == "sdft") //qinarui add it 2021-7-21
	{
        this->f_en.eband /= GlobalV::NPROC_IN_POOL;
        MPI_Allreduce(MPI_IN_PLACE, &this->f_en.eband, 1, MPI_DOUBLE, MPI_SUM, STO_WORLD);
    }
#endif
}

template<typename T, typename Device>
void ElecStatePW<T, Device>::rhoBandK(const psi::Psi<T, Device>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "rhoBandK");

    // moved by denghui to constructor at 20221110
    // used for plane wavefunction FFT3D to real space
    // static std::vector<T> wfcr;
    // wfcr.resize(this->basis->nmaxgr);
    // used for plane wavefunction FFT3D to real space, non-collinear spin case
    // static std::vector<std::complex<double>> wfcr_another_spin;
    // if (GlobalV::NSPIN == 4)
    //     wfcr_another_spin.resize(this->charge->nrxx);

    if (!this->init_rho) {
        this->init_rho_data();
    }
    int ik = psi.get_current_k();
    int npw = psi.get_current_nbas();
    int current_spin = 0;
    if (GlobalV::NSPIN == 2)
    {
        current_spin = this->klist->isk[ik];
    }
    int nbands = psi.get_nbands();
    const double threshold = ModuleBase::threshold_wg * this->wg(ik, 0);
    //  here we compute the band energy: the sum of the eigenvalues
    if (GlobalV::NSPIN == 4)
    {
        int npwx = npw / 2;
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            /// be care of when smearing_sigma is large, wg would less than 0
            ///
            if (std::fabs(this->wg(ik, ibnd)) < threshold) {
                continue;
            }

            this->basis->recip_to_real(this->ctx, &psi(ibnd,0), this->wfcr, ik);

            this->basis->recip_to_real(this->ctx, &psi(ibnd,npwx), this->wfcr_another_spin, ik);

            const auto w1 = static_cast<Real>(this->wg(ik, ibnd) / get_ucell_omega());

            // replaced by denghui at 20221110
            elecstate_pw_op()(this->ctx, GlobalV::DOMAG, GlobalV::DOMAG_Z, this->charge->nrxx, w1, this->rho, this->wfcr, this->wfcr_another_spin);
        }
    }
    else
    {
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///
            if (this->wg(ik, ibnd) < threshold) {
                continue;
            }

            this->basis->recip_to_real(this->ctx, &psi(ibnd,0), this->wfcr, ik);

            const auto w1 = static_cast<Real>(this->wg(ik, ibnd) / get_ucell_omega());

            if (w1 != 0.0)
            {
                // replaced by denghui at 20221110
                elecstate_pw_op()(this->ctx,  current_spin, this->charge->nrxx, w1,  this->rho,  this->wfcr);
            }

            // kinetic energy density
            if (get_xc_func_type() == 3)
            {
                for (int j = 0; j < 3; j++)
                {
                    setmem_complex_op()(this->ctx, this->wfcr, 0,  this->charge->nrxx);

                    meta_op()(this->ctx,
                              ik,
                              j,
                              npw,
                              this->basis->npwk_max,
                              static_cast<Real>(get_ucell_tpiba()),
                              this->basis->template get_gcar_data<Real>(),
                              this->basis->template get_kvec_c_data<Real>(),
                              &psi(ibnd, 0),
                              this->wfcr);

                    this->basis->recip_to_real(this->ctx, this->wfcr, this->wfcr, ik);

                    elecstate_pw_op()(this->ctx, current_spin, this->charge->nrxx, w1, this->kin_r, this->wfcr);
                }
            }
        }
    }
}

template class ElecStatePW<std::complex<float>, psi::DEVICE_CPU>;
template class ElecStatePW<std::complex<double>, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ElecStatePW<std::complex<float>, psi::DEVICE_GPU>;
template class ElecStatePW<std::complex<double>, psi::DEVICE_GPU>;
#endif 

} // namespace elecstate