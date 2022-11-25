#include "elecstate_pw.h"

#include "module_base/constants.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pw/global.h"
#include "module_base/timer.h"
#include "module_psi/include/device.h"

namespace elecstate
{

template<typename FPTYPE, typename Device>
ElecStatePW<FPTYPE, Device>::ElecStatePW(ModulePW::PW_Basis_K *wfc_basis_in, Charge* chg_in, K_Vectors *pkv_in) : basis(wfc_basis_in)  
{
    this->classname = "ElecStatePW";
    this->init_ks(chg_in, pkv_in, pkv_in->nks);
}

template<typename FPTYPE, typename Device>
ElecStatePW<FPTYPE, Device>::~ElecStatePW() 
{
    if (psi::device::get_device_type<Device>(this->ctx) == psi::GpuDevice) {
        delmem_var_op()(this->ctx, this->rho_data);
        if (XC_Functional::get_func_type() == 3) {
            delmem_var_op()(this->ctx, this->kin_r_data);
        }
    }
    delmem_complex_op()(this->ctx, this->wfcr);
    delmem_complex_op()(this->ctx, this->wfcr_another_spin);
}

template<typename FPTYPE, typename Device>
void ElecStatePW<FPTYPE, Device>::init_rho_data() 
{
    if (psi::device::get_device_type<Device>(this->ctx) == psi::GpuDevice) {
        this->rho = new FPTYPE*[this->charge->nspin];
        resmem_var_op()(this->ctx, this->rho_data, this->charge->nspin * this->charge->nrxx);
        for (int ii = 0; ii < this->charge->nspin; ii++) {
            this->rho[ii] = this->rho_data + ii * this->charge->nrxx;
        }
        if (XC_Functional::get_func_type() == 3) {
            this->kin_r = new FPTYPE*[this->charge->nspin];
            resmem_var_op()(this->ctx, this->kin_r_data, this->charge->nspin * this->charge->nrxx);
            for (int ii = 0; ii < this->charge->nspin; ii++) {
                this->kin_r[ii] = this->kin_r_data + ii * this->charge->nrxx;
            }
        }
    }
    else {
        this->rho = this->charge->rho;
        if (XC_Functional::get_func_type() == 3) {
            this->kin_r = this->charge->kin_r;
        }
    }
    resmem_complex_op()(this->ctx, this->wfcr, this->basis->nmaxgr);
    resmem_complex_op()(this->ctx, this->wfcr_another_spin, this->charge->nrxx);
    this->init_rho = true;
}

template<typename FPTYPE, typename Device>
void ElecStatePW<FPTYPE, Device>::psiToRho(const psi::Psi<std::complex<FPTYPE>, Device>& psi)
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
		if (XC_Functional::get_func_type() == 3)
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
    if (psi::device::get_device_type<Device>(this->ctx) == psi::GpuDevice) {
        for (int ii = 0; ii < GlobalV::NSPIN; ii++) {
            syncmem_var_d2h_op()(
                this->cpu_ctx, 
                this->ctx,
                this->charge->rho[ii],
                this->rho[ii],
                this->charge->nrxx);
        }
    }
    this->parallelK();
    ModuleBase::timer::tick("ElecStatePW", "psiToRho");
}

template<typename FPTYPE, typename Device>
void ElecStatePW<FPTYPE, Device>::updateRhoK(const psi::Psi<std::complex<FPTYPE>, Device>& psi)
{
    this->rhoBandK(psi);
}

template<typename FPTYPE, typename Device>
void ElecStatePW<FPTYPE, Device>::parallelK()
{
#ifdef __MPI
    this->charge->rho_mpi();
    if(GlobalV::ESOLVER_TYPE == "sdft") //qinarui add it 2021-7-21
	{
		this->eband /= GlobalV::NPROC_IN_POOL;
		MPI_Allreduce(MPI_IN_PLACE, &this->eband, 1, MPI_DOUBLE, MPI_SUM , STO_WORLD);
	}
#endif
}

template<typename FPTYPE, typename Device>
void ElecStatePW<FPTYPE, Device>::rhoBandK(const psi::Psi<std::complex<FPTYPE>, Device>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "rhoBandK");

    // moved by denghui to constructor at 20221110
    // used for plane wavefunction FFT3D to real space
    // static std::vector<std::complex<FPTYPE>> wfcr;
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
    //  here we compute the band energy: the sum of the eigenvalues
    if (GlobalV::NSPIN == 4)
    {
        int npwx = npw / 2;
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///
            if (this->wg(ik, ibnd) < ModuleBase::threshold_wg) {
                continue;
            }

            this->basis->recip_to_real(this->ctx, &psi(ibnd,0), this->wfcr, ik);

            this->basis->recip_to_real(this->ctx, &psi(ibnd,npwx), this->wfcr_another_spin, ik);

            const double w1 = this->wg(ik, ibnd) / GlobalC::ucell.omega;

            // replaced by denghui at 20221110
            elecstate_pw_op()(
                this->ctx, 
                GlobalV::DOMAG, GlobalV::DOMAG_Z, 
                this->charge->nrxx, w1, 
                this->rho, 
                this->wfcr, 
                this->wfcr_another_spin);
            // Increment the charge density in chr.rho for real space
            // for (int ir = 0; ir < this->charge->nrxx; ir++)
            // {
            //     this->charge->rho[0][ir] += w1 * (norm(wfcr[ir]) + norm(wfcr_another_spin[ir]));
            // }
            // // In this case, calculate the three components of the magnetization
            // if (GlobalV::DOMAG)
            // {
            //     for (int ir = 0; ir < this->charge->nrxx; ir++)
            //     {
            //         this->charge->rho[1][ir] += w1 * 2.0
            //                                     * (wfcr[ir].real() * wfcr_another_spin[ir].real()
            //                                        + wfcr[ir].imag() * wfcr_another_spin[ir].imag());
            //         this->charge->rho[2][ir] += w1 * 2.0
            //                                     * (wfcr[ir].real() * wfcr_another_spin[ir].imag()
            //                                        - wfcr_another_spin[ir].real() * wfcr[ir].imag());
            //         this->charge->rho[3][ir] += w1 * (norm(wfcr[ir]) - norm(wfcr_another_spin[ir]));
            //     }
            // }
            // else if (GlobalV::DOMAG_Z)
            // {
            //     for (int ir = 0; ir < this->charge->nrxx; ir++)
            //     {
            //         this->charge->rho[1][ir] = 0;
            //         this->charge->rho[2][ir] = 0;
            //         this->charge->rho[3][ir] += w1 * (norm(wfcr[ir]) - norm(wfcr_another_spin[ir]));
            //     }
            // }
            // else
            //     for (int is = 1; is < 4; is++)
            //     {
            //         for (int ir = 0; ir < this->charge->nrxx; ir++)
            //             this->charge->rho[is][ir] = 0;
            //     }
        }
    }
    else
    {
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///
            if (this->wg(ik, ibnd) < ModuleBase::threshold_wg) {
                continue;
            }

            this->basis->recip_to_real(this->ctx, &psi(ibnd,0), this->wfcr, ik);

            const FPTYPE w1 = this->wg(ik, ibnd) / GlobalC::ucell.omega;

            if (w1 != 0.0)
            {
                // replaced by denghui at 20221110
                elecstate_pw_op()(
                    this->ctx, 
                    current_spin, this->charge->nrxx, 
                    w1, 
                    this->rho, 
                    this->wfcr);
                // for (int ir = 0; ir < this->charge->nrxx; ir++)
                // {
                //     this->rho[current_spin][ir] += w1 * norm(wfcr[ir]);
                // }
            }

            // kinetic energy density
            if (XC_Functional::get_func_type() == 3)
            {
                for (int j = 0; j < 3; j++)
                {
                    ModuleBase::GlobalFunc::ZEROS(this->wfcr, this->charge->nrxx);
                    for (int ig = 0; ig < npw; ig++)
                    {
                        double fact
                            = this->basis->getgpluskcar(ik, ig)[j] * GlobalC::ucell.tpiba;
                        wfcr[ig] = psi(ibnd, ig) * complex<double>(0.0, fact);
                    }

                    this->basis->recip2real(this->wfcr, this->wfcr, ik);
                    
                    for (int ir = 0; ir < this->charge->nrxx; ir++)
                    {
                        this->kin_r[current_spin][ir] += w1 * norm(wfcr[ir]);
                    }
                }
            }
        }
    }
}

template class ElecStatePW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ElecStatePW<double, psi::DEVICE_GPU>;
#endif 

} // namespace elecstate