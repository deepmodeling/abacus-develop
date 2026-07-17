#include "elecstate_pw.h"

#include "source_base/constants.h"
#include "source_base/libm/libm.h"
#include "source_base/math_ylmreal.h"
#include "source_base/module_device/device.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/vnl_pw.h"

namespace elecstate {

template <typename T, typename Device>
ElecStatePW<T, Device>::ElecStatePW(ModulePW::PW_Basis_K* wfc_basis_in,
                                    Charge* chr_in,
                                    K_Vectors* pkv_in,
                                    UnitCell* ucell_in,
                                    pseudopot_cell_vnl* ppcell_in,
                                    ModulePW::PW_Basis* rhopw_in,
                                    ModulePW::PW_Basis_Big* bigpw_in)
    : basis(wfc_basis_in)
{
    this->classname = "ElecStatePW";
    this->rhopw_smooth = rhopw_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->init_ks(chr_in, pkv_in, pkv_in->get_nks(), bigpw_in);
}

template<typename T, typename Device>
ElecStatePW<T, Device>::~ElecStatePW() 
{
    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        delmem_var_op()(this->rho_data);
        delete[] this->rho;

        if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
        {
            delmem_complex_op()(this->rhog_data);
            delete[] this->rhog;
        }
        if (XC_Functional::get_func_type() == 3 || PARAM.inp.out_elf[0] > 0)
        {
            delmem_var_op()(this->kin_r_data);
            delete[] this->kin_r;
        }
    }
    if (PARAM.globalv.use_uspp)
    {
        delmem_var_h_op()(this->becsum);
    }
    delmem_complex_op()(this->wfcr);
    delmem_complex_op()(this->wfcr_another_spin);
}

template<typename T, typename Device>
void ElecStatePW<T, Device>::init_rho_data() 
{
    if (this->init_rho)
    {
        return;
    }

    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        this->rho = new Real*[this->charge->nspin];
        resmem_var_op()(this->rho_data, this->charge->nspin * this->charge->nrxx);
        for (int ii = 0; ii < this->charge->nspin; ii++)
        {
            this->rho[ii] = this->rho_data + ii * this->charge->nrxx;
        }
        if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
        {
            this->rhog = new T*[this->charge->nspin];
            resmem_complex_op()(this->rhog_data, this->charge->nspin * this->charge->rhopw->npw);
            for (int ii = 0; ii < this->charge->nspin; ii++)
            {
                this->rhog[ii] = this->rhog_data + ii * this->charge->rhopw->npw;
            }
        }
        if (XC_Functional::get_func_type() == 3 || PARAM.inp.out_elf[0] > 0)
        {
            this->kin_r = new Real*[this->charge->nspin];
            resmem_var_op()(this->kin_r_data, this->charge->nspin * this->charge->nrxx);
            for (int ii = 0; ii < this->charge->nspin; ii++) {
                this->kin_r[ii] = this->kin_r_data + ii * this->charge->nrxx;
            }
        }
    }
    else
    {
        this->rho = reinterpret_cast<Real **>(this->charge->rho);
        if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
        {
            this->rhog = reinterpret_cast<T**>(this->charge->rhog);
        }
        if (XC_Functional::get_func_type() == 3 || PARAM.inp.out_elf[0] > 0)
        {
            this->kin_r = reinterpret_cast<Real **>(this->charge->kin_r);
        }
    }
    resmem_complex_op()(this->wfcr, this->basis->nmaxgr, "ElecSPW::wfcr");
    resmem_complex_op()(this->wfcr_another_spin, this->basis->nrxx, "ElecSPW::wfcr_a");
    this->init_rho = true;
}

template<typename T, typename Device>
void ElecStatePW<T, Device>::psiToRho(const psi::Psi<T, Device>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "psiToRho");
    ModuleBase::timer::start("ElecStatePW", "psiToRho");

    this->init_rho_data();

    for(int is=0; is<PARAM.inp.nspin; is++)
	{
        // denghui replaced at 20221110
		// ModuleBase::GlobalFunc::ZEROS(this->rho[is], this->charge->nrxx);
        setmem_var_op()(this->rho[is], 0,  this->charge->nrxx);
        if (XC_Functional::get_func_type() == 3)
        {
            // ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
            setmem_var_op()(this->kin_r[is], 0,  this->charge->nrxx);
        }
        if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
        {
            setmem_complex_op()(this->rhog[is], 0, this->charge->rhopw->npw);
        }
    }

    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        this->updateRhoK(psi);
    }

    this->add_usrho(psi);

    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        for (int ii = 0; ii < PARAM.inp.nspin; ii++)
        {
            castmem_var_d2h_op()(this->charge->rho[ii], this->rho[ii], this->charge->nrxx);
            if (XC_Functional::get_func_type() == 3)
            {
                castmem_var_d2h_op()(this->charge->kin_r[ii], this->kin_r[ii], this->charge->nrxx);
            }
        }
    }
    this->parallelK();
    ModuleBase::timer::end("ElecStatePW", "psiToRho");
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
    // if (PARAM.inp.nspin == 4)
    //     wfcr_another_spin.resize(this->charge->nrxx);

    this->init_rho_data();
    int ik = psi.get_current_k();
    int npw = psi.get_current_ngk();
    int current_spin = 0;
    if (PARAM.inp.nspin == 2)
    {
        current_spin = this->klist->isk[ik];
    }
    int nbands = psi.get_nbands();
    //  here we compute the band energy: the sum of the eigenvalues
    if (PARAM.inp.nspin == 4)
    {
        int npwx = npw / 2;
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            /// be care of when smearing_sigma is large, wg would less than 0
            ///

            this->basis->recip_to_real(this->ctx, &psi(ibnd,0), this->wfcr, ik);

            this->basis->recip_to_real(this->ctx, &psi(ibnd,npwx), this->wfcr_another_spin, ik);

            const auto w1 = static_cast<Real>(this->wg(ik, ibnd) / ucell->omega);

            if (w1 != 0.0)
            {
                // replaced by denghui at 20221110
                elecstate_pw_op()(this->ctx,
                                  PARAM.globalv.domag,
                                  PARAM.globalv.domag_z,
                                  this->basis->nrxx,
                                  w1,
                                  this->rho,
                                  this->wfcr,
                                  this->wfcr_another_spin);
            }
        }
    }
    else
    {
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///

            this->basis->recip_to_real(this->ctx, &psi(ibnd,0), this->wfcr, ik);

            const auto w1 = static_cast<Real>(this->wg(ik, ibnd) / ucell->omega);

            if (w1 != 0.0)
            {
                // replaced by denghui at 20221110
                elecstate_pw_op()(this->ctx, current_spin, this->basis->nrxx, w1, this->rho, this->wfcr);
            }

            // kinetic energy density
            if (XC_Functional::get_func_type() == 3)
            {
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
    }
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::cal_becsum(const psi::Psi<T, Device>& psi)
{
    const T one{1, 0};
    const T zero{0, 0};
    const int npol = psi.get_npol();
    const int npwx = psi.get_nbasis() / npol;
    const int nbands = psi.get_nbands() * npol;
    const int nkb = this->ppcell->nkb;
    this->vkb = this->ppcell->template get_vkb_data<Real>();
    const int nh_tot = this->ppcell->nhm * (this->ppcell->nhm + 1) / 2;
    // becsum on CPU (forces_us / stress_us use CPU dgemm)
    resmem_var_h_op()(becsum, nh_tot * ucell->nat * PARAM.inp.nspin, "ElecState<PW>::becsum");
    setmem_var_h_op()(becsum, 0, nh_tot * ucell->nat * PARAM.inp.nspin);

    // becp: device buffer for gemm, then D2H for host loops
    T* becp = nullptr;
    std::vector<T> becp_host;
    if (nkb > 0)
    {
        resmem_complex_op()(becp, nbands * nkb, "ElecState<PW>::becp");
        becp_host.resize(nbands * nkb);
    }

    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        const T* psi_now = psi.get_pointer();
        const int currect_spin = this->klist->isk[ik];
        const int npw = psi.get_current_ngk();

        // get |beta>
        if (this->ppcell->nkb > 0)
        {
            this->ppcell->getvnl(this->ctx, *ucell,ik, this->vkb);
        }

        // becp = <beta|psi> (device gemm)
        char transa = 'C';
        char transb = 'N';
        if (this->ppcell->nkb > 0)
        {
            if (nbands == 1)
            {
                int inc = 1;
                gemv_op()(transa,
                          npw,
                          this->ppcell->nkb,
                          &one,
                          this->vkb,
                          this->ppcell->vkbnc,
                          psi_now,
                          inc,
                          &zero,
                          becp,
                          inc);
            }
            else
            {
                gemm_op()(transa,
                          transb,
                          this->ppcell->nkb,
                          nbands,
                          npw,
                          &one,
                          this->vkb,
                          this->ppcell->vkbnc,
                          psi_now,
                          npwx,
                          &zero,
                          becp,
                          this->ppcell->nkb);
            }
            // D2H: device becp → host becp_host
            syncmem_complex_d2h_op()(becp_host.data(), becp, nbands * nkb);
        }
        Parallel_Reduce::reduce_pool(becp_host.data(), this->ppcell->nkb * nbands);

        // sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
        for (int it = 0; it < ucell->ntype; it++)
        {
            Atom* atom = &ucell->atoms[it];
            if (atom->ncpp.tvanp)
            {
                const int nh_atom = atom->ncpp.nh;
                // auxk1, auxk2: host buffers for filling loops
                std::vector<T> auxk1_host(nbands * nh_atom);
                std::vector<T> auxk2_host(nbands * nh_atom);
                // device buffers allocated once per atom type
                T *aux_gk = nullptr, *auxk1 = nullptr, *auxk2 = nullptr;
                resmem_complex_op()(auxk1, nbands * nh_atom, "ElecState<PW>::auxk1");
                resmem_complex_op()(auxk2, nbands * nh_atom, "ElecState<PW>::auxk2");
                resmem_complex_op()(aux_gk, nh_atom * nh_atom * npol * npol, "ElecState<PW>::aux_gk");
                for (int ia = 0; ia < atom->na; ia++)
                {
                    const int iat = ucell->itia2iat(it, ia);
                    if (PARAM.inp.noncolin)
                    {
                        // noncolinear case
                    }
                    else
                    {
                        for (int ih = 0; ih < nh_atom; ih++)
                        {
                            const int ikb = this->ppcell->indv_ijkb0[iat] + ih;
                            for (int ib = 0; ib < nbands; ib++)
                            {
                                auxk1_host[ih * nbands + ib] = becp_host[ib * this->ppcell->nkb + ikb];
                                auxk2_host[ih * nbands + ib]
                                    = becp_host[ib * this->ppcell->nkb + ikb] * static_cast<Real>(this->wg(ik, ib));
                            }
                        }

                        // device gemm: aux_gk = auxk1^H * auxk2
                        syncmem_complex_h2d_op()(auxk1, auxk1_host.data(), nbands * nh_atom);
                        syncmem_complex_h2d_op()(auxk2, auxk2_host.data(), nbands * nh_atom);

                        char transa2 = 'C';
                        char transb2 = 'N';
                        gemm_op()(transa2,
                                  transb2,
                                  nh_atom,
                                  nh_atom,
                                  nbands,
                                  &one,
                                  auxk1,
                                  nbands,
                                  auxk2,
                                  nbands,
                                  &zero,
                                  aux_gk,
                                  nh_atom);

                        // D2H: device aux_gk → host
                        std::vector<T> aux_gk_host(nh_atom * nh_atom);
                        syncmem_complex_d2h_op()(aux_gk_host.data(), aux_gk, nh_atom * nh_atom);

                        // copy output from GEMM into desired format
                        int ijh = 0;
                        const int index = currect_spin * ucell->nat * nh_tot + iat * nh_tot;
                        for (int ih = 0; ih < nh_atom; ih++)
                        {
                            for (int jh = ih; jh < nh_atom; jh++)
                            {
                                if (ih == jh)
                                {
                                    becsum[index + ijh] += std::real(aux_gk_host[ih * nh_atom + jh]);
                                }
                                else
                                {
                                    becsum[index + ijh] += 2.0 * std::real(aux_gk_host[ih * nh_atom + jh]);
                                }
                                ijh++;
                            }
                        }
                    }
                }
                delmem_complex_op()(auxk1);
                delmem_complex_op()(auxk2);
                delmem_complex_op()(aux_gk);
            }
        }
    }
    delmem_complex_op()(becp);
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::add_usrho(const psi::Psi<T, Device>& psi)
{
    if (PARAM.globalv.use_uspp)
    {
        this->cal_becsum(psi);
    }

    // transform soft charge to recip space using smooth grids
    if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            // D2H rho → CPU buffer (real2recip reads rho via CPU loops)
            std::vector<Real> rho_host(this->rhopw_smooth->nrxx);
            syncmem_var_d2h_op()(rho_host.data(), this->rho[is], this->rhopw_smooth->nrxx);
            // D2H rhog → CPU temp
            std::vector<T> rhog_host(this->rhopw_smooth->npw);
            syncmem_complex_d2h_op()(rhog_host.data(), this->rhog[is], this->rhopw_smooth->npw);
            // CPU real2recip
            this->rhopw_smooth->real2recip(rho_host.data(), rhog_host.data());
            // H2D rhog back to device
            syncmem_complex_h2d_op()(this->rhog[is], rhog_host.data(), this->rhopw_smooth->npw);
        }
    }

    // \sum_lm Q_lm(r) \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
    // add to the charge density in reciprocal space the part which is due to the US augmentation.
    if (PARAM.globalv.use_uspp)
    {
        this->addusdens_g(becsum, rhog);
    }
    // transform back to real space using dense grids
    if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            // D2H rhog → CPU temp (recip2real reads rhog via CPU loops)
            std::vector<T> rhog_host(this->charge->rhopw->npw);
            syncmem_complex_d2h_op()(rhog_host.data(), this->rhog[is], this->charge->rhopw->npw);
            // D2H rho → CPU buffer
            std::vector<Real> rho_host(this->charge->rhopw->nrxx);
            syncmem_var_d2h_op()(rho_host.data(), this->rho[is], this->charge->rhopw->nrxx);
            // CPU recip2real
            this->charge->rhopw->recip2real(rhog_host.data(), rho_host.data());
            // H2D rho back to GPU
            syncmem_var_h2d_op()(this->rho[is], rho_host.data(), this->charge->rhopw->nrxx);
        }
    }
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::addusdens_g(const Real* becsum, T** rhog)
{
    const T one{1, 0};
    const T zero{0, 0};
    const int npw = this->charge->rhopw->npw;
    const int lmaxq = this->ppcell->lmaxq;
    const int nh_tot = this->ppcell->nhm * (this->ppcell->nhm + 1) / 2;
    Structure_Factor* psf = this->ppcell->psf;
    const std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;

    // ---------- all on CPU ----------
    std::vector<double> qmod_host(npw);
    std::vector<std::complex<double>> qgm_host(npw);
    for (int ig = 0; ig < npw; ig++)
    {
        qmod_host[ig] = static_cast<double>(this->charge->rhopw->gcar[ig].norm() * ucell->tpiba);
    }

    // ylmk0: compute on device then D2H
    Real* ylmk0 = nullptr;
    resmem_var_op()(ylmk0, npw * lmaxq * lmaxq, "ElecState<PW>::ylmk0");
    Real* g = nullptr;
    resmem_var_op()(g, npw * 3, "ElecState<PW>::g");
    syncmem_var_h2d_op()(g, reinterpret_cast<Real*>(this->charge->rhopw->gcar), npw * 3);
    ModuleBase::YlmReal::Ylm_Real(this->ctx, lmaxq * lmaxq, npw, g, ylmk0);
    delmem_var_op()(g);
    std::vector<Real> ylmk0_host(npw * lmaxq * lmaxq);
    syncmem_var_d2h_op()(ylmk0_host.data(), ylmk0, npw * lmaxq * lmaxq);
    std::vector<double> ylmk0_double(npw * lmaxq * lmaxq);
    for (int i = 0; i < npw * lmaxq * lmaxq; i++) ylmk0_double[i] = static_cast<double>(ylmk0_host[i]);

    for (int it = 0; it < ucell->ntype; it++)
    {
        Atom* atom = &ucell->atoms[it];
        if (atom->ncpp.tvanp)
        {
            const int nij = atom->ncpp.nh * (atom->ncpp.nh + 1) / 2;

            // skk, tbecsum: CPU vectors
            std::vector<std::complex<double>> skk_host(atom->na * npw);
            std::vector<std::complex<double>> tbecsum_host(PARAM.inp.nspin * atom->na * nij);
            for (int ia = 0; ia < atom->na; ia++)
            {
                const int iat = ucell->itia2iat(it, ia);
                for (int is = 0; is < PARAM.inp.nspin; is++)
                {
                    for (int ij = 0; ij < nij; ij++)
                    {
                        tbecsum_host[is * atom->na * nij + ia * nij + ij]
                            = static_cast<std::complex<double>>(becsum[is * ucell->nat * nh_tot + iat * nh_tot + ij]);
                    }
                }
                for (int ig = 0; ig < npw; ig++)
                {
                    double arg = this->charge->rhopw->gcar[ig] * atom->tau[ia];
                    skk_host[ia * npw + ig] = ModuleBase::libm::exp(ci_tpi * arg);
                }
            }

            for (int is = 0; is < PARAM.inp.nspin; is++)
            {
                // CPU BLAS: aux2 = skk * tbecsum^T
                std::vector<std::complex<double>> aux2_host(nij * npw);
                const std::complex<double> one_d(1, 0), zero_d(0, 0);
                char transa = 'N';
                char transb = 'T';
                zgemm_(&transa, &transb, &npw, &nij, &atom->na,
                       &one_d, skk_host.data(), &npw,
                       &tbecsum_host[is * atom->na * nij], &nij,
                       &zero_d, aux2_host.data(), &npw);

                // D2H rhog
                std::vector<T> rhog_host(npw);
                syncmem_complex_d2h_op()(rhog_host.data(), rhog[is], npw);

                int ijh = 0;
                for (int ih = 0; ih < atom->ncpp.nh; ih++)
                {
                    for (int jh = ih; jh < atom->ncpp.nh; jh++)
                    {
                        // CPU radial_fft_q (template version, DEVICE_CPU)
                        this->ppcell->template radial_fft_q<double, base_device::DEVICE_CPU>(
                            nullptr, npw, ih, jh, it, qmod_host.data(), ylmk0_double.data(), qgm_host.data());
                        for (int ig = 0; ig < npw; ig++)
                        {
                            rhog_host[ig] += static_cast<T>(qgm_host[ig] * aux2_host[ijh * npw + ig]);
                        }
                        ijh++;
                    }
                }
                // H2D rhog back
                syncmem_complex_h2d_op()(rhog[is], rhog_host.data(), npw);
            }
        }
    }

    delmem_var_op()(ylmk0);
}

template class ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ElecStatePW<std::complex<float>, base_device::DEVICE_GPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>;
#endif 

} // namespace elecstate
