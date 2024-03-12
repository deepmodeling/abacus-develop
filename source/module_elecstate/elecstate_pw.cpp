#include "elecstate_pw.h"

#include "elecstate_getters.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/math_ylmreal.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_psi/kernels/device.h"

namespace elecstate {

template <typename T, typename Device>
ElecStatePW<T, Device>::ElecStatePW(ModulePW::PW_Basis_K* wfc_basis_in,
                                    Charge* chg_in,
                                    K_Vectors* pkv_in,
                                    UnitCell* ucell_in,
                                    pseudopot_cell_vnl* ppcell_in,
                                    ModulePW::PW_Basis* rhodpw_in,
                                    ModulePW::PW_Basis* rhopw_in,
                                    ModulePW::PW_Basis_Big* bigpw_in)
    : basis(wfc_basis_in)
{
    this->classname = "ElecStatePW";
    this->rhopw_smooth = rhopw_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->init_ks(chg_in, pkv_in, pkv_in->nks, rhodpw_in, bigpw_in);
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
    delmem_var_op()(this->ctx, becsum);
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
    if (GlobalV::use_uspp)
    {
        this->add_usrho(psi);
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

            this->basis->recip_to_real(this->ctx, &psi(ibnd,0), this->wfcr, ik);

            this->basis->recip_to_real(this->ctx, &psi(ibnd,npwx), this->wfcr_another_spin, ik);

            const auto w1 = static_cast<Real>(this->wg(ik, ibnd) / get_ucell_omega());

            if (w1 != 0.0)
            {
                // replaced by denghui at 20221110
                elecstate_pw_op()(this->ctx, GlobalV::DOMAG, GlobalV::DOMAG_Z, this->charge->nrxx, w1, this->rho, this->wfcr, this->wfcr_another_spin);
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

template <typename T, typename Device>
void ElecStatePW<T, Device>::add_usrho(const psi::Psi<T, Device>& psi)
{
    const T one{1, 0};
    const T zero{0, 0};
    const int npol = psi.npol;
    const int npwx = psi.get_nbasis() / npol;
    const int nbands = psi.get_nbands() * npol;
    const int nkb = this->ppcell->nkb;
    this->vkb = this->ppcell->template get_vkb_data<Real>();
    T* becp = nullptr;
    resmem_complex_op()(this->ctx, becp, nbands * nkb, "ElecState<PW>::becp");
    const int nh_tot = this->ppcell->nhm * (this->ppcell->nhm + 1) / 2;
    resmem_var_op()(this->ctx, becsum, nh_tot * ucell->nat * GlobalV::NSPIN, "ElecState<PW>::becsum");
    setmem_var_op()(this->ctx, becsum, 0, nh_tot * ucell->nat * GlobalV::NSPIN);

    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        const T* psi_now = psi.get_pointer();
        const int currect_spin = this->klist->isk[ik];
        const int npw = psi.get_current_nbas();

        // get |beta>
        if (this->ppcell->nkb > 0)
        {
            this->ppcell->getvnl(this->ctx, ik, this->vkb);
        }

        // becp = <beta|psi>
        char transa = 'C';
        char transb = 'N';
        if (nbands == 1)
        {
            int inc = 1;
            gemv_op()(this->ctx,
                      transa,
                      npw,
                      this->ppcell->nkb,
                      &one,
                      this->vkb,
                      this->ppcell->vkb.nc,
                      psi_now,
                      inc,
                      &zero,
                      becp,
                      inc);
        }
        else
        {
            gemm_op()(this->ctx,
                      transa,
                      transb,
                      this->ppcell->nkb,
                      nbands,
                      npw,
                      &one,
                      this->vkb,
                      this->ppcell->vkb.nc,
                      psi_now,
                      npwx,
                      &zero,
                      becp,
                      this->ppcell->nkb);
        }
        Parallel_Reduce::reduce_pool(becp, this->ppcell->nkb * nbands);

        // sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
        for (int it = 0; it < ucell->ntype; it++)
        {
            Atom* atom = &ucell->atoms[it];
            if (atom->ncpp.tvanp)
            {
                T *auxk1 = nullptr, *auxk2 = nullptr, *aux_gk = nullptr;
                resmem_complex_op()(this->ctx, auxk1, nbands * atom->ncpp.nh, "ElecState<PW>::auxk1");
                resmem_complex_op()(this->ctx, auxk2, nbands * atom->ncpp.nh, "ElecState<PW>::auxk2");
                resmem_complex_op()(this->ctx,
                                    aux_gk,
                                    atom->ncpp.nh * atom->ncpp.nh * npol * npol,
                                    "ElecState<PW>::aux_gk");
                for (int ia = 0; ia < atom->na; ia++)
                {
                    const int iat = ucell->itia2iat(it, ia);
                    if (GlobalV::NONCOLIN)
                    {
                        // noncolinear case
                    }
                    else
                    {
                        for (int ih = 0; ih < atom->ncpp.nh; ih++)
                        {
                            const int ikb = this->ppcell->indv_ijkb0[iat] + ih;
                            for (int ib = 0; ib < nbands; ib++)
                            {
                                auxk1[ih * nbands + ib] = becp[ib * this->ppcell->nkb + ikb];
                                auxk2[ih * nbands + ib]
                                    = becp[ib * this->ppcell->nkb + ikb] * static_cast<Real>(this->wg(ik, ib));
                            }
                        }

                        char transa = 'C';
                        char transb = 'N';
                        gemm_op()(this->ctx,
                                  transa,
                                  transb,
                                  atom->ncpp.nh,
                                  atom->ncpp.nh,
                                  nbands,
                                  &one,
                                  auxk1,
                                  nbands,
                                  auxk2,
                                  nbands,
                                  &zero,
                                  aux_gk,
                                  atom->ncpp.nh);
                    }

                    // copy output from GEMM into desired format
                    if (GlobalV::NONCOLIN && !atom->ncpp.has_so)
                    {
                    }
                    else if (GlobalV::NONCOLIN && atom->ncpp.has_so)
                    {
                    }
                    else
                    {
                        int ijh = 0;
                        const int index = currect_spin * ucell->nat * nh_tot + iat * nh_tot;
                        for (int ih = 0; ih < atom->ncpp.nh; ih++)
                        {
                            for (int jh = ih; jh < atom->ncpp.nh; jh++)
                            {
                                // nondiagonal terms summed and collapsed into a
                                // single index (matrix is symmetric wrt (ih,jh))
                                if (ih == jh)
                                {
                                    becsum[index + ijh] += std::real(aux_gk[ih * atom->ncpp.nh + jh]);
                                }
                                else
                                {
                                    becsum[index + ijh] += 2.0 * std::real(aux_gk[ih * atom->ncpp.nh + jh]);
                                }
                                ijh++;
                            }
                        }
                    }
                }
                delmem_complex_op()(this->ctx, auxk1);
                delmem_complex_op()(this->ctx, auxk2);
                delmem_complex_op()(this->ctx, aux_gk);
            }
        }
    }
    delmem_complex_op()(this->ctx, becp);

    // transform soft charge to recip space using smooth grids
    T* rhog = nullptr;
    resmem_complex_op()(this->ctx, rhog, this->charge->rhopw->npw * GlobalV::NSPIN, "ElecState<PW>::rhog");
    setmem_complex_op()(this->ctx, rhog, 0, this->charge->rhopw->npw * GlobalV::NSPIN);
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        this->rhopw_smooth->real2recip(this->rho[is], &rhog[is * this->charge->rhopw->npw]);
    }

    // \sum_lm Q_lm(r) \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
    // add to the charge density in reciprocal space the part which is due to the US augmentation.
    this->addusdens_g(becsum, rhog);

    // transform back to real space using dense grids
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        this->charge->rhopw->recip2real(&rhog[is * this->charge->rhopw->npw], this->rho[is]);
    }

    delmem_complex_op()(this->ctx, rhog);
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::addusdens_g(const Real* becsum, T* rhog)
{
    const T one{1, 0};
    const T zero{0, 0};
    const int npw = this->charge->rhopw->npw;
    const int lmaxq = this->ppcell->lmaxq;
    const int nh_tot = this->ppcell->nhm * (this->ppcell->nhm + 1) / 2;
    Structure_Factor* psf = this->ppcell->psf;
    const std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;

    Real* qmod = nullptr;
    resmem_var_op()(this->ctx, qmod, npw, "ElecState<PW>::qmod");
    T* qgm = nullptr;
    resmem_complex_op()(this->ctx, qgm, npw, "ElecState<PW>::qgm");
    Real* ylmk0 = nullptr;
    resmem_var_op()(this->ctx, ylmk0, npw * lmaxq * lmaxq, "ElecState<PW>::ylmk0");
    Real* g = reinterpret_cast<Real*>(this->charge->rhopw->gcar);

    ModuleBase::YlmReal::Ylm_Real(this->ctx, lmaxq * lmaxq, npw, g, ylmk0);

    for (int ig = 0; ig < npw; ig++)
    {
        qmod[ig] = static_cast<Real>(this->charge->rhopw->gcar[ig].norm() * ucell->tpiba);
    }

    for (int it = 0; it < ucell->ntype; it++)
    {
        Atom* atom = &ucell->atoms[it];
        if (atom->ncpp.tvanp)
        {
            // nij = max number of (ih,jh) pairs per atom type nt
            const int nij = atom->ncpp.nh * (atom->ncpp.nh + 1) / 2;

            T *skk = nullptr, *aux2 = nullptr, *tbecsum = nullptr;
            resmem_complex_op()(this->ctx, skk, atom->na * npw, "ElecState<PW>::skk");
            resmem_complex_op()(this->ctx, aux2, nij * npw, "ElecState<PW>::aux2");
            resmem_complex_op()(this->ctx, tbecsum, GlobalV::NSPIN * atom->na * nij, "ElecState<PW>::tbecsum");
            for (int ia = 0; ia < atom->na; ia++)
            {
                const int iat = ucell->itia2iat(it, ia);
                for (int is = 0; is < GlobalV::NSPIN; is++)
                {
                    for (int ij = 0; ij < nij; ij++)
                    {
                        tbecsum[is * atom->na * nij + ia * nij + ij]
                            = static_cast<T>(becsum[is * ucell->nat * nh_tot + iat * nh_tot + ij]);
                    }
                }
                for (int ig = 0; ig < npw; ig++)
                {
                    double arg = this->charge->rhopw->gcar[ig] * atom->tau[ia];
                    skk[ia * npw + ig] = static_cast<T>(ModuleBase::libm::exp(ci_tpi * arg));
                }
            }

            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                // sum over atoms
                char transa = 'N';
                char transb = 'T';
                gemm_op()(this->ctx,
                          transa,
                          transb,
                          npw,
                          nij,
                          atom->na,
                          &one,
                          skk,
                          npw,
                          &tbecsum[is * atom->na * nij],
                          nij,
                          &zero,
                          aux2,
                          npw);

                // sum over lm indices of Q_{lm}
                int ijh = 0;
                for (int ih = 0; ih < atom->ncpp.nh; ih++)
                {
                    for (int jh = ih; jh < atom->ncpp.nh; jh++)
                    {
                        this->ppcell->radial_fft_q(this->ctx, npw, ih, jh, it, qmod, ylmk0, qgm);
                        for (int ig = 0; ig < npw; ig++)
                        {
                            rhog[is * npw + ig] += qgm[ig] * aux2[ijh * npw + ig];
                        }
                        ijh++;
                    }
                }
            }
            delmem_complex_op()(this->ctx, skk);
            delmem_complex_op()(this->ctx, aux2);
            delmem_complex_op()(this->ctx, tbecsum);
        }
    }

    delmem_var_op()(this->ctx, qmod);
    delmem_complex_op()(this->ctx, qgm);
    delmem_var_op()(this->ctx, ylmk0);
}

template class ElecStatePW<std::complex<float>, psi::DEVICE_CPU>;
template class ElecStatePW<std::complex<double>, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ElecStatePW<std::complex<float>, psi::DEVICE_GPU>;
template class ElecStatePW<std::complex<double>, psi::DEVICE_GPU>;
#endif 

} // namespace elecstate
