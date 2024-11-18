#include "sto_iter.h"

#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_elecstate/kernels/elecstate_op.h"

template <typename T, typename Device>
Stochastic_Iter<T, Device>::Stochastic_Iter()
{
    change = false;
    mu0 = 0;
    method = 2;
}

template <typename T, typename Device>
Stochastic_Iter<T, Device>::~Stochastic_Iter()
{
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::dot(const int& n, const Real* x, const int& incx, const Real* y, const int& incy, Real& result)
{
    Real* result_device = nullptr;
    resmem_var_op()(this->ctx, result_device, 1);
    container::kernels::blas_dot<Real, ct_Device>()(n, p_che->coef_real, 1, spolyv, 1, result_device);
    syncmem_var_d2h_op()(cpu_ctx, this->ctx, &result, result_device, 1);
    delmem_var_op()(this->ctx, result_device);
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::init(K_Vectors* pkv_in,
                                      ModulePW::PW_Basis_K* wfc_basis,
                                      Stochastic_WF<T, Device>& stowf,
                                      StoChe<Real, Device>& stoche,
                                      hamilt::HamiltSdftPW<T, Device>* p_hamilt_sto)
{
    p_che = stoche.p_che;
    spolyv = stoche.spolyv;
    spolyv_cpu = stoche.spolyv_cpu;
    nchip = stowf.nchip;
    targetne = PARAM.inp.nelec;
    this->pkv = pkv_in;
    this->method = stoche.method_sto;
    this->p_hamilt_sto = p_hamilt_sto;
    this->stofunc.set_E_range(&stoche.emin_sto, &stoche.emax_sto);
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::orthog(const int& ik, psi::Psi<T, Device>& psi, Stochastic_WF<T, Device>& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter", "orthog");
    ModuleBase::timer::tick("Stochastic_Iter", "orthog");
    // orthogonal part
    if (PARAM.inp.nbands > 0)
    {
        const int nchipk = stowf.nchip[ik];
        const int npw = psi.get_current_nbas();
        const int npwx = psi.get_nbasis();
        stowf.chi0->fix_k(ik);
        stowf.chiortho->fix_k(ik);
        T *wfgin = stowf.chi0->get_pointer(), *wfgout = stowf.chiortho->get_pointer();
        cpymem_complex_op()(this->ctx, this->ctx, wfgout, wfgin, npwx * nchipk);
        // for (int ig = 0; ig < npwx * nchipk; ++ig)
        // {
        //     wfgout[ig] = wfgin[ig];
        // }

        // orthogonal part
        T* sum = nullptr;
        resmem_complex_op()(this->ctx, sum, PARAM.inp.nbands * nchipk);
        char transC = 'C';
        char transN = 'N';

        // sum(b<NBANDS, a<nchi) = < psi_b | chi_a >
        hsolver::gemm_op<T, Device>()(ctx,
                                      transC,
                                      transN,
                                      PARAM.inp.nbands,
                                      nchipk,
                                      npw,
                                      &ModuleBase::ONE,
                                      &psi(ik, 0, 0),
                                      npwx,
                                      wfgout,
                                      npwx,
                                      &ModuleBase::ZERO,
                                      sum,
                                      PARAM.inp.nbands);
        Parallel_Reduce::reduce_pool(sum, PARAM.inp.nbands * nchipk);

        // psi -= psi * sum
        hsolver::gemm_op<T, Device>()(ctx,
                                      transN,
                                      transN,
                                      npw,
                                      nchipk,
                                      PARAM.inp.nbands,
                                      &ModuleBase::NEG_ONE,
                                      &psi(ik, 0, 0),
                                      npwx,
                                      sum,
                                      PARAM.inp.nbands,
                                      &ModuleBase::ONE,
                                      wfgout,
                                      npwx);
        delmem_complex_op()(this->ctx, sum);
    }
    ModuleBase::timer::tick("Stochastic_Iter", "orthog");
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::checkemm(const int& ik,
                                          const int istep,
                                          const int iter,
                                          Stochastic_WF<T, Device>& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter", "checkemm");
    ModuleBase::timer::tick("Stochastic_Iter", "checkemm");
    // iter = 1,2,...   istep = 0,1,2,...
    //  if( istep%PARAM.inp.initsto_freq != 0 )    return;
    const int npw = stowf.ngk[ik];
    const int nks = stowf.nks;
    if (istep == 0)
    {
        if (iter > 5)
        {
            return;
        }
    }
    else
    {
        if (iter > 1)
        {
            return;
        }
    }

    const int norder = p_che->norder;
    T* pchi;
    int ntest = 1;

    if (nchip[ik] < ntest)
    {
        ntest = nchip[ik];
    }

    for (int ichi = 0; ichi < ntest; ++ichi)
    {
        if (PARAM.inp.nbands > 0)
        {
            pchi = &stowf.chiortho->operator()(ik, ichi, 0);
        }
        else
        {
            pchi = &stowf.chi0->operator()(ik, ichi, 0);
        }
        while (true)
        {
            bool converge;
            auto hchi_norm = std::bind(&hamilt::HamiltSdftPW<T, Device>::hPsi_norm,
                                       p_hamilt_sto,
                                       std::placeholders::_1,
                                       std::placeholders::_2,
                                       std::placeholders::_3);
            converge
                = p_che->checkconverge(hchi_norm, pchi, npw, stowf.npwx, *p_hamilt_sto->emax, *p_hamilt_sto->emin, 5.0);

            if (!converge)
            {
                change = true;
            }
            else
            {
                break;
            }
        }
    }
    if (ik == nks - 1)
    {
#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, p_hamilt_sto->emax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, p_hamilt_sto->emin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &change, 1, MPI_CHAR, MPI_LOR, MPI_COMM_WORLD);
#endif
        if (change)
        {
            GlobalV::ofs_running << "New Emax Ry" << *p_hamilt_sto->emax << " ; new Emin " << *p_hamilt_sto->emin
                                 << " Ry" << std::endl;
        }
        change = false;
    }
    ModuleBase::timer::tick("Stochastic_Iter", "checkemm");
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::check_precision(const double ref, const double thr, const std::string info)
{
    //==============================
    // precision check
    //==============================
    double error = 0;
    if (this->method == 1)
    {
        Real last_coef = 0;
        Real last_spolyv = 0;
        syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, &last_coef, &p_che->coef_real[p_che->norder - 1], 1);
        syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, &last_spolyv, &spolyv[p_che->norder - 1], 1);
        error = last_coef * last_spolyv;
    }
    else
    {
        const int norder = p_che->norder;
        // double last_coef = p_che->coef_real[norder - 1];
        // double last_spolyv = spolyv[norder * norder - 1];
        Real last_coef = 0;
        Real last_spolyv = 0;
        syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, &last_coef, &p_che->coef_real[norder - 1], 1);
        syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, &last_spolyv, &spolyv[norder * norder - 1], 1);
        Real dot1 = 0, dot2 = 0;
        this->dot(norder, p_che->coef_real, 1, spolyv + norder * (norder - 1), 1, dot1);
        this->dot(norder, p_che->coef_real, 1, spolyv + norder - 1, norder, dot2);
        error = last_coef * (dot1 + dot2 - last_coef * last_spolyv);
    }

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    double relative_error = std::abs(error / ref);
    GlobalV::ofs_running << info << "Relative Chebyshev Precision: " << relative_error * 1e9 << "E-09" << std::endl;
    if (relative_error > thr)
    {
        std::stringstream ss;
        ss << relative_error;
        std::string fractxt, tartxt;
        ss >> fractxt;
        ss.clear();
        ss << thr;
        ss >> tartxt;
        std::string warningtxt = "( " + info + " relative Chebyshev error = " + fractxt + " > threshold = " + tartxt
                                 + " ) Maybe you should increase the parameter \"nche_sto\" for more accuracy.";
        ModuleBase::WARNING("Stochastic_Chebychev", warningtxt);
    }
    //===============================
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::itermu(const int iter, elecstate::ElecState* pes)
{
    ModuleBase::TITLE("Stochastic_Iter", "itermu");
    ModuleBase::timer::tick("Stochastic_Iter", "itermu");
    double dmu;
    if (iter == 1)
    {
        dmu = 2;
        th_ne = 0.1 * PARAM.inp.scf_thr * PARAM.inp.nelec;
        // std::cout<<"th_ne "<<th_ne<<std::endl;
    }
    else
    {
        dmu = 0.1;
        th_ne = 1e-2 * PARAM.inp.scf_thr * PARAM.inp.nelec;
        th_ne = std::min(th_ne, 1e-5);
    }
    this->stofunc.mu = mu0 - dmu;
    double ne1 = calne(pes);
    double mu1 = this->stofunc.mu;

    this->stofunc.mu = mu0 + dmu;
    double ne2 = calne(pes);
    double mu2 = this->stofunc.mu;
    double Dne = th_ne + 1;
    double ne3;
    double mu3;

    while (ne1 > targetne)
    {
        mu2 = mu1;
        mu1 -= dmu;
        this->stofunc.mu = mu1;
        ne1 = calne(pes);
        // std::cout<<"Reset mu1 from "<<mu1+dmu<<" to "<<mu1<<std::endl;
        dmu *= 2;
    }
    while (ne2 < targetne)
    {
        mu1 = mu2;
        mu2 += dmu;
        this->stofunc.mu = mu2;
        ne2 = calne(pes);
        // cout<<"Reset mu2 from "<<mu2-dmu<<" to "<<mu2<<endl;
        dmu *= 2;
    }
    int count = 0;
    while (Dne > th_ne)
    {
        mu3 = (mu2 + mu1) / 2;
        this->stofunc.mu = mu3;
        ne3 = calne(pes);
        if (ne3 < targetne)
        {
            mu1 = mu3;
        }
        else if (ne3 > targetne)
        {
            mu2 = mu3;
        }
        Dne = std::abs(targetne - ne3);

        count++;
        if (count > 60)
        {
            std::cout << "Fermi energy cannot be converged. Set THNE to " << th_ne << std::endl;
            th_ne *= 1e1;
            if (th_ne > 1e1)
            {
                ModuleBase::WARNING_QUIT("Stochastic_Iter",
                                         "Cannot converge feimi energy. Please retry with different random number");
            }
        }
    }
    pes->eferm.ef = this->stofunc.mu = mu0 = mu3;
    GlobalV::ofs_running << "Converge fermi energy = " << mu3 << " Ry in " << count << " steps." << std::endl;
    this->check_precision(targetne, 10 * PARAM.inp.scf_thr, "Ne");

    // Set wf.wg
    if (PARAM.inp.nbands > 0)
    {
        for (int ikk = 0; ikk < this->pkv->get_nks(); ++ikk)
        {
            double* en = &pes->ekb(ikk, 0);
            for (int iksb = 0; iksb < PARAM.inp.nbands; ++iksb)
            {
                pes->wg(ikk, iksb) = stofunc.fd(en[iksb]) * this->pkv->wk[ikk];
            }
        }
    }
    ModuleBase::timer::tick("Stochastic_Iter", "itermu");
    return;
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::calPn(const int& ik, Stochastic_WF<T, Device>& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter", "calPn");
    ModuleBase::timer::tick("Stochastic_Iter", "calPn");

    const int norder = p_che->norder;
    const int nchip_ik = nchip[ik];
    const int npw = stowf.ngk[ik];
    const int npwx = stowf.npwx;
    if (ik == 0)
    {
        if (this->method == 1)
        {
            ModuleBase::GlobalFunc::ZEROS(spolyv_cpu, norder);
        }
        else
        {
            setmem_var_op()(this->ctx, spolyv, 0, norder * norder);
        }
    }
    T* pchi;
    if (PARAM.inp.nbands > 0)
    {
        stowf.chiortho->fix_k(ik);
        pchi = stowf.chiortho->get_pointer();
    }
    else
    {
        stowf.chi0->fix_k(ik);
        pchi = stowf.chi0->get_pointer();
    }

    auto hchi_norm = std::bind(&hamilt::HamiltSdftPW<T, Device>::hPsi_norm,
                               p_hamilt_sto,
                               std::placeholders::_1,
                               std::placeholders::_2,
                               std::placeholders::_3);
    if (this->method == 1)
    {
        p_che->tracepolyA(hchi_norm, pchi, npw, npwx, nchip_ik);
        for (int i = 0; i < norder; ++i)
        {
            spolyv_cpu[i] += p_che->polytrace[i] * this->pkv->wk[ik];
        }
        if(ik == this->pkv->get_nks() - 1)
        {
            syncmem_var_h2d_op()(this->ctx, cpu_ctx, spolyv, spolyv_cpu, norder);
        }
    }
    else
    {
        p_che->calpolyvec_complex(hchi_norm, pchi, stowf.chiallorder[ik].get_pointer(), npw, npwx, nchip_ik);
        const Real* vec_all = (Real*)stowf.chiallorder[ik].get_pointer();
        const char trans = 'T';
        const char normal = 'N';
        const Real one = 1;
        const int LDA = npwx * nchip_ik * 2;
        const int M = npwx * nchip_ik * 2; // Do not use kv.ngk[ik]
        const int N = norder;
        const Real kweight = this->pkv->wk[ik];
        
        hsolver::gemm_op<Real, Device>()(this->ctx, trans, normal, N, N, M, &kweight, vec_all, LDA, vec_all, LDA, &one, spolyv, N);
        // dgemm_(&trans, &normal, &N, &N, &M, &kweight, vec_all, &LDA, vec_all, &LDA, &one, spolyv, &N);
    }
    ModuleBase::timer::tick("Stochastic_Iter", "calPn");
    return;
}

template <typename T, typename Device>
double Stochastic_Iter<T, Device>::calne(elecstate::ElecState* pes)
{
    ModuleBase::TITLE("Stochastic_Iter", "calne");
    ModuleBase::timer::tick("Stochastic_Iter", "calne");
    double totne = 0;
    KS_ne = 0;
    const int norder = p_che->norder;
    double sto_ne;
    if (this->method == 1)
    {
        // Note: spolyv contains kv.wk[ik]
        auto nfd = std::bind(&Sto_Func<double>::nfd, &this->stofunc, std::placeholders::_1);
        p_che->calcoef_real(nfd);
        this->dot(norder, p_che->coef_real, 1, spolyv, 1, sto_ne);
    }
    else
    {
        auto nroot_fd = std::bind(&Sto_Func<double>::nroot_fd, &this->stofunc, std::placeholders::_1);
        p_che->calcoef_real(nroot_fd);
        sto_ne = vTMv<Real, Device>(p_che->coef_real, spolyv, norder);
    }
    if (PARAM.inp.nbands > 0)
    {
        for (int ikk = 0; ikk < this->pkv->get_nks(); ++ikk)
        {
            double* en = &pes->ekb(ikk, 0);
            for (int iksb = 0; iksb < PARAM.inp.nbands; ++iksb)
            {
                KS_ne += stofunc.fd(en[iksb]) * this->pkv->wk[ikk];
            }
        }
    }
    KS_ne /= GlobalV::NPROC_IN_POOL;
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &KS_ne, 1, MPI_DOUBLE, MPI_SUM, STO_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sto_ne, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    totne = KS_ne + sto_ne;
    ModuleBase::timer::tick("Stochastic_Iter", "calne");
    return totne;
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::calHsqrtchi(Stochastic_WF<T, Device>& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter", "calHsqrtchi");
    ModuleBase::timer::tick("Stochastic_Iter", "calHsqrtchi");
    auto nroot_fd = std::bind(&Sto_Func<double>::nroot_fd, &this->stofunc, std::placeholders::_1);
    p_che->calcoef_real(nroot_fd);
    for (int ik = 0; ik < this->pkv->get_nks(); ++ik)
    {
        this->calTnchi_ik(ik, stowf);
    }
    ModuleBase::timer::tick("Stochastic_Iter", "calHsqrtchi");
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::sum_stoband(Stochastic_WF<T, Device>& stowf,
                                             elecstate::ElecStatePW<T, Device>* pes,
                                             hamilt::Hamilt<T, Device>* pHamilt,
                                             ModulePW::PW_Basis_K* wfc_basis)
{
    ModuleBase::TITLE("Stochastic_Iter", "sum_stoband");
    ModuleBase::timer::tick("Stochastic_Iter", "sum_stoband");
    int nrxx = wfc_basis->nrxx;
    int npwx = wfc_basis->npwk_max;
    const int norder = p_che->norder;

    //---------------cal demet-----------------------
    Real stodemet;
    if (this->method == 1)
    {
        auto nfdlnfd = std::bind(&Sto_Func<double>::nfdlnfd, &this->stofunc, std::placeholders::_1);
        p_che->calcoef_real(nfdlnfd);
        this->dot(norder, p_che->coef_real, 1, spolyv, 1, stodemet);
    }
    else
    {
        auto nroot_fdlnfd = std::bind(&Sto_Func<double>::n_root_fdlnfd, &this->stofunc, std::placeholders::_1);
        p_che->calcoef_real(nroot_fdlnfd);
        stodemet = -vTMv<Real, Device>(p_che->coef_real, spolyv, norder);
    }

    if (PARAM.inp.nbands > 0)
    {
        for (int ikk = 0; ikk < this->pkv->get_nks(); ++ikk)
        {
            double* enb = &pes->ekb(ikk, 0);
            // number of electrons in KS orbitals
            for (int iksb = 0; iksb < PARAM.inp.nbands; ++iksb)
            {
                pes->f_en.demet += stofunc.fdlnfd(enb[iksb]) * this->pkv->wk[ikk];
            }
        }
    }
    pes->f_en.demet /= GlobalV::NPROC_IN_POOL;
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &pes->f_en.demet, 1, MPI_DOUBLE, MPI_SUM, STO_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &stodemet, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    pes->f_en.demet += stodemet;
    this->check_precision(pes->f_en.demet, 1e-4, "TS");
    pes->f_en.demet *= Occupy::gaussian_parameter;

    //--------------------cal eband------------------------
    double sto_eband = 0;
    if (this->method == 1)
    {
        auto nxfd = std::bind(&Sto_Func<double>::nxfd, &this->stofunc, std::placeholders::_1);
        p_che->calcoef_real(nxfd);
        this->dot(norder, p_che->coef_real, 1, spolyv, 1, sto_eband);
    }
    else
    {
        for (int ik = 0; ik < this->pkv->get_nks(); ++ik)
        {
            const int nchip_ik = nchip[ik];
            if (this->pkv->get_nks() > 1)
            {
                pHamilt->updateHk(ik); // can be merged with calTnchi_ik, but it does not nearly cost time.
                stowf.shchi->fix_k(ik);
            }
            const int npw = this->pkv->ngk[ik];
            const double kweight = this->pkv->wk[ik];
            T* hshchi = nullptr;
            resmem_complex_op()(this->ctx, hshchi, nchip_ik * npwx);
            T* tmpin = stowf.shchi->get_pointer();
            T* tmpout = hshchi;
            p_hamilt_sto->hPsi(tmpin, tmpout, nchip_ik);
            for (int ichi = 0; ichi < nchip_ik; ++ichi)
            {
                sto_eband += kweight * p_che->ddot_real(tmpin, tmpout, npw);
                tmpin += npwx;
                tmpout += npwx;
            }
            delmem_complex_op()(this->ctx, hshchi);
        }
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &sto_eband, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    pes->f_en.eband += sto_eband;
    //---------------------cal rho-------------------------
    double* sto_rho = new double[nrxx];

    double dr3 = GlobalC::ucell.omega / wfc_basis->nxyz;
    double tmprho, tmpne;
    T outtem;
    double sto_ne = 0;
    ModuleBase::GlobalFunc::ZEROS(sto_rho, nrxx);

    T* porter = nullptr;
    resmem_complex_op()(this->ctx, porter, nrxx);
    double out2;

    double* ksrho = nullptr;
    if (PARAM.inp.nbands > 0 && GlobalV::MY_STOGROUP == 0)
    {
        ksrho = new double[nrxx];
        ModuleBase::GlobalFunc::DCOPY(pes->charge->rho[0], ksrho, nrxx);
        setmem_var_op()(this->ctx, pes->rho[0], 0, nrxx);
        // ModuleBase::GlobalFunc::ZEROS(pes->charge->rho[0], nrxx);
    }

    for (int ik = 0; ik < this->pkv->get_nks(); ++ik)
    {
        const int nchip_ik = nchip[ik];
        int current_spin = 0;
        if (PARAM.inp.nspin == 2)
        {
            current_spin = this->pkv->isk[ik];
        }
        stowf.shchi->fix_k(ik);
        T* tmpout = stowf.shchi->get_pointer();
        for (int ichi = 0; ichi < nchip_ik; ++ichi)
        {
            wfc_basis->recip_to_real(this->ctx, tmpout, porter, ik);
            const auto w1 = static_cast<Real>(this->pkv->wk[ik]);
            elecstate::elecstate_pw_op<Real, Device>()(this->ctx, current_spin, nrxx, w1, pes->rho, porter);
            // for (int ir = 0; ir < nrxx; ++ir)
            // {
            //     pes->charge->rho[0][ir] += norm(porter[ir]) * this->pkv->wk[ik];
            // }
            tmpout += npwx;
        }
    }
    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single") {
        for (int ii = 0; ii < PARAM.inp.nspin; ii++) {
            castmem_var_d2h_op()(this->cpu_ctx, this->ctx, pes->charge->rho[ii], pes->rho[ii], nrxx);
        }
    }
    delmem_complex_op()(this->ctx, porter);
#ifdef __MPI
    // temporary, rho_mpi should be rewrite as a tool function! Now it only treats pes->charge->rho
    pes->charge->rho_mpi();
#endif
    for (int ir = 0; ir < nrxx; ++ir)
    {
        tmprho = pes->charge->rho[0][ir] / GlobalC::ucell.omega;
        sto_rho[ir] = tmprho;
        sto_ne += tmprho;
    }
    sto_ne *= dr3;

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &sto_ne, 1, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sto_ne, 1, MPI_DOUBLE, MPI_SUM, PARAPW_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, sto_rho, nrxx, MPI_DOUBLE, MPI_SUM, PARAPW_WORLD);
#endif
    double factor = targetne / (KS_ne + sto_ne);
    if (std::abs(factor - 1) > 1e-10)
    {
        GlobalV::ofs_running << "Renormalize rho from ne = " << sto_ne + KS_ne << " to targetne = " << targetne
                             << std::endl;
    }
    else
    {
        factor = 1;
    }

    if (GlobalV::MY_STOGROUP == 0)
    {
        if (PARAM.inp.nbands > 0)
        {
            ModuleBase::GlobalFunc::DCOPY(ksrho, pes->charge->rho[0], nrxx);
        }
        else
        {
            ModuleBase::GlobalFunc::ZEROS(pes->charge->rho[0], nrxx);
        }
    }

    if (GlobalV::MY_STOGROUP == 0)
    {
        for (int is = 0; is < 1; ++is)
        {
            for (int ir = 0; ir < nrxx; ++ir)
            {
                pes->charge->rho[is][ir] += sto_rho[ir];
                pes->charge->rho[is][ir] *= factor;
            }
        }
    }
    delete[] sto_rho;
    delete[] ksrho;
    ModuleBase::timer::tick("Stochastic_Iter", "sum_stoband");
    return;
}

template <typename T, typename Device>
void Stochastic_Iter<T, Device>::calTnchi_ik(const int& ik, Stochastic_WF<T, Device>& stowf)
{
    const int npw = stowf.ngk[ik];
    const int npwx = stowf.npwx;
    stowf.shchi->fix_k(ik);
    T* out = stowf.shchi->get_pointer();
    T* pchi;
    if (PARAM.inp.nbands > 0)
    {
        stowf.chiortho->fix_k(ik);
        pchi = stowf.chiortho->get_pointer();
    }
    else
    {
        stowf.chi0->fix_k(ik);
        pchi = stowf.chi0->get_pointer();
    }
    if (this->method == 2)
    {
        const char transa = 'N';
        const T one = 1;
        const int inc = 1;
        const T zero = 0;
        const int LDA = npwx * nchip[ik];
        const int M = npwx * nchip[ik];
        const int N = p_che->norder;
        T* coef_real = nullptr;
        resmem_complex_op()(this->ctx, coef_real, N);
        castmem_d2z_op()(this->ctx, this->ctx, coef_real, p_che->coef_real, p_che->norder);
        gemv_op()(this->ctx, transa, M, N, &one, stowf.chiallorder[ik].get_pointer(), LDA, coef_real, inc, &zero, out, inc);
        // zgemv_(&transa, &M, &N, &one, stowf.chiallorder[ik].get_pointer(), &LDA, coef_real, &inc, &zero, out, &inc);
        delmem_complex_op()(this->ctx, coef_real);
    }
    else
    {
        if (this->pkv->get_nks() > 1)
        {
            p_hamilt_sto->updateHk(ik); // necessary, because itermu should be called before this function
        }
        auto hchi_norm = std::bind(&hamilt::HamiltSdftPW<T, Device>::hPsi_norm,
                                   p_hamilt_sto,
                                   std::placeholders::_1,
                                   std::placeholders::_2,
                                   std::placeholders::_3);
        p_che->calfinalvec_real(hchi_norm, pchi, out, npw, npwx, nchip[ik]);
    }
}

template class Stochastic_Iter<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stochastic_Iter<std::complex<double>, base_device::DEVICE_GPU>;
#endif
