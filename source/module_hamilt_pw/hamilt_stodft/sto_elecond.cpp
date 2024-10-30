#include "sto_elecond.h"

#include "module_parameter/parameter.h"
#include "module_base/complexmatrix.h"
#include "module_base/constants.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "sto_tool.h"

#include <chrono>

#define TWOSQRT2LN2 2.354820045030949 // FWHM = 2sqrt(2ln2) * \sigma
#define FACTOR 1.839939223835727e7

Sto_EleCond::Sto_EleCond(UnitCell* p_ucell_in,
                         K_Vectors* p_kv_in,
                         elecstate::ElecState* p_elec_in,
                         ModulePW::PW_Basis_K* p_wfcpw_in,
                         psi::Psi<std::complex<double>>* p_psi_in,
                         pseudopot_cell_vnl* p_ppcell_in,
                         hamilt::Hamilt<std::complex<double>>* p_hamilt_in,
                         StoChe<double>& stoche,
                         Stochastic_WF* p_stowf_in)
    : EleCond(p_ucell_in, p_kv_in, p_elec_in, p_wfcpw_in, p_psi_in, p_ppcell_in)
{
    this->p_hamilt = p_hamilt_in;
    this->p_hamilt_sto = static_cast<hamilt::HamiltSdftPW<std::complex<double>>*>(p_hamilt_in);
    this->p_stowf = p_stowf_in;
    this->nbands_ks = p_psi_in->get_nbands();
    this->nbands_sto = p_stowf_in->nchi;
    this->stofunc.set_E_range(&stoche.emin_sto, &stoche.emax_sto);
}

void Sto_EleCond::decide_nche(const double dt,
                              const double cond_thr,
                              const int& fd_nche,
                              double try_emin,
                              double try_emax)
{
    int nche_guess = 1000;
    ModuleBase::Chebyshev<double> chet(nche_guess);
    const double mu = this->p_elec->eferm.ef;
    this->stofunc.mu = mu;
    int& nbatch = this->cond_dtbatch;
    auto ncos = std::bind(&Sto_Func<double>::ncos, &this->stofunc, std::placeholders::_1);
    auto n_sin = std::bind(&Sto_Func<double>::n_sin, &this->stofunc, std::placeholders::_1);
    // try to find nbatch
    if (nbatch == 0)
    {
        for (int test_nbatch = 128; test_nbatch >= 1; test_nbatch /= 2)
        {
            nbatch = test_nbatch;
            this->stofunc.t = 0.5 * dt * nbatch;
            chet.calcoef_pair(ncos, n_sin);
            double minerror = std::abs(chet.coef_complex[nche_guess - 1] / chet.coef_complex[0]);
            if (minerror < cond_thr)
            {
                for (int i = 1; i < nche_guess; ++i)
                {
                    double error = std::abs(chet.coef_complex[i] / chet.coef_complex[0]);
                    if (error < cond_thr)
                    {
                        // To make nche to around 100
                        nbatch = ceil(double(test_nbatch) / i * 100.0);
                        std::cout << "set cond_dtbatch to " << nbatch << std::endl;
                        break;
                    }
                }
                break;
            }
        }
    }

    // first try to find nche
    this->stofunc.t = 0.5 * dt * nbatch;
    auto getnche = [&](int& nche) {
        chet.calcoef_pair(ncos, n_sin);
        for (int i = 1; i < nche_guess; ++i)
        {
            double error = std::abs(chet.coef_complex[i] / chet.coef_complex[0]);
            if (error < cond_thr)
            {
                nche = i + 1;
                break;
            }
        }
    };
    int nche_old = 0;
    getnche(nche_old);

    int nche_new = 0;
loop:
    // re-set Emin & Emax both in p_hamilt_sto & stofunc
    check_che(std::max(nche_old * 2, fd_nche),
              try_emin,
              try_emax,
              this->nbands_sto,
              this->p_kv,
              this->p_stowf,
              this->p_hamilt_sto);

    // second try to find nche with new Emin & Emax
    getnche(nche_new);

    if (nche_new > nche_old * 2)
    {
        nche_old = nche_new;
        try_emin = *p_hamilt_sto->emin;
        try_emax = *p_hamilt_sto->emax;
        goto loop;
    }

    std::cout << "set N order of Chebyshev for KG as " << nche_new << std::endl;
    std::ofstream cheofs("Chebycoef");
    for (int i = 1; i < nche_guess; ++i)
    {
        double error = std::abs(chet.coef_complex[i] / chet.coef_complex[0]);
        cheofs << std::setw(5) << i << std::setw(20) << error << std::endl;
    }
    cheofs.close();

    if (nche_new >= 1000)
    {
        ModuleBase::WARNING_QUIT("ESolver_SDFT_PW", "N order of Chebyshev for KG will be larger than 1000!");
    }

    this->cond_nche = nche_new;
    this->fd_nche = fd_nche;
}

void Sto_EleCond::cal_jmatrix(const psi::Psi<std::complex<float>>& kspsi_all,
                              const psi::Psi<std::complex<float>>& vkspsi,
                              const double* en,
                              const double* en_all,
                              std::complex<double>* leftfact,
                              std::complex<double>* rightfact,
                              const psi::Psi<std::complex<double>>& leftchi,
                              psi::Psi<std::complex<double>>& rightchi,
                              psi::Psi<std::complex<double>>& left_hchi,
                              psi::Psi<std::complex<double>>& batch_vchi,
                              psi::Psi<std::complex<double>>& batch_vhchi,
#ifdef __MPI
                              psi::Psi<std::complex<float>>& chi_all,
                              psi::Psi<std::complex<float>>& hchi_all,
                              void* gatherinfo_ks,
                              void* gatherinfo_sto,
#endif
                              const int& bsize_psi,
                              std::vector<std::complex<float>>& j1,
                              std::vector<std::complex<float>>& j2,
                              hamilt::Velocity& velop,
                              const int& ik,
                              const std::complex<double>& factor,
                              const int bandinfo[6])
{
    ModuleBase::timer::tick("Sto_EleCond", "cal_jmatrix");
    const char transa = 'C';
    const char transb = 'N';
    const std::complex<float> float_factor = static_cast<std::complex<float>>(factor);
    const std::complex<float> conjfactor = std::conj(float_factor);
    const float mu = static_cast<float>(this->p_elec->eferm.ef);
    const std::complex<float> zero = 0.0;
    const int npwx = this->p_wfcpw->npwk_max;
    const int npw = this->p_wfcpw->npwk[ik];
    const int ndim = 3;
    const int perbands_ks = bandinfo[0];
    const int perbands_sto = bandinfo[1];
    const int perbands = bandinfo[2];
    const int allbands_ks = bandinfo[3];
    const int allbands_sto = bandinfo[4];
    const int allbands = bandinfo[5];
    const int dim_jmatrix = perbands_ks * allbands_sto + perbands_sto * allbands;

    psi::Psi<std::complex<double>> right_hchi(1, perbands_sto, npwx, p_kv->ngk.data());
    psi::Psi<std::complex<float>> f_rightchi(1, perbands_sto, npwx, p_kv->ngk.data());
    psi::Psi<std::complex<float>> f_right_hchi(1, perbands_sto, npwx, p_kv->ngk.data());

    this->p_hamilt_sto->hPsi(leftchi.get_pointer(), left_hchi.get_pointer(), perbands_sto);
    this->p_hamilt_sto->hPsi(rightchi.get_pointer(), right_hchi.get_pointer(), perbands_sto);
    convert_psi(rightchi, f_rightchi);
    convert_psi(right_hchi, f_right_hchi);
    right_hchi.resize(1, 1, 1);

    psi::Psi<std::complex<float>>* rightchi_all = &f_rightchi;
    psi::Psi<std::complex<float>>* righthchi_all = &f_right_hchi;
    std::vector<std::complex<double>> vec_rightf_all;
    std::complex<double>* rightf_all = rightfact;
#ifdef __MPI
    info_gatherv* ks_fact = static_cast<info_gatherv*>(gatherinfo_ks);
    info_gatherv* sto_npwx = static_cast<info_gatherv*>(gatherinfo_sto);
    rightchi_all = gatherchi(f_rightchi, chi_all, npwx, sto_npwx->nrecv, sto_npwx->displs, perbands_sto);
    righthchi_all = gatherchi(f_right_hchi, hchi_all, npwx, sto_npwx->nrecv, sto_npwx->displs, perbands_sto);
    if (PARAM.inp.bndpar > 1 && rightfact != nullptr)
    {
        vec_rightf_all.resize(allbands_ks);
        rightf_all = vec_rightf_all.data();
        MPI_Allgatherv(rightfact,
                       perbands_ks,
                       MPI_DOUBLE_COMPLEX,
                       rightf_all,
                       ks_fact->nrecv,
                       ks_fact->displs,
                       MPI_DOUBLE_COMPLEX,
                       PARAPW_WORLD);
    }
#endif

    psi::Psi<std::complex<float>> f_batch_vchi(1, bsize_psi * ndim, npwx, p_kv->ngk.data());
    psi::Psi<std::complex<float>> f_batch_vhchi(1, bsize_psi * ndim, npwx, p_kv->ngk.data());
    std::vector<std::complex<float>> tmpj(ndim * allbands_sto * perbands_sto);

    // 1. (<\psi|J|\chi>)^T
    // (allbands_sto, perbands_ks)
    if (perbands_ks > 0)
    {
        for (int id = 0; id < ndim; ++id)
        {
            const int idnb = id * perbands_ks;
            const int jbais = 0;
            std::complex<float>* j1mat = &j1[id * dim_jmatrix];
            std::complex<float>* j2mat = &j2[id * dim_jmatrix];
            //(<\psi|v|\chi>)^T
            cgemm_(&transa,
                   &transb,
                   &allbands_sto,
                   &perbands_ks,
                   &npw,
                   &conjfactor,
                   rightchi_all->get_pointer(),
                   &npwx,
                   &vkspsi(idnb, 0),
                   &npwx,
                   &zero,
                   j1mat,
                   &allbands_sto);

            //(<\psi|Hv|\chi>)^T
            // for(int i = 0 ; i < perbands_ks ; ++i)
            // {
            //     double* evmat = &j2(id, 0 + i * allbands_sto);
            //     double* vmat = &j1(id, 0 + i * allbands_sto);
            //     double ei = en[i];
            //     for(int j = 0 ; j < allbands_sto ; ++j)
            //     {
            //         evmat[j] = vmat[j] * ei;
            //     }
            // }

            //(<\psi|vH|\chi>)^T
            cgemm_(&transa,
                   &transb,
                   &allbands_sto,
                   &perbands_ks,
                   &npw,
                   &conjfactor,
                   righthchi_all->get_pointer(),
                   &npwx,
                   &vkspsi(idnb, 0),
                   &npwx,
                   &zero,
                   j2mat,
                   &allbands_sto);
        }
    }

    int remain = perbands_sto;
    int startnb = 0;
    while (remain > 0)
    {
        int tmpnb = std::min(remain, bsize_psi);
        // v|\chi>
        velop.act(&leftchi, tmpnb, &leftchi(0, startnb, 0), batch_vchi.get_pointer());
        convert_psi(batch_vchi, f_batch_vchi);
        // v|H\chi>
        velop.act(&leftchi, tmpnb, &left_hchi(0, startnb, 0), batch_vhchi.get_pointer());
        convert_psi(batch_vhchi, f_batch_vhchi);
        // 2. <\chi|J|\psi>
        // (perbands_sto, allbands_ks)
        if (allbands_ks > 0)
        {
            for (int id = 0; id < ndim; ++id)
            {
                const int idnb = id * tmpnb;
                const int jbais = perbands_ks * allbands_sto + startnb;
                std::complex<float>* j1mat = &j1[id * dim_jmatrix + jbais];
                std::complex<float>* j2mat = &j2[id * dim_jmatrix + jbais];
                //<\chi|v|\psi>
                cgemm_(&transa,
                       &transb,
                       &tmpnb,
                       &allbands_ks,
                       &npw,
                       &float_factor,
                       &f_batch_vchi(idnb, 0),
                       &npwx,
                       kspsi_all.get_pointer(),
                       &npwx,
                       &zero,
                       j1mat,
                       &perbands_sto);

                //<\chi|vH|\psi> = \epsilon * <\chi|v|\psi>
                // for(int i = 0 ; i < allbands_ks ; ++i)
                // {
                //     double* evmat = &j2(id, jbais + i * allbands_sto);
                //     double* vmat = &j1(id, jbais + i * allbands_sto);
                //     double ei = en_all[i];
                //     for(int j = 0 ; j < tmpnb ; ++j)
                //     {
                //         evmat[j] = vmat[j] * en_all[i];
                //     }
                // }

                //<\chi|Hv|\psi>
                cgemm_(&transa,
                       &transb,
                       &tmpnb,
                       &allbands_ks,
                       &npw,
                       &float_factor,
                       &f_batch_vhchi(idnb, 0),
                       &npwx,
                       kspsi_all.get_pointer(),
                       &npwx,
                       &zero,
                       j2mat,
                       &perbands_sto);
            }
        }

        // 3. <\chi|J|\chi>
        // (perbands_sto, allbands_sto)
        for (int id = 0; id < ndim; ++id)
        {
            const int idnb = id * tmpnb;
            const int jbais = perbands_ks * allbands_sto + perbands_sto * allbands_ks + startnb;
            std::complex<float>* j1mat = &j1[id * dim_jmatrix + jbais];
            std::complex<float>* j2mat = &j2[id * dim_jmatrix + jbais];
            std::complex<float>* tmpjmat = &tmpj[id * allbands_sto * perbands_sto + startnb];
            //<\chi|v|\chi>
            cgemm_(&transa,
                   &transb,
                   &tmpnb,
                   &allbands_sto,
                   &npw,
                   &float_factor,
                   &f_batch_vchi(idnb, 0),
                   &npwx,
                   rightchi_all->get_pointer(),
                   &npwx,
                   &zero,
                   j1mat,
                   &perbands_sto);

            //<\chi|Hv|\chi>
            cgemm_(&transa,
                   &transb,
                   &tmpnb,
                   &allbands_sto,
                   &npw,
                   &float_factor,
                   &f_batch_vhchi(idnb, 0),
                   &npwx,
                   rightchi_all->get_pointer(),
                   &npwx,
                   &zero,
                   j2mat,
                   &perbands_sto);

            //<\chi|vH|\chi>
            cgemm_(&transa,
                   &transb,
                   &tmpnb,
                   &allbands_sto,
                   &npw,
                   &float_factor,
                   &f_batch_vchi(idnb, 0),
                   &npwx,
                   righthchi_all->get_pointer(),
                   &npwx,
                   &zero,
                   tmpjmat,
                   &perbands_sto);
        }

        remain -= tmpnb;
        startnb += tmpnb;
        if (remain == 0) {
            break;
}
    }

    for (int id = 0; id < ndim; ++id)
    {
        for (int i = 0; i < perbands_ks; ++i)
        {
            const float ei = static_cast<float>(en[i]);
            const int jst = i * allbands_sto;
            std::complex<float>* j2mat = j2.data() + id * dim_jmatrix + jst;
            std::complex<float>* j1mat = j1.data() + id * dim_jmatrix + jst;
            if (leftfact == nullptr)
            {
                for (int j = 0; j < allbands_sto; ++j)
                {
                    j2mat[j] = 0.5f * j2mat[j] + (0.5f * ei - mu) * j1mat[j];
                }
            }
            else
            {
                const std::complex<float> jfac = static_cast<std::complex<float>>(leftfact[i]);
                for (int j = 0; j < allbands_sto; ++j)
                {
                    j2mat[j] = jfac * (0.5f * j2mat[j] + (0.5f * ei - mu) * j1mat[j]);
                    j1mat[j] *= jfac;
                }
            }
        }

        for (int i = 0; i < allbands_ks; ++i)
        {
            const float ei = static_cast<float>(en_all[i]);
            const int jst = perbands_ks * allbands_sto + i * perbands_sto;
            std::complex<float>* j2mat = j2.data() + id * dim_jmatrix + jst;
            std::complex<float>* j1mat = j1.data() + id * dim_jmatrix + jst;
            if (rightfact == nullptr)
            {
                for (int j = 0; j < perbands_sto; ++j)
                {
                    j2mat[j] = 0.5f * j2mat[j] + (0.5f * ei - mu) * j1mat[j];
                }
            }
            else
            {
                const std::complex<float> jfac = static_cast<std::complex<float>>(rightf_all[i]);
                for (int j = 0; j < perbands_sto; ++j)
                {
                    j2mat[j] = jfac * (0.5f * j2mat[j] + (0.5f * ei - mu) * j1mat[j]);
                    j1mat[j] *= jfac;
                }
            }
        }

        const int jst = perbands_ks * allbands_sto + perbands_sto * allbands_ks;
        const int ed = dim_jmatrix - jst;
        std::complex<float>* j2mat = j2.data() + id * dim_jmatrix + jst;
        std::complex<float>* j1mat = j1.data() + id * dim_jmatrix + jst;
        std::complex<float>* tmpjmat = tmpj.data() + id * allbands_sto * perbands_sto;

        for (int j = 0; j < ed; ++j)
        {
            j2mat[j] = 0.5f * (j2mat[j] + tmpjmat[j]) - mu * j1mat[j];
        }
    }

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, j1.data(), ndim * dim_jmatrix, MPI_COMPLEX, MPI_SUM, POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, j2.data(), ndim * dim_jmatrix, MPI_COMPLEX, MPI_SUM, POOL_WORLD);
#endif

    ModuleBase::timer::tick("Sto_EleCond", "cal_jmatrix");

    return;
}

void Sto_EleCond::sKG(const int& smear_type,
                      const double& fwhmin,
                      const double& wcut,
                      const double& dw_in,
                      const double& dt_in,
                      const bool& nonlocal,
                      const int& npart_sto)
{
    ModuleBase::TITLE("Sto_EleCond", "sKG");
    ModuleBase::timer::tick("Sto_EleCond", "sKG");
    std::cout << "Calculating conductivity...." << std::endl;
    // if (PARAM.inp.bndpar > 1)
    // {
    //     ModuleBase::WARNING_QUIT("ESolver_SDFT_PW", "sKG is not supported in parallel!");
    // }

    //------------------------------------------------------------------
    //                    Init
    //------------------------------------------------------------------
    // Parameters
    const int nbatch = this->cond_dtbatch;
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    double gamma = fwhmin / 2.0 / ModuleBase::Ry_to_eV;
    double dt = dt_in;              // unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    const double expfactor = 18.42; // exp(-18.42) = 1e-8
    int nt = 0;                     // set nt empirically
    if (smear_type == 1)
    {
        nt = ceil(sqrt(2 * expfactor) / sigma / dt);
    }
    else if (smear_type == 2)
    {
        nt = ceil(expfactor / gamma / dt);
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW::calcondw", "smear_type should be 0 or 1");
    }
    std::cout << "nw: " << nw << " ; dw: " << dw * ModuleBase::Ry_to_eV << " eV" << std::endl;
    std::cout << "nt: " << nt << " ; dt: " << dt << " a.u.(ry^-1)" << std::endl;
    assert(nw >= 1);
    assert(nt >= 1);
    const int ndim = 3;
    const int nk = p_kv->get_nks();
    const int npwx = p_wfcpw->npwk_max;
    const double tpiba = p_wfcpw->tpiba;
    psi::Psi<std::complex<double>>* stopsi;
    if (this->nbands_ks > 0)
    {
        stopsi = p_stowf->chiortho;
        // clean memories //Note shchi is different from \sqrt(fH_here)|chi>, since veffs are different
        p_stowf->shchi->resize(1, 1, 1);
        p_stowf->chi0->resize(1, 1, 1); // clean memories
    }
    else
    {
        stopsi = p_stowf->chi0;
        p_stowf->shchi->resize(1, 1, 1); // clean memories
    }
    const double dEcut = (wcut + fwhmin) / ModuleBase::Ry_to_eV;

    // response funtion
    std::vector<double> ct11(nt, 0);
    std::vector<double> ct12(nt, 0);
    std::vector<double> ct22(nt, 0);

    // Init Chebyshev
    ModuleBase::Chebyshev<double> che(fd_nche);
    ModuleBase::Chebyshev<double> chet(cond_nche);
    ModuleBase::Chebyshev<double> chemt(cond_nche);

    //------------------------------------------------------------------
    //                    Calculate
    //------------------------------------------------------------------

    // Prepare Chebyshev coefficients for exp(-i H/\hbar t)
    const double mu = this->p_elec->eferm.ef;
    this->stofunc.mu = mu;
    this->stofunc.t = 0.5 * dt * nbatch;
    auto ncos = std::bind(&Sto_Func<double>::ncos, &this->stofunc, std::placeholders::_1);
    auto nsin = std::bind(&Sto_Func<double>::nsin, &this->stofunc, std::placeholders::_1);
    auto n_sin = std::bind(&Sto_Func<double>::n_sin, &this->stofunc, std::placeholders::_1);
    chet.calcoef_pair(ncos, nsin);
    chemt.calcoef_pair(ncos, n_sin);
    std::vector<std::complex<double>> batchcoef, batchmcoef;
    if (nbatch > 1)
    {
        batchcoef.resize(cond_nche * nbatch);
        std::complex<double>* tmpcoef = batchcoef.data() + (nbatch - 1) * cond_nche;
        batchmcoef.resize(cond_nche * nbatch);
        std::complex<double>* tmpmcoef = batchmcoef.data() + (nbatch - 1) * cond_nche;
        for (int i = 0; i < cond_nche; ++i)
        {
            tmpcoef[i] = chet.coef_complex[i];
            tmpmcoef[i] = chemt.coef_complex[i];
        }
        for (int ib = 0; ib < nbatch - 1; ++ib)
        {
            tmpcoef = batchcoef.data() + ib * cond_nche;
            tmpmcoef = batchmcoef.data() + ib * cond_nche;
            this->stofunc.t = 0.5 * dt * (ib + 1);
            chet.calcoef_pair(ncos, nsin);
            chemt.calcoef_pair(ncos, n_sin);
            for (int i = 0; i < cond_nche; ++i)
            {
                tmpcoef[i] = chet.coef_complex[i];
                tmpmcoef[i] = chemt.coef_complex[i];
            }
        }
        this->stofunc.t = 0.5 * dt * nbatch;
    }

    // ik loop
    ModuleBase::timer::tick("Sto_EleCond", "kloop");
    hamilt::Velocity velop(p_wfcpw, p_kv->isk.data(), p_ppcell, p_ucell, nonlocal);
    for (int ik = 0; ik < nk; ++ik)
    {
        velop.init(ik);
        stopsi->fix_k(ik);
        p_psi->fix_k(ik);
        if (nk > 1)
        {
            this->p_hamilt->updateHk(ik);
        }
        const int npw = p_kv->ngk[ik];

        // get allbands_ks
        int cutib0 = 0;
        if (this->nbands_ks > 0)
        {
            double Emax_KS = std::max(*this->stofunc.Emin, this->p_elec->ekb(ik, this->nbands_ks - 1));
            for (cutib0 = this->nbands_ks - 1; cutib0 >= 0; --cutib0)
            {
                if (Emax_KS - this->p_elec->ekb(ik, cutib0) > dEcut)
                {
                    break;
                }
            }
            ++cutib0;
            double Emin_KS = (cutib0 < this->nbands_ks) ? this->p_elec->ekb(ik, cutib0) : *this->stofunc.Emin;
            double dE = *this->stofunc.Emax - Emin_KS + wcut / ModuleBase::Ry_to_eV;
            std::cout << "Emin_KS(" << cutib0 + 1 << "): " << Emin_KS * ModuleBase::Ry_to_eV
                      << " eV; Emax: " << *this->stofunc.Emax * ModuleBase::Ry_to_eV
                      << " eV; Recommended max dt: " << 2 * M_PI / dE << " a.u." << std::endl;
        }
        else
        {
            double dE = *this->stofunc.Emax - *this->stofunc.Emin + wcut / ModuleBase::Ry_to_eV;
            std::cout << "Emin: " << *this->stofunc.Emin * ModuleBase::Ry_to_eV
                      << " eV; Emax: " << *this->stofunc.Emax * ModuleBase::Ry_to_eV
                      << " eV; Recommended max dt: " << 2 * M_PI / dE << " a.u." << std::endl;
        }
        // Parallel for bands
        int allbands_ks = this->nbands_ks - cutib0;
        parallel_distribution paraks(allbands_ks, PARAM.inp.bndpar, GlobalV::MY_STOGROUP);
        int perbands_ks = paraks.num_per;
        int ib0_ks = paraks.start;
        ib0_ks += this->nbands_ks - allbands_ks;
        int perbands_sto = this->p_stowf->nchip[ik];
        int perbands = perbands_sto + perbands_ks;
        int allbands_sto = perbands_sto;
        int allbands = perbands;
#ifdef __MPI
        MPI_Allreduce(&perbands, &allbands, 1, MPI_INT, MPI_SUM, PARAPW_WORLD);
        allbands_sto = allbands - allbands_ks;
        info_gatherv ks_fact(perbands_ks, PARAM.inp.bndpar, 1, PARAPW_WORLD);
        info_gatherv sto_npwx(perbands_sto, PARAM.inp.bndpar, npwx, PARAPW_WORLD);
#endif
        const int bandsinfo[6]{perbands_ks, perbands_sto, perbands, allbands_ks, allbands_sto, allbands};
        double* en_all = nullptr;
        std::vector<double> en;
        if (allbands_ks > 0)
        {
            en_all = &(this->p_elec->ekb(ik, this->nbands_ks - allbands_ks));
        }
        if (perbands_ks > 0)
        {
            en.resize(perbands_ks);
            for (int ib = 0; ib < perbands_ks; ++ib)
            {
                en[ib] = this->p_elec->ekb(ik, ib0_ks + ib);
            }
        }

        //-----------------------------------------------------------
        //               ks conductivity
        //-----------------------------------------------------------
        if (GlobalV::MY_STOGROUP == 0 && allbands_ks > 0)
        {
            jjresponse_ks(ik, nt, dt, dEcut, this->p_elec->wg, velop, ct11.data(), ct12.data(), ct22.data());
        }

        //-----------------------------------------------------------
        //               sto conductivity
        //-----------------------------------------------------------
        //-------------------     allocate  -------------------------
        size_t ks_memory_cost = perbands_ks * npwx * sizeof(std::complex<float>);
        psi::Psi<std::complex<double>> kspsi(1, perbands_ks, npwx, p_kv->ngk.data());
        psi::Psi<std::complex<double>> vkspsi(1, perbands_ks * ndim, npwx, p_kv->ngk.data());
        std::vector<std::complex<double>> expmtmf_fact(perbands_ks), expmtf_fact(perbands_ks);
        psi::Psi<std::complex<float>> f_kspsi(1, perbands_ks, npwx, p_kv->ngk.data());
        ModuleBase::Memory::record("SDFT::kspsi", ks_memory_cost);
        psi::Psi<std::complex<float>> f_vkspsi(1, perbands_ks * ndim, npwx, p_kv->ngk.data());
        ModuleBase::Memory::record("SDFT::vkspsi", ks_memory_cost);
        psi::Psi<std::complex<float>>* kspsi_all = &f_kspsi;

        size_t sto_memory_cost = perbands_sto * npwx * sizeof(std::complex<double>);
        psi::Psi<std::complex<double>> sfchi(1, perbands_sto, npwx, p_kv->ngk.data());
        ModuleBase::Memory::record("SDFT::sfchi", sto_memory_cost);
        psi::Psi<std::complex<double>> smfchi(1, perbands_sto, npwx, p_kv->ngk.data());
        ModuleBase::Memory::record("SDFT::smfchi", sto_memory_cost);
#ifdef __MPI
        psi::Psi<std::complex<float>> chi_all, hchi_all, psi_all;
        if (PARAM.inp.bndpar > 1)
        {
            chi_all.resize(1, allbands_sto, npwx);
            hchi_all.resize(1, allbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::chi_all", allbands_sto * npwx * sizeof(std::complex<float>));
            ModuleBase::Memory::record("SDFT::hchi_all", allbands_sto * npwx * sizeof(std::complex<float>));
            psi_all.resize(1, allbands_ks, npwx);
            ModuleBase::Memory::record("SDFT::kspsi_all", allbands_ks * npwx * sizeof(std::complex<float>));
            for (int ib = 0; ib < allbands_ks; ++ib)
            {
                for (int ig = 0; ig < npw; ++ig)
                {
                    psi_all(0, ib, ig)
                        = static_cast<std::complex<float>>(p_psi[0](this->nbands_ks - allbands_ks + ib, ig));
                }
            }
            kspsi_all = &psi_all;
            f_kspsi.resize(1, 1, 1);
        }
#endif

        const int nbatch_psi = npart_sto;
        const int bsize_psi = ceil(double(perbands_sto) / nbatch_psi);
        psi::Psi<std::complex<double>> batch_vchi(1, bsize_psi * ndim, npwx, p_kv->ngk.data());
        psi::Psi<std::complex<double>> batch_vhchi(1, bsize_psi * ndim, npwx, p_kv->ngk.data());
        ModuleBase::Memory::record("SDFT::batchjpsi", 3 * bsize_psi * ndim * npwx * sizeof(std::complex<double>));

        //-------------------     sqrt(f)|psi>   sqrt(1-f)|psi>   ---------------
        if (perbands_ks > 0)
        {
            for (int ib = 0; ib < perbands_ks; ++ib)
            {
                for (int ig = 0; ig < npw; ++ig)
                {
                    kspsi(0, ib, ig) = p_psi[0](ib0_ks + ib, ig);
                }
                double fi = this->stofunc.fd(en[ib]);
                expmtmf_fact[ib] = 1 - fi;
                expmtf_fact[ib] = fi;
            }
            // v|\psi>
            velop.act(&kspsi, perbands_ks, kspsi.get_pointer(), vkspsi.get_pointer());
            // convert to complex<float>
            if (PARAM.inp.bndpar == 1)
            {
                convert_psi(kspsi, f_kspsi);
            }
            convert_psi(vkspsi, f_vkspsi);
            kspsi.resize(1, 1, 1);
            vkspsi.resize(1, 1, 1);
        }

        auto nroot_fd = std::bind(&Sto_Func<double>::nroot_fd, &this->stofunc, std::placeholders::_1);
        che.calcoef_real(nroot_fd);
        auto hchi_norm = std::bind(&hamilt::HamiltSdftPW<std::complex<double>>::hPsi_norm,
                                   p_hamilt_sto,
                                   std::placeholders::_1,
                                   std::placeholders::_2,
                                   std::placeholders::_3);
        che.calfinalvec_real(hchi_norm, stopsi->get_pointer(), sfchi.get_pointer(), npw, npwx, perbands_sto);

        auto nroot_mfd = std::bind(&Sto_Func<double>::nroot_mfd, &this->stofunc, std::placeholders::_1);
        che.calcoef_real(nroot_mfd);

        che.calfinalvec_real(hchi_norm, stopsi->get_pointer(), smfchi.get_pointer(), npw, npwx, perbands_sto);

        //------------------------  allocate ------------------------
        psi::Psi<std::complex<double>>& expmtsfchi = sfchi;
        psi::Psi<std::complex<double>>& expmtsmfchi = smfchi;
        psi::Psi<std::complex<double>> exptsfchi = expmtsfchi;
        ModuleBase::Memory::record("SDFT::exptsfchi", sto_memory_cost);
        psi::Psi<std::complex<double>> exptsmfchi = expmtsmfchi;
        ModuleBase::Memory::record("SDFT::exptsmfchi", sto_memory_cost);
        psi::Psi<std::complex<double>> poly_expmtsfchi, poly_expmtsmfchi;
        psi::Psi<std::complex<double>> poly_exptsfchi, poly_exptsmfchi;
        if (nbatch > 1)
        {
            poly_exptsfchi.resize(cond_nche, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_exptsfchi",
                                       sizeof(std::complex<double>) * cond_nche * perbands_sto * npwx);

            poly_exptsmfchi.resize(cond_nche, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_exptsmfchi",
                                       sizeof(std::complex<double>) * cond_nche * perbands_sto * npwx);

            poly_expmtsfchi.resize(cond_nche, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_expmtsfchi",
                                       sizeof(std::complex<double>) * cond_nche * perbands_sto * npwx);

            poly_expmtsmfchi.resize(cond_nche, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_expmtsmfchi",
                                       sizeof(std::complex<double>) * cond_nche * perbands_sto * npwx);
        }

        const int dim_jmatrix = perbands_ks * allbands_sto + perbands_sto * allbands;
        parallel_distribution parajmat(ndim * dim_jmatrix, GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL);
        std::vector<std::complex<float>> j1l(ndim * dim_jmatrix), j2l(ndim * dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j1l", sizeof(std::complex<float>) * ndim * dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j2l", sizeof(std::complex<float>) * ndim * dim_jmatrix);
        std::vector<std::complex<float>> j1r(ndim * dim_jmatrix), j2r(ndim * dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j1r", sizeof(std::complex<float>) * ndim * dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j2r", sizeof(std::complex<float>) * ndim * dim_jmatrix);
        psi::Psi<std::complex<double>> tmphchil(1, perbands_sto, npwx, p_kv->ngk.data());
        ModuleBase::Memory::record("SDFT::tmphchil/r", sto_memory_cost * 2);

        //------------------------  t loop  --------------------------
        std::cout << "ik=" << ik << ": ";
        auto start = std::chrono::high_resolution_clock::now();
        const int print_step = ceil(20.0 / nbatch) * nbatch;
        for (int it = 1; it < nt; ++it)
        {
            // evaluate time cost
            if (it - 1 == print_step)
            {
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                double timeTaken = duration.count();
                std::cout << "(Time left " << timeTaken * (double(nt - 1) / print_step * (nk - ik) - 1) << " s) "
                          << std::endl;
                std::cout << "nt: " << std::endl;
            }
            if ((it - 1) % print_step == 0 && it > 1)
            {
                std::cout << std::setw(8) << it - 1;
                if ((it - 1) % (print_step * 10) == 0)
                {
                    std::cout << std::endl;
                }
            }

            // time evolution exp(-iHt)|\psi_ks>
            // KS
            ModuleBase::timer::tick("Sto_EleCond", "evolution");
            for (int ib = 0; ib < perbands_ks; ++ib)
            {
                double eigen = en[ib];
                const std::complex<double> expmfactor = exp(ModuleBase::NEG_IMAG_UNIT * eigen * dt);
                expmtf_fact[ib] *= expmfactor;
                expmtmf_fact[ib] *= expmfactor;
            }
            // Sto
            if (nbatch == 1)
            {
                chemt.calfinalvec_complex(hchi_norm,
                                          expmtsfchi.get_pointer(),
                                          expmtsfchi.get_pointer(),
                                          npw,
                                          npwx,
                                          perbands_sto);
                chemt.calfinalvec_complex(hchi_norm,
                                          expmtsmfchi.get_pointer(),
                                          expmtsmfchi.get_pointer(),
                                          npw,
                                          npwx,
                                          perbands_sto);
                chet.calfinalvec_complex(hchi_norm,
                                         exptsfchi.get_pointer(),
                                         exptsfchi.get_pointer(),
                                         npw,
                                         npwx,
                                         perbands_sto);
                chet.calfinalvec_complex(hchi_norm,
                                         exptsmfchi.get_pointer(),
                                         exptsmfchi.get_pointer(),
                                         npw,
                                         npwx,
                                         perbands_sto);
            }
            else
            {
                std::complex<double>* tmppolyexpmtsfchi = poly_expmtsfchi.get_pointer();
                std::complex<double>* tmppolyexpmtsmfchi = poly_expmtsmfchi.get_pointer();
                std::complex<double>* tmppolyexptsfchi = poly_exptsfchi.get_pointer();
                std::complex<double>* tmppolyexptsmfchi = poly_exptsmfchi.get_pointer();
                std::complex<double>* stoexpmtsfchi = expmtsfchi.get_pointer();
                std::complex<double>* stoexpmtsmfchi = expmtsmfchi.get_pointer();
                std::complex<double>* stoexptsfchi = exptsfchi.get_pointer();
                std::complex<double>* stoexptsmfchi = exptsmfchi.get_pointer();
                if ((it - 1) % nbatch == 0)
                {
                    chet.calpolyvec_complex(hchi_norm, stoexptsfchi, tmppolyexptsfchi, npw, npwx, perbands_sto);
                    chet.calpolyvec_complex(hchi_norm, stoexptsmfchi, tmppolyexptsmfchi, npw, npwx, perbands_sto);
                    chemt.calpolyvec_complex(hchi_norm, stoexpmtsfchi, tmppolyexpmtsfchi, npw, npwx, perbands_sto);
                    chemt.calpolyvec_complex(hchi_norm, stoexpmtsmfchi, tmppolyexpmtsmfchi, npw, npwx, perbands_sto);
                }

                std::complex<double>* tmpcoef = batchcoef.data() + (it - 1) % nbatch * cond_nche;
                std::complex<double>* tmpmcoef = batchmcoef.data() + (it - 1) % nbatch * cond_nche;
                const char transa = 'N';
                const int LDA = perbands_sto * npwx;
                const int M = perbands_sto * npwx;
                const int N = cond_nche;
                const int inc = 1;
                zgemv_(&transa,
                       &M,
                       &N,
                       &ModuleBase::ONE,
                       tmppolyexptsfchi,
                       &LDA,
                       tmpcoef,
                       &inc,
                       &ModuleBase::ZERO,
                       stoexptsfchi,
                       &inc);
                zgemv_(&transa,
                       &M,
                       &N,
                       &ModuleBase::ONE,
                       tmppolyexptsmfchi,
                       &LDA,
                       tmpcoef,
                       &inc,
                       &ModuleBase::ZERO,
                       stoexptsmfchi,
                       &inc);
                zgemv_(&transa,
                       &M,
                       &N,
                       &ModuleBase::ONE,
                       tmppolyexpmtsfchi,
                       &LDA,
                       tmpmcoef,
                       &inc,
                       &ModuleBase::ZERO,
                       stoexpmtsfchi,
                       &inc);
                zgemv_(&transa,
                       &M,
                       &N,
                       &ModuleBase::ONE,
                       tmppolyexpmtsmfchi,
                       &LDA,
                       tmpmcoef,
                       &inc,
                       &ModuleBase::ZERO,
                       stoexpmtsmfchi,
                       &inc);
            }
            ModuleBase::timer::tick("Sto_EleCond", "evolution");

            // calculate i<\psi|sqrt(f) exp(-iHt/2)*J*exp(iHt/2) sqrt(1-f)|\psi>^+
            //         = i<\psi|sqrt(1-f) exp(-iHt/2)*J*exp(iHt/2) sqrt(f)|\psi>
            cal_jmatrix(*kspsi_all,
                        f_vkspsi,
                        en.data(),
                        en_all,
                        nullptr,
                        nullptr,
                        exptsmfchi,
                        exptsfchi,
                        tmphchil,
                        batch_vchi,
                        batch_vhchi,
#ifdef __MPI
                        chi_all,
                        hchi_all,
                        (void*)&ks_fact,
                        (void*)&sto_npwx,
#endif
                        bsize_psi,
                        j1l,
                        j2l,
                        velop,
                        ik,
                        ModuleBase::IMAG_UNIT,
                        bandsinfo);

            // calculate <\psi|sqrt(1-f) exp(iHt/2)*J*exp(-iHt/2) sqrt(f)|\psi>
            cal_jmatrix(*kspsi_all,
                        f_vkspsi,
                        en.data(),
                        en_all,
                        expmtmf_fact.data(),
                        expmtf_fact.data(),
                        expmtsmfchi,
                        expmtsfchi,
                        tmphchil,
                        batch_vchi,
                        batch_vhchi,
#ifdef __MPI
                        chi_all,
                        hchi_all,
                        (void*)&ks_fact,
                        (void*)&sto_npwx,
#endif
                        bsize_psi,
                        j1r,
                        j2r,
                        velop,
                        ik,
                        ModuleBase::ONE,
                        bandsinfo);

            // prepare for parallel
            int num_per = parajmat.num_per;
            int st_per = parajmat.start;
            // Re(i<psi|sqrt(f)j(1-f) exp(iHt)|psi><psi|j exp(-iHt)\sqrt(f)|psi>)
            // Im(l_ij*r_ji) = Re(-il_ij * r_ji) = Re( ((il)^+_ji)^* * r_ji)=Re(((il)^+_i)^* * r^+_i)
            // ddot_real = real(A_i^* * B_i)
            ModuleBase::timer::tick("Sto_EleCond", "ddot_real");
            ct11[it] += static_cast<double>(
                ModuleBase::GlobalFunc::ddot_real(num_per, j1l.data() + st_per, j1r.data() + st_per, false)
                * p_kv->wk[ik] / 2.0);
            double tmp12 = static_cast<double>(
                ModuleBase::GlobalFunc::ddot_real(num_per, j1l.data() + st_per, j2r.data() + st_per, false));

            double tmp21 = static_cast<double>(
                ModuleBase::GlobalFunc::ddot_real(num_per, j2l.data() + st_per, j1r.data() + st_per, false));

            ct12[it] -= 0.5 * (tmp12 + tmp21) * p_kv->wk[ik] / 2.0;

            ct22[it] += static_cast<double>(
                ModuleBase::GlobalFunc::ddot_real(num_per, j2l.data() + st_per, j2r.data() + st_per, false)
                * p_kv->wk[ik] / 2.0);

            ModuleBase::timer::tick("Sto_EleCond", "ddot_real");
        }
        std::cout << std::endl;
    } // ik loop
    ModuleBase::timer::tick("Sto_EleCond", "kloop");
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ct11.data(), nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct12.data(), nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct22.data(), nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    //------------------------------------------------------------------
    //                    Output
    //------------------------------------------------------------------
    if (GlobalV::MY_RANK == 0)
    {
        calcondw(nt, dt, smear_type, fwhmin, wcut, dw_in, ct11.data(), ct12.data(), ct22.data());
    }
    ModuleBase::timer::tick("Sto_EleCond", "sKG");
}

namespace GlobalTemp
{
const ModuleBase::matrix* veff;
}
