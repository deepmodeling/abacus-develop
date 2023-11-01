#include <chrono>

#include "esolver_sdft_pw.h"
#include "module_base/complexmatrix.h"
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw_sdft.h"

#define TWOSQRT2LN2 2.354820045030949 // FWHM = 2sqrt(2ln2) * \sigma
#define FACTOR 1.839939223835727e7
namespace ModuleESolver
{
struct parallel_distribution
{
    parallel_distribution(const int& num_all, const int& np, const int myrank)
    {
        int num_per = num_all / np;
        int st_per = num_per * myrank;
        int re = num_all % np;
        if (myrank < re)
        {
            ++num_per;
            st_per += myrank;
        }
        else
        {
            st_per += re;
        }
        this->start = st_per;
        this->num_per = num_per;
    }
    int start;
    int num_per;
};

#ifdef __MPI
struct info_gatherv
{
    info_gatherv(const int& ngroup_per, const int& np, const int& num_in_group, MPI_Comm comm_world)
    {
        nrecv = new int[np];
        displs = new int[np];
        MPI_Allgather(&ngroup_per, 1, MPI_INT, nrecv, 1, MPI_INT, comm_world);
        displs[0] = 0;
        for (int i = 1; i < np; ++i)
        {
            displs[i] = displs[i - 1] + nrecv[i - 1];
        }
        for (int i = 0; i < np; ++i)
        {
            nrecv[i] *= num_in_group;
            displs[i] *= num_in_group;
        }
    }
    ~info_gatherv()
    {
        delete[] nrecv;
        delete[] displs;
    }
    int* nrecv = nullptr;
    int* displs = nullptr;
};
#endif

psi::Psi<std::complex<double>>* gatherchi(psi::Psi<std::complex<double>>& chi,
                                          psi::Psi<std::complex<double>>& chi_all,
                                          const int& npwx,
                                          int* nrecv_sto,
                                          int* displs_sto,
                                          const int perbands_sto)
{
    psi::Psi<std::complex<double>>* p_chi;
    p_chi = &chi;
#ifdef __MPI
    if (GlobalV::NSTOGROUP > 1)
    {
        p_chi = &chi_all;
        ModuleBase::timer::tick("sKG", "bands_gather");
        MPI_Allgatherv(chi.get_pointer(),
                       perbands_sto * npwx,
                       MPI_DOUBLE_COMPLEX,
                       chi_all.get_pointer(),
                       nrecv_sto,
                       displs_sto,
                       MPI_DOUBLE_COMPLEX,
                       PARAPW_WORLD);
        ModuleBase::timer::tick("sKG", "bands_gather");
    }
#endif
    return p_chi;
}

void ESolver_SDFT_PW::check_che(const int nche_in)
{
    //------------------------------
    //      Convergence test
    //------------------------------
    bool change = false;
    const int nk = kv.nks;
    ModuleBase::Chebyshev<double> chetest(nche_in);
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    Stochastic_hchi& stohchi = stoiter.stohchi;
    int ntest0 = 5;
    stohchi.Emax = pw_wfc->gk_ecut * pw_wfc->tpiba2;
    if (GlobalV::NBANDS > 0)
    {
        double tmpemin = 1e10;
        for (int ik = 0; ik < nk; ++ik)
        {
            tmpemin = std::min(tmpemin, this->pelec->ekb(ik, GlobalV::NBANDS - 1));
        }
        stohchi.Emin = tmpemin;
    }
    else
    {
        stohchi.Emin = 0;
    }
    for (int ik = 0; ik < nk; ++ik)
    {
        this->p_hamilt->updateHk(ik);
        stoiter.stohchi.current_ik = ik;
        const int npw = kv.ngk[ik];
        std::complex<double>* pchi = nullptr;
        int ntest = std::min(ntest0, stowf.nchip[ik]);
        for (int i = 0; i < ntest; ++i)
        {
            if (INPUT.nbands_sto == 0)
            {
                pchi = new std::complex<double>[npw];
                for (int ig = 0; ig < npw; ++ig)
                {
                    double rr = std::rand() / double(RAND_MAX);
                    double arg = std::rand() / double(RAND_MAX);
                    pchi[ig] = std::complex<double>(rr * cos(arg), rr * sin(arg));
                }
            }
            else if (GlobalV::NBANDS > 0)
            {
                pchi = &stowf.chiortho[0](ik, i, 0);
            }
            else
            {
                pchi = &stowf.chi0[0](ik, i, 0);
            }
            while (1)
            {
                bool converge;
                converge = chetest.checkconverge(&stohchi,
                                                 &Stochastic_hchi::hchi_norm,
                                                 pchi,
                                                 npw,
                                                 stohchi.Emax,
                                                 stohchi.Emin,
                                                 2.0);

                if (!converge)
                {
                    change = true;
                }
                else
                {
                    break;
                }
            }
            if (INPUT.nbands_sto == 0)
            {
                delete[] pchi;
            }
        }

        if (ik == nk - 1)
        {
#ifdef __MPI
            MPI_Allreduce(MPI_IN_PLACE, &stohchi.Emax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &stohchi.Emin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
            stoiter.stofunc.Emax = stohchi.Emax;
            stoiter.stofunc.Emin = stohchi.Emin;
            GlobalV::ofs_running << "New Emax " << stohchi.Emax << " ; new Emin " << stohchi.Emin << std::endl;
            change = false;
        }
    }
}

int ESolver_SDFT_PW::set_cond_nche(const double dt, const int nbatch, const double cond_thr)
{
    int nche_guess = 1000;
    ModuleBase::Chebyshev<double> chet(nche_guess);
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    const double mu = this->pelec->eferm.ef;
    stoiter.stofunc.mu = mu;
    stoiter.stofunc.t = 0.5 * dt * nbatch;
    chet.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::n_sin);

    int nche;
    bool find = false;
    std::ofstream cheofs("Chebycoef");
    for (int i = 1; i < nche_guess; ++i)
    {
        double error = std::abs(chet.coef_complex[i] / chet.coef_complex[0]);
        cheofs << std::setw(5) << i << std::setw(20) << error << std::endl;
        if (!find && error < cond_thr)
        {
            nche = i + 1;
            std::cout << "set N order of Chebyshev for KG as " << nche << std::endl;
            find = true;
        }
    }
    cheofs.close();

    if (!find)
    {
        ModuleBase::WARNING_QUIT("ESolver_SDFT_PW", "N order of Chebyshev for KG will be larger than 1000!");
    }

    return nche;
}

void ESolver_SDFT_PW::cal_jmatrix(const psi::Psi<std::complex<double>>& kspsi,
                                  const psi::Psi<std::complex<double>>& vkspsi,
                                  const std::complex<double>* psi_all,
                                  const double* en,
                                  const double* en_all,
                                  std::complex<double>* leftfact,
                                  std::complex<double>* rightfact,
                                  const psi::Psi<std::complex<double>>& leftchi,
                                  psi::Psi<std::complex<double>>& rightchi,
                                  psi::Psi<std::complex<double>>& left_hchi,
                                  psi::Psi<std::complex<double>>& right_hchi,
                                  psi::Psi<std::complex<double>>& batch_vchi,
                                  psi::Psi<std::complex<double>>& batch_vhchi,
#ifdef __MPI
                                  psi::Psi<std::complex<double>>& chi_all,
                                  psi::Psi<std::complex<double>>& hchi_all,
                                  void* gatherinfo_ks,
                                  void* gatherinfo_sto,
#endif
                                  const int& bsize_psi,
                                  ModuleBase::ComplexMatrix& j1,
                                  ModuleBase::ComplexMatrix& j2,
                                  hamilt::Velocity& velop,
                                  const int& ik,
                                  const std::complex<double>& factor,
                                  const int bandinfo[6])
{
    ModuleBase::timer::tick(this->classname, "cal_jmatrix");
    const char transa = 'C';
    const char transb = 'N';
    const std::complex<double>& conjfactor = std::conj(factor);
    const double mu = this->pelec->eferm.ef;
    const int npwx = wf.npwx;
    const int npw = kv.ngk[ik];
    const int ndim = 3;
    const int perbands_ks = bandinfo[0];
    const int perbands_sto = bandinfo[1];
    const int perbands = bandinfo[2];
    const int allbands_ks = bandinfo[3];
    const int allbands_sto = bandinfo[4];
    const int allbands = bandinfo[5];
    const int dim_jmatrix = perbands_ks * allbands_sto + perbands_sto * allbands;
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;

    ModuleBase::ComplexMatrix tmpj(ndim, allbands_sto * perbands_sto);

    stoiter.stohchi.hchi(leftchi.get_pointer(), left_hchi.get_pointer(), perbands_sto);
    stoiter.stohchi.hchi(rightchi.get_pointer(), right_hchi.get_pointer(), perbands_sto);

    psi::Psi<std::complex<double>>* rightchi_all = &rightchi;
    psi::Psi<std::complex<double>>* righthchi_all = &right_hchi;
    std::complex<double>*tmprightf_all = nullptr, *rightf_all = rightfact;
    if (GlobalV::NSTOGROUP > 1 && rightfact != nullptr)
    {
        tmprightf_all = new std::complex<double>[allbands_ks];
        rightf_all = tmprightf_all;
    }
#ifdef __MPI
    info_gatherv* ks_fact = static_cast<info_gatherv*>(gatherinfo_ks);
    info_gatherv* sto_npwx = static_cast<info_gatherv*>(gatherinfo_sto);
    rightchi_all = gatherchi(rightchi, chi_all, npwx, sto_npwx->nrecv, sto_npwx->displs, perbands_sto);
    righthchi_all = gatherchi(right_hchi, hchi_all, npwx, sto_npwx->nrecv, sto_npwx->displs, perbands_sto);
    if (GlobalV::NSTOGROUP > 1 && rightfact != nullptr)
    {
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

    // 1. (<\psi|J|\chi>)^T
    // (allbands_sto, perbands_ks)
    if (perbands_ks > 0)
    {
        for (int id = 0; id < ndim; ++id)
        {
            const int idnb = id * perbands_ks;
            const int jbais = 0;
            //(<\psi|v|\chi>)^T
            zgemm_(&transa,
                   &transb,
                   &allbands_sto,
                   &perbands_ks,
                   &npw,
                   &conjfactor,
                   rightchi_all->get_pointer(),
                   &npwx,
                   &vkspsi(idnb, 0),
                   &npwx,
                   &ModuleBase::ZERO,
                   &j1(id, jbais),
                   &allbands_sto);

            //(<\psi|Hv|\chi>)^T
            // for(int i = 0 ; i < perbands_ks ; ++i)
            // {
            //     double* evmat = &j2(id, jbais + i * allbands_sto);
            //     double* vmat = &j1(id, jbais + i * allbands_sto);
            //     double ei = en[i];
            //     for(int j = 0 ; j < allbands_sto ; ++j)
            //     {
            //         evmat[j] = vmat[j] * ei;
            //     }
            // }

            //(<\psi|vH|\chi>)^T
            zgemm_(&transa,
                   &transb,
                   &allbands_sto,
                   &perbands_ks,
                   &npw,
                   &conjfactor,
                   righthchi_all->get_pointer(),
                   &npwx,
                   &vkspsi(idnb, 0),
                   &npwx,
                   &ModuleBase::ZERO,
                   &j2(id, jbais),
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
        // v|H\chi>
        velop.act(&leftchi, tmpnb, &left_hchi(0, startnb, 0), batch_vhchi.get_pointer());

        // 2. <\chi|J|\psi>
        // (perbands_sto, allbands_ks)
        if (allbands_ks > 0)
        {
            for (int id = 0; id < ndim; ++id)
            {
                const int idnb = id * tmpnb;
                const int jbais = perbands_ks * allbands_sto + startnb;
                //<\chi|v|\psi>
                zgemm_(&transa,
                       &transb,
                       &tmpnb,
                       &allbands_ks,
                       &npw,
                       &factor,
                       &batch_vchi(idnb, 0),
                       &npwx,
                       psi_all,
                       &npwx,
                       &ModuleBase::ZERO,
                       &j1(id, jbais),
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
                zgemm_(&transa,
                       &transb,
                       &tmpnb,
                       &allbands_ks,
                       &npw,
                       &factor,
                       &batch_vhchi(idnb, 0),
                       &npwx,
                       psi_all,
                       &npwx,
                       &ModuleBase::ZERO,
                       &j2(id, jbais),
                       &perbands_sto);
            }
        }

        // 3. <\chi|J|\chi>
        // (perbands_sto, allbands_sto)
        for (int id = 0; id < ndim; ++id)
        {
            const int idnb = id * tmpnb;
            const int jbais = perbands_ks * allbands_sto + perbands_sto * allbands_ks + startnb;
            //<\chi|v|\chi>
            zgemm_(&transa,
                   &transb,
                   &tmpnb,
                   &allbands_sto,
                   &npw,
                   &factor,
                   &batch_vchi(idnb, 0),
                   &npwx,
                   rightchi_all->get_pointer(),
                   &npwx,
                   &ModuleBase::ZERO,
                   &j1(id, jbais),
                   &perbands_sto);

            //<\chi|Hv|\chi>
            zgemm_(&transa,
                   &transb,
                   &tmpnb,
                   &allbands_sto,
                   &npw,
                   &factor,
                   &batch_vhchi(idnb, 0),
                   &npwx,
                   rightchi_all->get_pointer(),
                   &npwx,
                   &ModuleBase::ZERO,
                   &j2(id, jbais),
                   &perbands_sto);

            //<\chi|vH|\chi>
            zgemm_(&transa,
                   &transb,
                   &tmpnb,
                   &allbands_sto,
                   &npw,
                   &factor,
                   &batch_vchi(idnb, 0),
                   &npwx,
                   righthchi_all->get_pointer(),
                   &npwx,
                   &ModuleBase::ZERO,
                   &tmpj(id, startnb),
                   &perbands_sto);
        }

        remain -= tmpnb;
        startnb += tmpnb;
        if (remain == 0)
            break;
    }

    for (int id = 0; id < ndim; ++id)
    {
        for (int i = 0; i < perbands_ks; ++i)
        {
            double ei = en[i];
            const int jst = i * allbands_sto;
            const int jend = jst + allbands_sto;
            if (leftfact == nullptr)
            {
                for (int j = jst; j < jend; ++j)
                {
                    j2(id, j) = 0.5 * j2(id, j) + (0.5 * ei - mu) * j1(id, j);
                }
            }
            else
            {
                const std::complex<double> factor = leftfact[i];
                for (int j = jst; j < jend; ++j)
                {
                    j2(id, j) = factor * (0.5 * j2(id, j) + (0.5 * ei - mu) * j1(id, j));
                    j1(id, j) *= factor;
                }
            }
        }

        for (int i = 0; i < allbands_ks; ++i)
        {
            double ei = en_all[i];

            const int jst = perbands_ks * allbands_sto + i * perbands_sto;
            const int jend = jst + perbands_sto;
            if (rightfact == nullptr)
            {
                for (int j = jst; j < jend; ++j)
                {
                    j2(id, j) = 0.5 * j2(id, j) + (0.5 * ei - mu) * j1(id, j);
                }
            }
            else
            {
                const std::complex<double> factor = rightf_all[i];
                for (int j = jst; j < jend; ++j)
                {
                    j2(id, j) = factor * (0.5 * j2(id, j) + (0.5 * ei - mu) * j1(id, j));
                    j1(id, j) *= factor;
                }
            }
        }

        const int st = perbands_ks * allbands_sto + perbands_sto * allbands_ks;
        for (int j = st; j < dim_jmatrix; ++j)
        {
            j2(id, j) = 0.5 * (j2(id, j) + tmpj(id, j - st)) - mu * j1(id, j);
        }
    }

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, j1.c, ndim * dim_jmatrix, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, j2.c, ndim * dim_jmatrix, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
    delete[] tmprightf_all;

    ModuleBase::timer::tick(this->classname, "cal_jmatrix");

    return;
}

void ESolver_SDFT_PW::sKG(const int nche_KG,
                          const double fwhmin,
                          const double wcut,
                          const double dw_in,
                          const double dt_in,
                          const int nbatch,
                          const int npart_sto)
{
    ModuleBase::TITLE(this->classname, "sKG");
    ModuleBase::timer::tick(this->classname, "sKG");
    std::cout << "Calculating conductivity...." << std::endl;
    // if (GlobalV::NSTOGROUP > 1)
    // {
    //     ModuleBase::WARNING_QUIT("ESolver_SDFT_PW", "sKG is not supported in parallel!");
    // }

    //------------------------------------------------------------------
    //                    Init
    //------------------------------------------------------------------
    // Parameters
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    double dt = dt_in;                               // unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    const double expfactor = 18.42;                  // exp(-18.42) = 1e-8
    int nt = ceil(sqrt(2 * expfactor) / sigma / dt); // set nt empirically
    std::cout << "nw: " << nw << " ; dw: " << dw * ModuleBase::Ry_to_eV << " eV" << std::endl;
    std::cout << "nt: " << nt << " ; dt: " << dt << " a.u.(ry^-1)" << std::endl;
    assert(nw >= 1);
    assert(nt >= 1);
    const int ndim = 3;
    const int nk = kv.nks;
    const int npwx = wf.npwx;
    const double tpiba = GlobalC::ucell.tpiba;
    psi::Psi<std::complex<double>>* stopsi;
    if (GlobalV::NBANDS > 0)
    {
        stopsi = stowf.chiortho;
        stowf.shchi->resize(
            1,
            1,
            1); // clean memories //Note shchi is different from \sqrt(fH_here)|chi>, since veffs are different
        stowf.chi0->resize(1, 1, 1); // clean memories
    }
    else
    {
        stopsi = stowf.chi0;
        stowf.shchi->resize(1, 1, 1); // clean memories
    }
    const double dEcut = (wcut + fwhmin) / ModuleBase::Ry_to_eV;

    // response funtion
    double* ct11 = new double[nt];
    double* ct12 = new double[nt];
    double* ct22 = new double[nt];
    ModuleBase::GlobalFunc::ZEROS(ct11, nt);
    ModuleBase::GlobalFunc::ZEROS(ct12, nt);
    ModuleBase::GlobalFunc::ZEROS(ct22, nt);

    // Init Chebyshev
    ModuleBase::Chebyshev<double> che(this->nche_sto);
    ModuleBase::Chebyshev<double> chet(nche_KG);
    ModuleBase::Chebyshev<double> chemt(nche_KG);
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    Stochastic_hchi& stohchi = stoiter.stohchi;

    //------------------------------------------------------------------
    //                    Calculate
    //------------------------------------------------------------------

    // Prepare Chebyshev coefficients for exp(-i H/\hbar t)
    const double mu = this->pelec->eferm.ef;
    stoiter.stofunc.mu = mu;
    stoiter.stofunc.t = 0.5 * dt * nbatch;
    chet.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::nsin);
    chemt.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::n_sin);
    std::complex<double>*batchcoef = nullptr, *batchmcoef = nullptr;
    if (nbatch > 1)
    {
        batchcoef = new std::complex<double>[nche_KG * nbatch];
        std::complex<double>* tmpcoef = batchcoef + (nbatch - 1) * nche_KG;
        batchmcoef = new std::complex<double>[nche_KG * nbatch];
        std::complex<double>* tmpmcoef = batchmcoef + (nbatch - 1) * nche_KG;
        for (int i = 0; i < nche_KG; ++i)
        {
            tmpcoef[i] = chet.coef_complex[i];
            tmpmcoef[i] = chemt.coef_complex[i];
        }
        for (int ib = 0; ib < nbatch - 1; ++ib)
        {
            tmpcoef = batchcoef + ib * nche_KG;
            tmpmcoef = batchmcoef + ib * nche_KG;
            stoiter.stofunc.t = 0.5 * dt * (ib + 1);
            chet.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::nsin);
            chemt.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::n_sin);
            for (int i = 0; i < nche_KG; ++i)
            {
                tmpcoef[i] = chet.coef_complex[i];
                tmpmcoef[i] = chemt.coef_complex[i];
            }
        }
        stoiter.stofunc.t = 0.5 * dt * nbatch;
    }

    // ik loop
    ModuleBase::timer::tick(this->classname, "kloop");
    hamilt::Velocity velop(pw_wfc, kv.isk.data(), &GlobalC::ppcell, &GlobalC::ucell, INPUT.cond_nonlocal);
    for (int ik = 0; ik < nk; ++ik)
    {
        velop.init(ik);
        stopsi->fix_k(ik);
        psi->fix_k(ik);
        if (nk > 1)
        {
            this->p_hamilt->updateHk(ik);
        }
        stoiter.stohchi.current_ik = ik;
        const int npw = kv.ngk[ik];

        // get allbands_ks
        int cutib0 = 0;
        if (GlobalV::NBANDS > 1)
        {
            double Emax_KS = this->pelec->ekb(ik, GlobalV::NBANDS - 1);
            for (cutib0 = GlobalV::NBANDS - 2; cutib0 >= 0; --cutib0)
            {
                if (Emax_KS - this->pelec->ekb(ik, cutib0) > dEcut)
                {
                    break;
                }
            }
            cutib0++;
            double Emin_KS = this->pelec->ekb(ik, cutib0);
            double dE = stoiter.stofunc.Emax - Emin_KS + wcut / ModuleBase::Ry_to_eV;
            std::cout << "Emin_KS(" << cutib0 << "): " << Emin_KS * ModuleBase::Ry_to_eV
                      << " eV; Emax: " << stoiter.stofunc.Emax * ModuleBase::Ry_to_eV
                      << " eV; Recommended dt: " << 0.25 * M_PI / dE << " a.u." << std::endl;
        }
        // Parallel for bands
        int allbands_ks = GlobalV::NBANDS - cutib0;
        parallel_distribution paraks(allbands_ks, GlobalV::NSTOGROUP, GlobalV::MY_STOGROUP);
        int perbands_ks = paraks.num_per;
        int ib0_ks = paraks.start;
        ib0_ks += GlobalV::NBANDS - allbands_ks;
        int perbands_sto = this->stowf.nchip[ik];
        int perbands = perbands_sto + perbands_ks;
        int allbands_sto = perbands_sto;
        int allbands = perbands;
#ifdef __MPI
        MPI_Allreduce(&perbands, &allbands, 1, MPI_INT, MPI_SUM, PARAPW_WORLD);
        allbands_sto = allbands - allbands_ks;
        info_gatherv ks_fact(perbands_ks, GlobalV::NSTOGROUP, 1, PARAPW_WORLD);
        info_gatherv sto_npwx(perbands_sto, GlobalV::NSTOGROUP, npwx, PARAPW_WORLD);
#endif
        const int bandsinfo[6]{perbands_ks, perbands_sto, perbands, allbands_ks, allbands_sto, allbands};
        double* en = nullptr;
        const double* en_all = &(this->pelec->ekb(ik, GlobalV::NBANDS - allbands_ks));
        const std::complex<double>* psi_all = &(psi[0](GlobalV::NBANDS - allbands_ks, 0));
        if (perbands_ks > 0)
        {
            en = new double[perbands_ks];
            for (int ib = 0; ib < perbands_ks; ++ib)
            {
                en[ib] = this->pelec->ekb(ik, ib0_ks + ib);
            }
        }

        //-----------------------------------------------------------
        //               ks conductivity
        //-----------------------------------------------------------
        if (GlobalV::MY_STOGROUP == 0 && allbands_ks > 0)
            jjcorr_ks(ik, nt, dt, dEcut, this->pelec->wg, velop, ct11, ct12, ct22);

        //-----------------------------------------------------------
        //               sto conductivity
        //-----------------------------------------------------------
        //-------------------     allocate  -------------------------
        size_t ks_memory_cost = perbands_ks * npwx * sizeof(std::complex<double>);
        psi::Psi<std::complex<double>> kspsi(1, perbands_ks, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::kspsi", ks_memory_cost);
        psi::Psi<std::complex<double>> vkspsi(1, perbands_ks * ndim, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::vkspsi", ks_memory_cost);
        std::vector<std::complex<double>> expmtmf_fact(perbands_ks), expmtf_fact(perbands_ks);

        size_t sto_memory_cost = perbands_sto * npwx * sizeof(std::complex<double>);
        psi::Psi<std::complex<double>> sfchi(1, perbands_sto, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::sfchi", sto_memory_cost);
        psi::Psi<std::complex<double>> smfchi(1, perbands_sto, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::smfchi", sto_memory_cost);
        psi::Psi<std::complex<double>> tmphchil(1, perbands_sto, npwx, kv.ngk.data());
        psi::Psi<std::complex<double>> tmphchir(1, perbands_sto, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::tmphchi", sto_memory_cost * 2);
#ifdef __MPI
        psi::Psi<std::complex<double>> chi_all, hchi_all;
        if (GlobalV::NSTOGROUP > 1)
        {
            size_t memory_all = allbands * npwx * sizeof(std::complex<double>);
            chi_all.resize(1, allbands, npwx);
            hchi_all.resize(1, allbands, npwx);
            ModuleBase::Memory::record("SDFT::chi_all", memory_all);
            ModuleBase::Memory::record("SDFT::hchi_all", memory_all);
        }
#endif

        const int dim_jmatrix = perbands_ks * allbands_sto + perbands_sto * allbands;
        parallel_distribution parajmat(ndim * dim_jmatrix, GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL);
        ModuleBase::ComplexMatrix j1l(ndim, dim_jmatrix), j2l(ndim, dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j1l", sizeof(std::complex<double>) * ndim * dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j2l", sizeof(std::complex<double>) * ndim * dim_jmatrix);
        ModuleBase::ComplexMatrix j1r(ndim, dim_jmatrix), j2r(ndim, dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j1r", sizeof(std::complex<double>) * ndim * dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j2r", sizeof(std::complex<double>) * ndim * dim_jmatrix);

        const int nbatch_psi = npart_sto;
        const int bsize_psi = ceil(double(perbands_sto) / nbatch_psi);
        psi::Psi<std::complex<double>> batch_vchi(1, bsize_psi * ndim, npwx, kv.ngk.data());
        psi::Psi<std::complex<double>> batch_vhchi(1, bsize_psi * ndim, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::batchjpsi", bsize_psi * 2 * ndim * npwx * sizeof(std::complex<double>));

        //-------------------     sqrt(f)|psi>   sqrt(1-f)|psi>   ---------------
        for (int ib = 0; ib < perbands_ks; ++ib)
        {
            for (int ig = 0; ig < npw; ++ig)
            {
                kspsi(0, ib, ig) = psi[0](ib0_ks + ib, ig);
            }
            double fi = stoiter.stofunc.fd(en[ib]);
            expmtmf_fact[ib] = 1 - fi;
            expmtf_fact[ib] = fi;
        }
        che.calcoef_real(&stoiter.stofunc, &Sto_Func<double>::nroot_fd);
        che.calfinalvec_real(&stohchi,
                             &Stochastic_hchi::hchi_norm,
                             stopsi->get_pointer(),
                             sfchi.get_pointer(),
                             npw,
                             npwx,
                             perbands_sto);
        che.calcoef_real(&stoiter.stofunc, &Sto_Func<double>::nroot_mfd);
        che.calfinalvec_real(&stohchi,
                             &Stochastic_hchi::hchi_norm,
                             stopsi->get_pointer(),
                             smfchi.get_pointer(),
                             npw,
                             npwx,
                             perbands_sto);

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
            poly_exptsfchi.resize(nche_KG, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_exptsfchi",
                                       sizeof(std::complex<double>) * nche_KG * perbands_sto * npwx);

            poly_exptsmfchi.resize(nche_KG, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_exptsmfchi",
                                       sizeof(std::complex<double>) * nche_KG * perbands_sto * npwx);

            poly_expmtsfchi.resize(nche_KG, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_expmtsfchi",
                                       sizeof(std::complex<double>) * nche_KG * perbands_sto * npwx);

            poly_expmtsmfchi.resize(nche_KG, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_expmtsmfchi",
                                       sizeof(std::complex<double>) * nche_KG * perbands_sto * npwx);
        }
        // v|\psi>
        velop.act(&kspsi, perbands_ks, kspsi.get_pointer(), vkspsi.get_pointer());

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
                std::cout << "nt: ";
            }
            if ((it - 1) % print_step == 0 && it > 1)
                std::cout << it - 1 << " ";

            // time evolution exp(-iHt)|\psi_ks>
            // KS
            ModuleBase::timer::tick(this->classname, "evolution");
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
                chemt.calfinalvec_complex(&stohchi,
                                          &Stochastic_hchi::hchi_norm,
                                          expmtsfchi.get_pointer(),
                                          expmtsfchi.get_pointer(),
                                          npw,
                                          npwx,
                                          perbands_sto);
                chemt.calfinalvec_complex(&stohchi,
                                          &Stochastic_hchi::hchi_norm,
                                          expmtsmfchi.get_pointer(),
                                          expmtsmfchi.get_pointer(),
                                          npw,
                                          npwx,
                                          perbands_sto);
                chet.calfinalvec_complex(&stohchi,
                                         &Stochastic_hchi::hchi_norm,
                                         exptsfchi.get_pointer(),
                                         exptsfchi.get_pointer(),
                                         npw,
                                         npwx,
                                         perbands_sto);
                chet.calfinalvec_complex(&stohchi,
                                         &Stochastic_hchi::hchi_norm,
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
                    chet.calpolyvec_complex(&stohchi,
                                            &Stochastic_hchi::hchi_norm,
                                            stoexptsfchi,
                                            tmppolyexptsfchi,
                                            npw,
                                            npwx,
                                            perbands_sto);
                    chet.calpolyvec_complex(&stohchi,
                                            &Stochastic_hchi::hchi_norm,
                                            stoexptsmfchi,
                                            tmppolyexptsmfchi,
                                            npw,
                                            npwx,
                                            perbands_sto);
                    chemt.calpolyvec_complex(&stohchi,
                                             &Stochastic_hchi::hchi_norm,
                                             stoexpmtsfchi,
                                             tmppolyexpmtsfchi,
                                             npw,
                                             npwx,
                                             perbands_sto);
                    chemt.calpolyvec_complex(&stohchi,
                                             &Stochastic_hchi::hchi_norm,
                                             stoexpmtsmfchi,
                                             tmppolyexpmtsmfchi,
                                             npw,
                                             npwx,
                                             perbands_sto);
                }

                std::complex<double>* tmpcoef = batchcoef + (it - 1) % nbatch * nche_KG;
                std::complex<double>* tmpmcoef = batchmcoef + (it - 1) % nbatch * nche_KG;
                const char transa = 'N';
                const int LDA = perbands_sto * npwx;
                const int M = perbands_sto * npwx;
                const int N = nche_KG;
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
            ModuleBase::timer::tick(this->classname, "evolution");

            // calculate i<\psi|sqrt(f) exp(-iHt/2)*J*exp(iHt/2) sqrt(1-f)|\psi>^+
            //         = i<\psi|sqrt(1-f) exp(-iHt/2)*J*exp(iHt/2) sqrt(f)|\psi>
            cal_jmatrix(kspsi,
                        vkspsi,
                        psi_all,
                        en,
                        en_all,
                        nullptr,
                        nullptr,
                        exptsmfchi,
                        exptsfchi,
                        tmphchil,
                        tmphchir,
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
            cal_jmatrix(kspsi,
                        vkspsi,
                        psi_all,
                        en,
                        en_all,
                        expmtmf_fact.data(),
                        expmtf_fact.data(),
                        expmtsmfchi,
                        expmtsfchi,
                        tmphchil,
                        tmphchir,
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
            // Im(l_ij*r_ji)=Re(i l^*_ij*r^+_ij)=Re(i l^*_i*r^+_i)
            // ddot_real = real(A^*_i * B_i)
            ModuleBase::timer::tick(this->classname, "ddot_real");
            ct11[it]
                += ModuleBase::GlobalFunc::ddot_real(num_per, j1l.c + st_per, j1r.c + st_per, false) * kv.wk[ik] / 2.0;
            double tmp12 = ModuleBase::GlobalFunc::ddot_real(num_per, j1l.c + st_per, j2r.c + st_per, false);
            double tmp21 = ModuleBase::GlobalFunc::ddot_real(num_per, j2l.c + st_per, j1r.c + st_per, false);
            ct12[it] -= 0.5 * (tmp12 + tmp21) * kv.wk[ik] / 2.0;
            ct22[it]
                += ModuleBase::GlobalFunc::ddot_real(num_per, j2l.c + st_per, j2r.c + st_per, false) * kv.wk[ik] / 2.0;
            ModuleBase::timer::tick(this->classname, "ddot_real");
        }
        std::cout << std::endl;
        delete[] en;
    } // ik loop
    ModuleBase::timer::tick(this->classname, "kloop");
    delete[] batchcoef;
    delete[] batchmcoef;

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ct11, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct12, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct22, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    //------------------------------------------------------------------
    //                    Output
    //------------------------------------------------------------------
    if (GlobalV::MY_RANK == 0)
    {
        calcondw(nt, dt, fwhmin, wcut, dw_in, ct11, ct12, ct22);
    }
    delete[] ct11;
    delete[] ct12;
    delete[] ct22;
    ModuleBase::timer::tick(this->classname, "sKG");
}

void ESolver_SDFT_PW::caldos(const int nche_dos,
                             const double sigmain,
                             const double emin,
                             const double emax,
                             const double de,
                             const int npart)
{
    ModuleBase::TITLE(this->classname, "caldos");
    ModuleBase::timer::tick(this->classname, "caldos");
    std::cout << "=========================" << std::endl;
    std::cout << "###Calculating Dos....###" << std::endl;
    std::cout << "=========================" << std::endl;
    ModuleBase::Chebyshev<double> che(nche_dos);
    const int nk = kv.nks;
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    Stochastic_hchi& stohchi = stoiter.stohchi;
    const int npwx = wf.npwx;

    double* spolyv = nullptr;
    std::complex<double>* allorderchi = nullptr;
    if (stoiter.method == 1)
    {
        spolyv = new double[nche_dos];
        ModuleBase::GlobalFunc::ZEROS(spolyv, nche_dos);
    }
    else
    {
        spolyv = new double[nche_dos * nche_dos];
        ModuleBase::GlobalFunc::ZEROS(spolyv, nche_dos * nche_dos);
        int nchip_new = ceil((double)this->stowf.nchip_max / npart);
        allorderchi = new std::complex<double>[nchip_new * npwx * nche_dos];
    }
    ModuleBase::timer::tick(this->classname, "Tracepoly");
    std::cout << "1. TracepolyA:" << std::endl;
    for (int ik = 0; ik < nk; ik++)
    {
        std::cout << "ik: " << ik + 1 << std::endl;
        if (nk > 1)
        {
            this->p_hamilt->updateHk(ik);
        }
        stohchi.current_ik = ik;
        const int npw = kv.ngk[ik];
        const int nchipk = this->stowf.nchip[ik];

        std::complex<double>* pchi;
        if (GlobalV::NBANDS > 0)
        {
            stowf.chiortho->fix_k(ik);
            pchi = stowf.chiortho->get_pointer();
        }
        else
        {
            stowf.chi0->fix_k(ik);
            pchi = stowf.chi0->get_pointer();
        }
        if (stoiter.method == 1)
        {
            che.tracepolyA(&stohchi, &Stochastic_hchi::hchi_norm, pchi, npw, npwx, nchipk);
            for (int i = 0; i < nche_dos; ++i)
            {
                spolyv[i] += che.polytrace[i] * kv.wk[ik] / 2;
            }
        }
        else
        {
            int N = nche_dos;
            double kweight = kv.wk[ik] / 2;
            char trans = 'T';
            char normal = 'N';
            double one = 1;
            for (int ipart = 0; ipart < npart; ++ipart)
            {
                int nchipk_new = nchipk / npart;
                int start_nchipk = ipart * nchipk_new + nchipk % npart;
                if (ipart < nchipk % npart)
                {
                    nchipk_new++;
                    start_nchipk = ipart * nchipk_new;
                }
                ModuleBase::GlobalFunc::ZEROS(allorderchi, nchipk_new * npwx * nche_dos);
                std::complex<double>* tmpchi = pchi + start_nchipk * npwx;
                che.calpolyvec_complex(&stohchi,
                                       &Stochastic_hchi::hchi_norm,
                                       tmpchi,
                                       allorderchi,
                                       npw,
                                       npwx,
                                       nchipk_new);
                double* vec_all = (double*)allorderchi;
                int LDA = npwx * nchipk_new * 2;
                int M = npwx * nchipk_new * 2;
                dgemm_(&trans, &normal, &N, &N, &M, &kweight, vec_all, &LDA, vec_all, &LDA, &one, spolyv, &N);
            }
        }
    }
    if (stoiter.method == 2)
        delete[] allorderchi;

    std::ofstream ofsdos;
    int ndos = int((emax - emin) / de) + 1;
    stoiter.stofunc.sigma = sigmain / ModuleBase::Ry_to_eV;
    ModuleBase::timer::tick(this->classname, "Tracepoly");

    std::cout << "2. Dos:" << std::endl;
    ModuleBase::timer::tick(this->classname, "DOS Loop");
    int n10 = ndos / 10;
    int percent = 10;
    double* sto_dos = new double[ndos];
    double* ks_dos = new double[ndos];
    double* error = new double[ndos];
    for (int ie = 0; ie < ndos; ++ie)
    {
        double tmpks = 0;
        double tmpsto = 0;
        stoiter.stofunc.targ_e = (emin + ie * de) / ModuleBase::Ry_to_eV;
        if (stoiter.method == 1)
        {
            che.calcoef_real(&stoiter.stofunc, &Sto_Func<double>::ngauss);
            tmpsto = BlasConnector::dot(nche_dos, che.coef_real, 1, spolyv, 1);
        }
        else
        {
            che.calcoef_real(&stoiter.stofunc, &Sto_Func<double>::nroot_gauss);
            tmpsto = stoiter.vTMv(che.coef_real, spolyv, nche_dos);
        }
        if (GlobalV::NBANDS > 0)
        {
            for (int ik = 0; ik < nk; ++ik)
            {
                double* en = &(this->pelec->ekb(ik, 0));
                for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
                {
                    tmpks += stoiter.stofunc.gauss(en[ib]) * kv.wk[ik] / 2;
                }
            }
        }
        tmpks /= GlobalV::NPROC_IN_POOL;

        double tmperror = 0;
        if (stoiter.method == 1)
        {
            tmperror = che.coef_real[nche_dos - 1] * spolyv[nche_dos - 1];
        }
        else
        {
            const int norder = nche_dos;
            double last_coef = che.coef_real[norder - 1];
            double last_spolyv = spolyv[norder * norder - 1];
            tmperror = last_coef
                       * (BlasConnector::dot(norder, che.coef_real, 1, spolyv + norder * (norder - 1), 1)
                          + BlasConnector::dot(norder, che.coef_real, 1, spolyv + norder - 1, norder)
                          - last_coef * last_spolyv);
        }

        if (ie % n10 == n10 - 1)
        {
            std::cout << percent << "%"
                      << " ";
            percent += 10;
        }
        sto_dos[ie] = tmpsto;
        ks_dos[ie] = tmpks;
        error[ie] = tmperror;
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ks_dos, ndos, MPI_DOUBLE, MPI_SUM, STO_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, sto_dos, ndos, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, error, ndos, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (GlobalV::MY_RANK == 0)
    {
        std::string dosfile = GlobalV::global_out_dir + "DOS1_smearing.dat";
        ofsdos.open(dosfile.c_str());
        double maxerror = 0;
        double sum = 0;
        ofsdos << std::setw(8) << "## E(eV) " << std::setw(20) << "dos(eV^-1)" << std::setw(20) << "sum"
               << std::setw(20) << "Error(eV^-1)" << std::endl;
        for (int ie = 0; ie < ndos; ++ie)
        {
            double tmperror = 2.0 * std::abs(error[ie]);
            if (maxerror < tmperror)
                maxerror = tmperror;
            double dos = 2.0 * (ks_dos[ie] + sto_dos[ie]) / ModuleBase::Ry_to_eV;
            sum += dos;
            ofsdos << std::setw(8) << emin + ie * de << std::setw(20) << dos << std::setw(20) << sum * de
                   << std::setw(20) << tmperror << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Finish DOS" << std::endl;
        std::cout << std::scientific << "DOS max absolute Chebyshev Error: " << maxerror << std::endl;
        ofsdos.close();
    }
    delete[] sto_dos;
    delete[] ks_dos;
    delete[] error;
    delete[] spolyv;
    ModuleBase::timer::tick(this->classname, "DOS Loop");
    ModuleBase::timer::tick(this->classname, "caldos");
    return;
}

} // namespace ModuleESolver

namespace GlobalTemp
{

const ModuleBase::matrix* veff;

}