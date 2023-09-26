#include <chrono>

#include "./esolver_sdft_pw.h"
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
void parallelks(const int& allbands, int& perbands_ks, int& ib0_ks)
{
    perbands_ks = allbands / GlobalV::NSTOGROUP;
    ib0_ks = perbands_ks * GlobalV::MY_STOGROUP;
    if (GlobalV::MY_STOGROUP < allbands % GlobalV::NSTOGROUP)
    {
        ++perbands_ks;
        ib0_ks += GlobalV::MY_STOGROUP;
    }
    else
    {
        ib0_ks += allbands % GlobalV::NSTOGROUP;
    }
}

psi::Psi<std::complex<double>>* gatherpsi(psi::Psi<std::complex<double>>& psi,
                                          psi::Psi<std::complex<double>>& tmppsi_all,
                                          const int& npwx,
                                          int* nrecv_ks, int* displs_ks,
                                          int* nrecv_sto, int* displs_sto,
                                          const int bandinfo[6] )
{
    psi::Psi<std::complex<double>>* p_psi;
    p_psi = &psi;
#ifdef __MPI
    if (GlobalV::NSTOGROUP > 1)
    {
        const int perbands_ks = bandinfo[0];
        const int perbands_sto = bandinfo[1];
        const int allbands_ks = bandinfo[3];
        p_psi = &tmppsi_all;
        ModuleBase::timer::tick("sKG", "bands_gather");
        MPI_Allgatherv(&psi(0, 0),
                       perbands_ks * npwx,
                       MPI_DOUBLE_COMPLEX,
                       &p_psi[0](0, 0),
                       nrecv_ks,
                       displs_ks,
                       MPI_DOUBLE_COMPLEX,
                       PARAPW_WORLD);
        MPI_Allgatherv(&psi(perbands_ks, 0),
                       perbands_sto * npwx,
                       MPI_DOUBLE_COMPLEX,
                       &p_psi[0](allbands_ks, 0),
                       nrecv_sto,
                       displs_sto,
                       MPI_DOUBLE_COMPLEX,
                       PARAPW_WORLD);
         ModuleBase::timer::tick("sKG", "bands_gather");
    }
#endif
    return p_psi;

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
    stohchi.Emax = pw_wfc->gk_ecut*pw_wfc->tpiba2;
    if (GlobalV::NBANDS > 0)
    {
        double tmpemin = 1e10;
        for(int ik = 0 ; ik < nk ;++ik)
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
            std::cout<<"set N order of Chebyshev for KG as "<<nche<<std::endl;
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

void ESolver_SDFT_PW::cal_j(const psi::Psi<std::complex<double>>& psi_in,
                            psi::Psi<std::complex<double>>& j1psi,
                            psi::Psi<std::complex<double>>& j2psi,
                            hamilt::Velocity& velop,
                            const int& start_band,
                            const int& nbands,
                            const int& npw)
{
    ModuleBase::timer::tick(this->classname, "cal_j");
    const int npwx = wf.npwx;
    const int ndim = 3;
    const double mu = this->pelec->eferm.ef;
    psi::Psi<std::complex<double>> hpsi(1, nbands, npwx, kv.ngk.data());

    // H|psi>
    psi::Range psi_range(1, 0, start_band, start_band + nbands - 1);
    hamilt::Operator<std::complex<double>>::hpsi_info info_psi0(&psi_in, psi_range, hpsi.get_pointer());
    this->p_hamilt->ops->hPsi(info_psi0);

    // v|\psi>
    velop.act(&psi_in, nbands, &psi_in(0, start_band, 0), j1psi.get_pointer());

    // H|v\psi>
    psi::Range j1_range(1, 0, 0, ndim*nbands - 1);
    hamilt::Operator<std::complex<double>>::hpsi_info info_j1psi(&j1psi, j1_range, j2psi.get_pointer());
    this->p_hamilt->ops->hPsi(info_j1psi);

    //|Hv\psi> + v|H\psi>
    velop.act(this->psi, nbands, hpsi.get_pointer(), j2psi.get_pointer(), true);

    //(Hv+vH-mu)|\psi>
    for (int ib = 0; ib < nbands * 3; ++ib)
    {
        for (int ig = 0; ig < npw; ++ig)
        {
            j2psi(0, ib, ig) = j2psi(0, ib, ig) / 2.0 - mu * j1psi(0, ib, ig);
        }
    }
    
    ModuleBase::timer::tick(this->classname, "cal_j");
    return;
}

void ESolver_SDFT_PW::cal_jmatrix(const psi::Psi<std::complex<double>>& leftv,
                                  const psi::Psi<std::complex<double>>& rightv,
                                  psi::Psi<std::complex<double>>& batchj1psi,
                                  psi::Psi<std::complex<double>>& batchj2psi,
                                  ModuleBase::ComplexMatrix& j1, 
                                  ModuleBase::ComplexMatrix& j2,
                                  hamilt::Velocity& velop,
                                  const int& ik,
                                  const std::complex<double>& factor,
                                  const int bandinfo[6],
                                  const int& bsize_psi)
{
    ModuleBase::timer::tick(this->classname, "cal_jmatrix");
    const char transa = 'C';
    const char transb = 'N';
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
    
    //calculate <leftv|J|rightv>
    int remain = perbands_ks;
    int startnb = 0;
    if (perbands_ks > 0)
    {
        while (remain > 0)
        {
            int tmpnb = std::min(remain, bsize_psi);
            // calculate J|leftv>
            cal_j(leftv, batchj1psi, batchj2psi, velop, startnb, tmpnb, npw);
            // calculate i<leftv| J |rightv>
            for (int id = 0; id < ndim; ++id)
            {
                const int idnb = id * tmpnb;
                const int jbais = startnb;
                zgemm_(&transa,
                       &transb,
                       &tmpnb,
                       &allbands_sto,
                       &npw,
                       &factor,
                       &batchj1psi(idnb, 0),
                       &npwx,
                       &(rightv(allbands_ks, 0)),
                       &npwx,
                       &ModuleBase::ZERO,
                       &j1(id, jbais),
                       &perbands_ks);
                zgemm_(&transa,
                       &transb,
                       &tmpnb,
                       &allbands_sto,
                       &npw,
                       &factor,
                       &batchj2psi(idnb, 0),
                       &npwx,
                       &(rightv(allbands_ks, 0)),
                       &npwx,
                       &ModuleBase::ZERO,
                       &j2(id, jbais),
                       &perbands_ks);
            }
            remain -= tmpnb;
            startnb += tmpnb;
            if (remain == 0)
                break;
        }
    }
    remain = perbands_sto;
    startnb = 0;
    while (remain > 0)
    {
        int tmpnb = std::min(remain, bsize_psi);
        // calculate J|leftv>
        cal_j(leftv, batchj1psi, batchj2psi, velop, perbands_ks + startnb, tmpnb, npw);
       // calculate i<leftv| J |rightv>
        for (int id = 0; id < ndim; ++id)
        {
            const int idnb = id * tmpnb;
            const int jbais = perbands_ks * allbands_sto + startnb;
            zgemm_(&transa,
                   &transb,
                   &tmpnb,
                   &allbands,
                   &npw,
                   &factor,
                   &batchj1psi(idnb, 0),
                   &npwx,
                   rightv.get_pointer(),
                   &npwx,
                   &ModuleBase::ZERO,
                   &j1(id, jbais),
                   &perbands_sto);
            zgemm_(&transa,
                   &transb,
                   &tmpnb,
                   &allbands,
                   &npw,
                   &factor,
                   &batchj2psi(idnb, 0),
                   &npwx,
                   rightv.get_pointer(),
                   &npwx,
                   &ModuleBase::ZERO,
                   &j2(id, jbais),
                   &perbands_sto);
        }
        remain -= tmpnb;
        startnb += tmpnb;
        if (remain == 0)
            break;
    }

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, j1.c, ndim * dim_jmatrix, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, j2.c, ndim * dim_jmatrix, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif

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
        stowf.shchi->resize(1, 1, 1); // clean memories //Note shchi is different from \sqrt(fH_here)|chi>, since veffs are different
        stowf.chi0->resize(1, 1, 1);  // clean memories
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
    std::complex<double>* batchcoef = nullptr, * batchmcoef = nullptr;
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

        //get allbands_ks
        int cutib0 = 0;
        if( GlobalV::NBANDS > 1)
        {
            double Emax_KS = this->pelec->ekb(ik, GlobalV::NBANDS - 1);
            for(cutib0 = GlobalV::NBANDS -2 ; cutib0 >= 0 ; --cutib0)
            {
                if(Emax_KS - this->pelec->ekb(ik, cutib0) > dEcut)
                {
                    break;
                }
            }
            cutib0++;
            double Emin_KS = this->pelec->ekb(ik, cutib0);
            double dE = stoiter.stofunc.Emax - Emin_KS + wcut/ModuleBase::Ry_to_eV;
            std::cout<<"Emin_KS("<<cutib0<<"): "<<Emin_KS*ModuleBase::Ry_to_eV
                    <<" eV; Emax: "<<stoiter.stofunc.Emax*ModuleBase::Ry_to_eV
                    <<" eV; Recommended dt: "<<0.25*M_PI/dE<<" a.u."<<std::endl;
        }
        //Parallel for bands
        int allbands_ks = GlobalV::NBANDS - cutib0;
        int perbands_ks, ib0_ks;
        parallelks(allbands_ks, perbands_ks, ib0_ks);
        ib0_ks += GlobalV::NBANDS - allbands_ks;
        int perbands_sto = this->stowf.nchip[ik];
        int perbands = perbands_sto + perbands_ks;
        int allbands_sto = perbands_sto;
        int allbands = perbands;
#ifdef __MPI
        MPI_Allreduce(&perbands, &allbands, 1, MPI_INT, MPI_SUM, PARAPW_WORLD);
        allbands_sto = allbands - allbands_ks;
        const int nstogroup = GlobalV::NSTOGROUP;
        int nrecv_ks[nstogroup];
        int nrecv_sto[nstogroup];
        int displs_ks[nstogroup];
        int displs_sto[nstogroup];
        MPI_Allgather(&perbands_ks, 1, MPI_INT, nrecv_ks, 1, MPI_INT, PARAPW_WORLD);
        MPI_Allgather(&perbands_sto, 1, MPI_INT, nrecv_sto, 1, MPI_INT, PARAPW_WORLD);
        displs_ks[0] = 0;
        displs_sto[0] = 0;
        for (int i = 1; i < nstogroup; ++i)
        {
            displs_ks[i] = displs_ks[i - 1] + nrecv_ks[i - 1];
            displs_sto[i] = displs_sto[i - 1] + nrecv_sto[i - 1];
        }
        for (int i = 0; i < nstogroup; ++i)
        {
            nrecv_ks[i] *= npwx;
            nrecv_sto[i] *= npwx;
            displs_ks[i] *= npwx;
            displs_sto[i] *= npwx;
        }
#endif
        const int bandsinfo[6]{perbands_ks, perbands_sto, perbands, allbands_ks, allbands_sto, allbands};
        double* en = nullptr;
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
        size_t memory_cost = perbands * npwx * sizeof(std::complex<double>);
        psi::Psi<std::complex<double>> sfpsi(1, perbands, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::sfpsi", memory_cost);
        psi::Psi<std::complex<double>> smfpsi(1, perbands, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::smfpsi", memory_cost);
#ifdef __MPI
        psi::Psi<std::complex<double>> tmppsi_all;
        if (GlobalV::NSTOGROUP > 1)
        {
            size_t memory_all = allbands * npwx * sizeof(std::complex<double>);
            tmppsi_all.resize(1, allbands, npwx);
            ModuleBase::Memory::record("SDFT::tmppsi_all", memory_all);
        }
#endif

        const int dim_jmatrix = perbands_ks * allbands_sto + perbands_sto * allbands;
        ModuleBase::ComplexMatrix j1l(ndim, dim_jmatrix), j2l(ndim, dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j1l", sizeof(std::complex<double>) * ndim * dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j2l", sizeof(std::complex<double>) * ndim * dim_jmatrix);
        ModuleBase::ComplexMatrix j1r(ndim, dim_jmatrix), j2r(ndim, dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j1r", sizeof(std::complex<double>) * ndim * dim_jmatrix);
        ModuleBase::Memory::record("SDFT::j2r", sizeof(std::complex<double>) * ndim * dim_jmatrix);

        const int nbatch_psi = npart_sto;
        const int bsize_psi = std::min((int)ceil(perbands / nbatch_psi), std::max(perbands_sto, perbands_ks));
        psi::Psi<std::complex<double>> batchj1psi(1, bsize_psi * ndim, npwx, kv.ngk.data());
        psi::Psi<std::complex<double>> batchj2psi(1, bsize_psi * ndim, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::batchjpsi", bsize_psi * (2*ndim+1) * npwx * sizeof(std::complex<double>));

        //-------------------     sqrt(f)|psi>   sqrt(1-f)|psi>   ---------------
        for (int ib = 0; ib < perbands_ks; ++ib)
        {
            double fi = stoiter.stofunc.fd(en[ib]);
            double sfi = sqrt(fi);
            double smfi = sqrt(1-fi);
            for (int ig = 0; ig < npw; ++ig)
            {
                sfpsi(ib, ig) = psi[0](ib0_ks + ib, ig) * sfi;
                smfpsi(ib, ig) = psi[0](ib0_ks + ib, ig) * smfi;
            }
        }
        std::complex<double>* stosfpsi = sfpsi.get_pointer() + perbands_ks * npwx;
        std::complex<double>* stosmfpsi = smfpsi.get_pointer() + perbands_ks * npwx;
        che.calcoef_real(&stoiter.stofunc, &Sto_Func<double>::nroot_fd);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_norm, stopsi->get_pointer(), stosfpsi, npw, npwx, perbands_sto);
        che.calcoef_real(&stoiter.stofunc, &Sto_Func<double>::nroot_mfd);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_norm, stopsi->get_pointer(), stosmfpsi, npw, npwx, perbands_sto);
    

        //------------------------  allocate ------------------------
        psi::Psi<std::complex<double>>& expmtsfpsi = sfpsi;
        psi::Psi<std::complex<double>>& expmtsmfpsi = smfpsi;
        psi::Psi<std::complex<double>> exptsfpsi = expmtsfpsi;
        ModuleBase::Memory::record("SDFT::exptsfpsi", memory_cost);
        psi::Psi<std::complex<double>> exptsmfpsi = expmtsmfpsi;
        ModuleBase::Memory::record("SDFT::exptsmfpsi", memory_cost);
        psi::Psi<std::complex<double>> poly_expmtsfpsi, poly_expmtsmfpsi;
        psi::Psi<std::complex<double>> poly_exptsfpsi, poly_exptsmfpsi;
        if (nbatch > 1)
        {
             poly_exptsfpsi.resize(nche_KG, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_exptsfpsi", sizeof(std::complex<double>) * nche_KG * perbands_sto * npwx);

            poly_exptsmfpsi.resize(nche_KG, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_exptsmfpsi", sizeof(std::complex<double>) * nche_KG * perbands_sto * npwx);

            poly_expmtsfpsi.resize(nche_KG, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_expmtsfpsi", sizeof(std::complex<double>) * nche_KG * perbands_sto * npwx);

            poly_expmtsmfpsi.resize(nche_KG, perbands_sto, npwx);
            ModuleBase::Memory::record("SDFT::poly_expmtsmfpsi", sizeof(std::complex<double>) * nche_KG * perbands_sto * npwx);
        }

        //------------------------  t loop  --------------------------
        std::cout << "ik=" << ik << ": ";
        auto start = std::chrono::high_resolution_clock::now();
        const int print_step = ceil(20 / nbatch) * nbatch;
        for (int it = 1; it < nt; ++it)
        {
            // evaluate time cost
            if (it - 1 == print_step)
            {
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                double timeTaken = duration.count();
                std::cout << "(Time left " << timeTaken * (double(nt - 1) / print_step * (nk - ik) - 1) << " s) "<<std::endl;
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
                const std::complex<double> expfactor = exp(0.5*ModuleBase::IMAG_UNIT * eigen * dt);
                const std::complex<double> expmfactor = 1.0 / expfactor;
                for (int ig = 0; ig < npw; ++ig)
                {
                    expmtsfpsi(ib, ig) *= expmfactor;
                    expmtsmfpsi(ib, ig) *= expmfactor;
                    exptsfpsi(ib, ig) *= expfactor;
                    exptsmfpsi(ib, ig) *= expfactor;
                }
            }
            // Sto
            if (nbatch == 1)
            {
                chemt.calfinalvec_complex(&stohchi,
                                         &Stochastic_hchi::hchi_norm,
                                         &expmtsfpsi(perbands_ks, 0),
                                         &expmtsfpsi(perbands_ks, 0),
                                         npw,
                                         npwx,
                                         perbands_sto);
                chemt.calfinalvec_complex(&stohchi,
                                         &Stochastic_hchi::hchi_norm,
                                         &expmtsmfpsi(perbands_ks, 0),
                                         &expmtsmfpsi(perbands_ks, 0),
                                         npw,
                                         npwx,
                                         perbands_sto);
                chet.calfinalvec_complex(&stohchi,
                                         &Stochastic_hchi::hchi_norm,
                                         &exptsfpsi(perbands_ks, 0),
                                         &exptsfpsi(perbands_ks, 0),
                                         npw,
                                         npwx,
                                         perbands_sto);
                chet.calfinalvec_complex(&stohchi,
                                         &Stochastic_hchi::hchi_norm,
                                         &exptsmfpsi(perbands_ks, 0),
                                         &exptsmfpsi(perbands_ks, 0),
                                         npw,
                                         npwx,
                                         perbands_sto);
            }
            else
            {
                std::complex<double>* tmppolyexpmtsfpsi = poly_expmtsfpsi.get_pointer();
                std::complex<double>* tmppolyexpmtsmfpsi = poly_expmtsmfpsi.get_pointer();
                std::complex<double>* tmppolyexptsfpsi = poly_exptsfpsi.get_pointer();
                std::complex<double>* tmppolyexptsmfpsi = poly_exptsmfpsi.get_pointer();
                std::complex<double>* stoexpmtsfpsi = &expmtsfpsi(perbands_ks, 0);
                std::complex<double>* stoexpmtsmfpsi = &expmtsmfpsi(perbands_ks, 0);
                std::complex<double>* stoexptsfpsi = &exptsfpsi(perbands_ks, 0);
                std::complex<double>* stoexptsmfpsi = &exptsmfpsi(perbands_ks, 0);
                if ((it - 1) % nbatch == 0)
                {
                    chet.calpolyvec_complex(&stohchi,
                                            &Stochastic_hchi::hchi_norm,
                                            stoexptsfpsi,
                                            tmppolyexptsfpsi,
                                            npw,
                                            npwx,
                                            perbands_sto);
                    chet.calpolyvec_complex(&stohchi,
                                            &Stochastic_hchi::hchi_norm,
                                            stoexptsmfpsi,
                                            tmppolyexptsmfpsi,
                                            npw,
                                            npwx,
                                            perbands_sto);
                    chemt.calpolyvec_complex(&stohchi,
                                            &Stochastic_hchi::hchi_norm,
                                            stoexpmtsfpsi,
                                            tmppolyexpmtsfpsi,
                                            npw,
                                            npwx,
                                            perbands_sto);
                    chemt.calpolyvec_complex(&stohchi,
                                            &Stochastic_hchi::hchi_norm,
                                            stoexpmtsmfpsi,
                                            tmppolyexpmtsmfpsi,
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
                zgemv_(&transa, &M, &N, &ModuleBase::ONE, tmppolyexptsfpsi, &LDA, tmpcoef, &inc, &ModuleBase::ZERO, stoexptsfpsi, &inc);
                zgemv_(&transa, &M, &N, &ModuleBase::ONE, tmppolyexptsmfpsi, &LDA, tmpcoef, &inc, &ModuleBase::ZERO, stoexptsmfpsi, &inc);
                zgemv_(&transa, &M, &N, &ModuleBase::ONE, tmppolyexpmtsfpsi, &LDA, tmpmcoef, &inc, &ModuleBase::ZERO, stoexpmtsfpsi, &inc);
                zgemv_(&transa, &M, &N, &ModuleBase::ONE, tmppolyexpmtsmfpsi, &LDA, tmpmcoef, &inc, &ModuleBase::ZERO, stoexpmtsmfpsi, &inc);
            }
            ModuleBase::timer::tick(this->classname, "evolution");

            psi::Psi<std::complex<double>>* allexp = &exptsfpsi;
#ifdef __MPI
            allexp = gatherpsi(exptsfpsi, tmppsi_all, npwx, nrecv_ks, displs_ks, nrecv_sto, displs_sto, bandsinfo);
#endif 
            //calculate <\psi| exp(iHt)*J*exp(-iHt) |\psi>
            cal_jmatrix(exptsmfpsi, *allexp, batchj1psi, batchj2psi, j1l, j2l, velop, ik, ModuleBase::IMAG_UNIT, bandsinfo, bsize_psi);
            
            allexp = &expmtsfpsi;
#ifdef __MPI
            allexp = gatherpsi(expmtsfpsi, tmppsi_all, npwx, nrecv_ks, displs_ks, nrecv_sto, displs_sto, bandsinfo);
#endif 
            cal_jmatrix(expmtsmfpsi, *allexp, batchj1psi, batchj2psi, j1r, j2r, velop, ik, ModuleBase::ONE, bandsinfo, bsize_psi);

            // prepare for parallel
            int totnum = ndim * dim_jmatrix;
            int num_per = totnum / GlobalV::NPROC_IN_POOL;
            int st_per = num_per * GlobalV::RANK_IN_POOL;
            int re = totnum % GlobalV::NPROC_IN_POOL;
            if (GlobalV::RANK_IN_POOL < re)
            {
                ++num_per;
                st_per += GlobalV::RANK_IN_POOL;
            }
            else
            {
                st_per += re;
            }
            // Re(i<psi|sqrt(f)j(1-f) exp(iHt)|psi><psi|j exp(-iHt)\sqrt(f)|psi>)
            // Im(l_ij*r_ji)=Re(i l^*_ij*r^+_ij)=Re(i l^*_i*r^+_i)
            // ddot_real = real(A^*_i * B_i)
            ModuleBase::timer::tick(this->classname, "ddot_real");
            ct11[it] += ModuleBase::GlobalFunc::ddot_real(num_per, j1l.c + st_per, j1r.c + st_per, false) * kv.wk[ik] / 2.0;
            double tmp12 = ModuleBase::GlobalFunc::ddot_real(num_per, j1l.c + st_per, j2r.c + st_per, false);
            double tmp21 = ModuleBase::GlobalFunc::ddot_real(num_per, j2l.c + st_per, j1r.c + st_per, false);
            ct12[it] -= 0.5 * (tmp12 + tmp21) * kv.wk[ik] / 2.0;
            ct22[it] += ModuleBase::GlobalFunc::ddot_real(num_per, j2l.c + st_per, j2r.c + st_per, false) * kv.wk[ik] / 2.0;
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