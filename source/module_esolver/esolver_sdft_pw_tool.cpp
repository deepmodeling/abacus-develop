#include "./esolver_sdft_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw_sdft.h"
#include "module_base/timer.h"
#include "module_base/constants.h"
#include "module_base/vector3.h"
#include "module_base/complexmatrix.h"
#include "module_base/global_variable.h"
#include "module_base/global_function.h"
#include "module_base/memory.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/velocity_pw.h"

#define TWOSQRT2LN2 2.354820045030949 //FWHM = 2sqrt(2ln2) * \sigma
#define FACTOR 1.839939223835727e7
namespace ModuleESolver
{
void parallelks(int &ksbandper, int &startband)
{
    ksbandper = GlobalV::NBANDS / GlobalV::NSTOGROUP;
    startband = ksbandper * GlobalV::MY_STOGROUP;
    if(GlobalV::MY_STOGROUP < GlobalV::NBANDS % GlobalV::NSTOGROUP)
    {  
        ++ksbandper;
        startband += GlobalV::MY_STOGROUP;
    }
    else
    {
        startband += GlobalV::NBANDS % GlobalV::NSTOGROUP;
    }
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
    int ntest = 2;
    for (int ik = 0;ik < nk; ++ik)
	{
        this->p_hamilt->updateHk(ik);
        stoiter.stohchi.current_ik = ik;
        const int npw = kv.ngk[ik];
        std::complex<double> *pchi = new std::complex<double> [npw];
        for(int i = 0; i < ntest ; ++i)
        {
            for(int ig = 0; ig < npw; ++ig)
            {
                double rr = std::rand()/double(RAND_MAX);
                double arg = std::rand()/double(RAND_MAX);
                pchi[ig] = std::complex<double>(rr * cos(arg), rr * sin(arg));
            }
            while(1)
            {
                bool converge;
                converge= chetest.checkconverge(&stohchi, &Stochastic_hchi::hchi_norm, 
	        	    	pchi, npw, stohchi.Emax, stohchi.Emin, 5.0);
                
                if(!converge)
	        	{
                    change = true;
	        	}
                else
	        	{
                    break;
	        	}
            }
        }
        delete[] pchi;

        if(ik == nk-1)
        {
            stoiter.stofunc.Emax = stohchi.Emax;
            stoiter.stofunc.Emin = stohchi.Emin;
#ifdef __MPI
            MPI_Allreduce(MPI_IN_PLACE, &stoiter.stofunc.Emax, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &stoiter.stofunc.Emin, 1, MPI_DOUBLE, MPI_MIN , MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &change, 1, MPI_CHAR, MPI_LOR , MPI_COMM_WORLD);
#endif
            stohchi.Emin = stoiter.stofunc.Emin;
            stohchi.Emax = stoiter.stofunc.Emax;
            if(change)
	        {
	        	GlobalV::ofs_running<<"New Emax "<<stohchi.Emax<<" ; new Emin "<<stohchi.Emin<<std::endl;
	        }
            change = false;
        }
    }
}

void ESolver_SDFT_PW::sKG(const int nche_KG, const double fwhmin, const double wcut, 
                          const double dw_in, const double dt_in, const int nbatch)
{
     ModuleBase::TITLE(this->classname,"sKG");
    ModuleBase::timer::tick(this->classname,"sKG");
    cout<<"Calculating conductivity...."<<endl;
    
    int nw = ceil(wcut/dw_in);
    double dw =  dw_in / ModuleBase::Ry_to_eV; //converge unit in eV to Ry 
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    double dt = dt_in; //unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    int nt = ceil(sqrt(20)/sigma/dt);
    cout<<"nw: "<<nw<<" ; dw: "<<dw*ModuleBase::Ry_to_eV<<" eV"<<endl;
    cout<<"nt: "<<nt<<" ; dt: "<<dt<<" a.u.(ry^-1)"<<endl;
    assert(nw >= 1);
    assert(nt >= 1);
    const int ndim = 3;

    ModuleBase::Chebyshev<double> che(this->nche_sto);
    ModuleBase::Chebyshev<double> chet(nche_KG);
    ModuleBase::Chebyshev<double> chet2(nche_KG);
    const int npwx = wf.npwx;
    const double tpiba = GlobalC::ucell.tpiba;
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    Stochastic_hchi& stohchi = stoiter.stohchi;
    const int nk = kv.nks;

    //------------------------------------------------------------------
    //                    Calculate
    //------------------------------------------------------------------
    const double mu = this->pelec->eferm.ef;
    stoiter.stofunc.mu = mu;
    double * ct11 = new double[nt];
    double * ct12 = new double[nt];
    double * ct22 = new double[nt];
    ModuleBase::GlobalFunc::ZEROS(ct11,nt);
    ModuleBase::GlobalFunc::ZEROS(ct12,nt);
    ModuleBase::GlobalFunc::ZEROS(ct22,nt);
    stoiter.stofunc.t = dt * nbatch;
    chet.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::nsin);
    chet2.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::n_sin);
    cout<<"Relative Chebyshev precision: "<<abs(chet.coef_complex[nche_KG-1]/chet.coef_complex[0])<<endl;
    ofstream cheofs("Chebycoef");
    for(int i  = 0 ; i < nche_KG ; ++i)
    {
        cheofs<<setw(5)<<i<<setw(20)<<abs(chet.coef_complex[i]/chet.coef_complex[0])<<endl;
    }
    cheofs.close();
    std::complex<double>* batchcoef = nullptr;
    if(nbatch > 1)
    {
        batchcoef = new std::complex<double> [nche_KG * nbatch];
        std::complex<double>* tmpcoef = batchcoef + (nbatch - 1) * nche_KG;
        for(int i = 0 ; i < nche_KG ; ++i)
        {
            tmpcoef[i] = chet.coef_complex[i];
        }
        for(int ib = 0 ; ib < nbatch - 1; ++ib)
        {
            tmpcoef = batchcoef + ib * nche_KG;
            stoiter.stofunc.t = dt * (ib + 1);
            chet.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::nsin);
            for(int i = 0 ; i < nche_KG ; ++i)
            {
                tmpcoef[i] = chet.coef_complex[i];
            }
        }
        stoiter.stofunc.t = dt * nbatch;
    }

    
    //parallel KS orbitals
    int ksbandper,startband;
    parallelks(ksbandper,startband);
    ModuleBase::timer::tick(this->classname,"kloop");
    hamilt::Velocity velop(pw_wfc, kv.isk.data(), &GlobalC::ppcell, &GlobalC::ucell, INPUT.cond_nonlocal);
    for (int ik = 0;ik < nk;++ik)
	{
        velop.init(ik);
        if(nk > 1) 
        {
            this->p_hamilt->updateHk(ik);
        }
        stoiter.stohchi.current_ik = ik;
        const int npw = kv.ngk[ik];

        int nchip = this->stowf.nchip[ik];
        int totbands_per = nchip + ksbandper;
        int totbands_ks = GlobalV::NBANDS;
        int totbands_sto = nchip;
        int totbands = totbands_per;
#ifdef __MPI
        MPI_Allreduce(&totbands_per,&totbands,1,MPI_INT,MPI_SUM,PARAPW_WORLD);
        totbands_sto = totbands - totbands_ks;
        const int nstogroup = GlobalV::NSTOGROUP;
        int nrecv_ks[nstogroup];
        int nrecv_sto[nstogroup];
        int displs_ks[nstogroup];
        int displs_sto[nstogroup];
        MPI_Allgather(&ksbandper,1,MPI_INT,nrecv_ks,1,MPI_INT,PARAPW_WORLD);
        MPI_Allgather(&nchip,1,MPI_INT,nrecv_sto,1,MPI_INT,PARAPW_WORLD);
        displs_ks[0] = 0;
        displs_sto[0] = 0;
        for(int i = 1; i < nstogroup ; ++i)
        {
            displs_ks[i] = displs_ks[i-1] + nrecv_ks[i-1];
            displs_sto[i] = displs_sto[i-1] + nrecv_sto[i-1];
        }
        for(int i = 0 ; i < nstogroup ; ++i)
        {
            nrecv_ks[i] *= npwx;
            nrecv_sto[i] *= npwx;
            displs_ks[i] *= npwx;
            displs_sto[i] *= npwx;
        }
#endif
        
        //-----------------------------------------------------------
        //               ks conductivity
        //-----------------------------------------------------------
        if(GlobalV::MY_STOGROUP == 0 && totbands_ks > 0)
            jjcorr_ks(ik, nt, dt, this->pelec->wg, velop, ct11,ct12,ct22);
        
        //-----------------------------------------------------------
        //               sto conductivity
        //-----------------------------------------------------------
        //before loop

        //|psi>
        size_t memory_cost = totbands_per * npwx * sizeof(std::complex<double>);
        psi::Psi<std::complex<double>> psi0(1, totbands_per, npwx, kv.ngk.data()); //|psi>
        ModuleBase::Memory::record("SDFT::psi0", memory_cost);
        psi::Psi<std::complex<double>> sfpsi0(1, totbands_per, npwx, kv.ngk.data()); // sqrt(f)|psi>
        ModuleBase::Memory::record("SDFT::sfpsi0", memory_cost);
        psi::Psi<std::complex<double>> hpsi0(1, totbands_per, npwx, kv.ngk.data()); // h|psi>
        ModuleBase::Memory::record("SDFT::hpsi0", memory_cost);
        psi::Psi<std::complex<double>> hsfpsi0(1, totbands_per, npwx, kv.ngk.data()); // h*sqrt(f)|psi>
        ModuleBase::Memory::record("SDFT::hsfpsi0", memory_cost);
        //j|psi> j1=p  j2=(Hp+pH)/2 - mu*p
        memory_cost = ndim * totbands_per * npwx * sizeof(std::complex<double>);
        psi::Psi<std::complex<double>> j1psi(1, ndim * totbands_per, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::j1psi", memory_cost);
        psi::Psi<std::complex<double>> j2psi(1, ndim * totbands_per, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::j2psi", memory_cost);
        //(1-f)*j*sqrt(f)|psi>
        psi::Psi<std::complex<double>> j1sfpsi(1, ndim * totbands_per, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::psi0", memory_cost);
        psi::Psi<std::complex<double>> j2sfpsi(1, ndim * totbands_per, npwx, kv.ngk.data());
        ModuleBase::Memory::record("SDFT::psi0", memory_cost);
        double* en;
        if(ksbandper > 0)   en = new double [ksbandper];
        for(int ib = 0 ; ib < ksbandper ; ++ib)
        {
            en[ib] = this->pelec->ekb(ik, ib+startband);
        }
        for(int ib = 0 ; ib < totbands_per ; ++ib )
        {
            std::complex<double> *tmp, *tmp2;
            double fac = 1;
            
            if(ib < ksbandper)
            {
                tmp = &(this->psi[0](ik,ib+startband,0));
                tmp2 = &(this->psi[0](ik,ib+startband,0));
                fac = sqrt(stoiter.stofunc.fd(en[ib]));
            }
            else
            {
                if(totbands_ks > 0)
                    tmp = &stowf.chiortho[ik](ib-ksbandper,0);
                else
                    tmp = &stowf.chi0[ik](ib-ksbandper,0);
                tmp2 = &stowf.shchi[ik](ib-ksbandper,0);
            }
        
            for(int ig = 0; ig < npw ; ++ig)
            {
                psi0(ib,ig) = tmp[ig];
                sfpsi0(ib,ig) = fac * tmp2[ig];
            }
        }
        

        //init j1psi,j2psi,j1sfpsi,j2sfpsi
        velop.act(this->psi, totbands_per, psi0.get_pointer(), j1psi.get_pointer());
        velop.act(this->psi, totbands_per, sfpsi0.get_pointer(), j1sfpsi.get_pointer());

        // this->p_hamilt->hPsi(psi0.get_pointer(), hpsi0.get_pointer(), totbands_per*npwx);
        // this->p_hamilt->hPsi(sfpsi0.get_pointer(), hsfpsi0.get_pointer(), totbands_per*npwx);
        // this->p_hamilt->hPsi(j1psi.get_pointer(), j2psi.get_pointer(), ndim*totbands_per*npwx);
        // this->p_hamilt->hPsi(j1sfpsi.get_pointer(), j2sfpsi.get_pointer(), ndim*totbands_per*npwx);
        psi::Range allbands(1,0,0,totbands_per-1);
        hamilt::Operator<std::complex<double>>::hpsi_info info_psi0(&psi0, allbands, hpsi0.get_pointer());
        this->p_hamilt->ops->hPsi(info_psi0);

        hamilt::Operator<std::complex<double>>::hpsi_info info_sfpsi0(&sfpsi0, allbands, hsfpsi0.get_pointer());
        this->p_hamilt->ops->hPsi(info_sfpsi0);

        psi::Range allndimbands(1,0,0,ndim*totbands_per-1);
        hamilt::Operator<std::complex<double>>::hpsi_info info_j1psi(&j1psi, allndimbands, j2psi.get_pointer());
        this->p_hamilt->ops->hPsi(info_j1psi);

        hamilt::Operator<std::complex<double>>::hpsi_info info_j1sfpsi(&j1sfpsi, allndimbands, j2sfpsi.get_pointer());
        this->p_hamilt->ops->hPsi(info_j1sfpsi);

        /*
        // stohchi.hchi_norm(psi0.get_pointer(), hpsi0.get_pointer(), totbands_per);
        // stohchi.hchi_norm(sfpsi0.get_pointer(), hsfpsi0.get_pointer(), totbands_per);
        // stohchi.hchi_norm(j1psi.get_pointer(), j2psi.get_pointer(), ndim*totbands_per);
        // stohchi.hchi_norm(j1sfpsi.get_pointer(), j2sfpsi.get_pointer(), ndim*totbands_per);
        // double Ebar = (stohchi.Emin + stohchi.Emax)/2;
	    // double DeltaE = (stohchi.Emax - stohchi.Emin)/2;
	    // for(int ib = 0 ; ib < totbands_per ; ++ib)
	    // {
	    // 	for(int ig = 0; ig < npw; ++ig)
	    // 	{
        //         hpsi0(ib,ig) = hpsi0(ib,ig) * DeltaE + Ebar * psi0(ib,ig);
	    // 		hsfpsi0(ib,ig) = hsfpsi0(ib,ig) * DeltaE + Ebar * sfpsi0(ib,ig);
	    // 	}
	    // }
        // for(int id = 0 ; id < ndim ; ++id)
        // {
        //     for(int ib = 0 ; ib < totbands_per ; ++ib)
	    //     {
	    //     	for(int ig = 0; ig < npw; ++ig)
	    //     	{
        //             const int idib = id * totbands_per + ib;
        //             j2psi(0,idib,ig) = j2psi(0,idib,ig) * DeltaE + Ebar * j1psi(0,idib,ig);
	    //     		j2sfpsi(0,idib,ig) = j2sfpsi(0,idib,ig) * DeltaE + Ebar * j1sfpsi(0,idib,ig);
	    //     	}
	    //     }
        // }*/

        velop.act(this->psi, totbands_per, hpsi0.get_pointer(), j2psi.get_pointer(), true);
        velop.act(this->psi, totbands_per, hsfpsi0.get_pointer(), j2sfpsi.get_pointer(), true);

        for(int idib = 0 ; idib < ndim * totbands_per; ++idib)
        {
            for(int ig = 0; ig < npw ; ++ig)
            {
                j2psi(0,idib,ig) = j2psi(0,idib,ig)/2.0 - mu * j1psi(0,idib,ig);
                j2sfpsi(0,idib,ig) = j2sfpsi(0,idib,ig)/2.0 - mu * j1sfpsi(0,idib,ig);
            }
        }

        //(1-f)
        che.calcoef_real(&stoiter.stofunc,&Sto_Func<double>::n_fd);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_norm, j1sfpsi.get_pointer(), j1sfpsi.get_pointer(), npw, npwx, totbands_per*ndim);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_norm, j2sfpsi.get_pointer(), j2sfpsi.get_pointer(), npw, npwx, totbands_per*ndim);
        
        psi::Psi<std::complex<double>> *p_j1psi = &j1psi; 
        psi::Psi<std::complex<double>> *p_j2psi = &j2psi; 
#ifdef __MPI
        psi::Psi<std::complex<double>> j1psi_tot,j2psi_tot;
        if (GlobalV::NSTOGROUP > 1)
        {
            j1psi_tot.resize(1,ndim*totbands,npwx);
            ModuleBase::Memory::record("SDFT::j1psi_tot", ndim*totbands*npwx);
            j2psi_tot.resize(1,ndim*totbands,npwx);
            ModuleBase::Memory::record("SDFT::j2psi_tot", ndim*totbands*npwx);
            for(int id = 0 ; id < ndim ; ++id)
            {
                const int idnb_per = id * totbands_per;
                const int idnb = id * totbands;
                MPI_Allgatherv(&j1psi(idnb_per,0), ksbandper * npwx, MPI_DOUBLE_COMPLEX, 
                            &j1psi_tot(idnb,0), nrecv_ks, displs_ks, MPI_DOUBLE_COMPLEX, PARAPW_WORLD);
                MPI_Allgatherv(&j2psi(idnb_per,0), ksbandper * npwx, MPI_DOUBLE_COMPLEX, 
                            &j2psi_tot(idnb,0), nrecv_ks, displs_ks, MPI_DOUBLE_COMPLEX, PARAPW_WORLD);
                MPI_Allgatherv(&j1psi(idnb_per + ksbandper,0), nchip * npwx, MPI_DOUBLE_COMPLEX, 
                            &j1psi_tot(idnb + totbands_ks, 0), nrecv_sto, displs_sto, MPI_DOUBLE_COMPLEX, PARAPW_WORLD);
                MPI_Allgatherv(&j2psi(idnb_per + ksbandper,0), nchip * npwx, MPI_DOUBLE_COMPLEX, 
                            &j2psi_tot(idnb + totbands_ks,0), nrecv_sto, displs_sto, MPI_DOUBLE_COMPLEX, PARAPW_WORLD);
            }
            p_j1psi = &j1psi_tot;
            p_j2psi = &j2psi_tot;
        }
#endif

        //loop of t
        psi::Psi<std::complex<double>> exppsi(1,totbands_per,npwx);
        ModuleBase::Memory::record("SDFT::exppsi", sizeof(std::complex<double>) * totbands_per * npwx);
        psi::Psi<std::complex<double>> expsfpsi(1,totbands_per,npwx);
        ModuleBase::Memory::record("SDFT::expsfpsi", sizeof(std::complex<double>) * totbands_per * npwx);
        psi::Psi<std::complex<double>> *poly_psi = nullptr, *poly_sfpsi = nullptr;
        if(nbatch > 1)
        {
            poly_psi = new psi::Psi<std::complex<double>>(nche_KG, nchip, npwx);
            ModuleBase::Memory::record("SDFT::poly_psi", sizeof(std::complex<double>) * nche_KG * nchip * npwx);
            poly_sfpsi = new psi::Psi<std::complex<double>>(nche_KG, nchip, npwx);
            ModuleBase::Memory::record("SDFT::poly_sfpsi", sizeof(std::complex<double>) * nche_KG * nchip * npwx);
        }
        for(int ib = 0; ib < totbands_per; ++ib)
        {
            for(int ig = 0 ; ig < npw ; ++ig)
            {
                exppsi(ib,ig) = psi0(ib,ig);
                expsfpsi(ib,ig) = sfpsi0(ib,ig);
            }
        }
        cout<<"ik="<<ik<<": ";
        for (int it = 1 ;it < nt ; ++it)
        {
            if(it%20==0) cout<<it<<" ";
            ModuleBase::timer::tick(this->classname,"evolution_ks");
            for(int ib = 0; ib < ksbandper; ++ib)
            {
                for(int ig = 0 ; ig < npw ; ++ig)
                {
                    double eigen = en[ib];
                    //exp(iHdt)|kschi>
                    exppsi(ib,ig) *= exp(ModuleBase::IMAG_UNIT*eigen*dt);
                    //exp(-iHdt)|shkschi>
                    expsfpsi(ib,ig) *= exp(ModuleBase::NEG_IMAG_UNIT*eigen*dt);
                }
            }
            ModuleBase::timer::tick(this->classname,"evolution_ks");
            
            ModuleBase::timer::tick(this->classname,"evolution_sto");
            if(nbatch == 1)
            {
                //exp(iHdt)|chi>
                chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_norm, &exppsi(ksbandper,0), &exppsi(ksbandper,0), npw, npwx, nchip);
                //exp(-iHdt)|shchi>
                chet2.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_norm, &expsfpsi(ksbandper,0), &expsfpsi(ksbandper,0), npw, npwx, nchip);
            }
            else
            {
                std::complex<double>* tmppolypsi = poly_psi->get_pointer();
                std::complex<double>* tmppolysfpsi = poly_sfpsi->get_pointer();
                std::complex<double>* psisto_c = &exppsi(ksbandper,0);
                std::complex<double>* sfpsisto_c = &expsfpsi(ksbandper,0);
                if((it - 1) % nbatch == 0)
                {
                    //exp(iHdt)|chi>
                    chet.calpolyvec_complex(&stohchi, &Stochastic_hchi::hchi_norm, psisto_c, tmppolypsi, npw, npwx, nchip);
                    //exp(-iHdt)|shchi>
                    chet2.calpolyvec_complex(&stohchi, &Stochastic_hchi::hchi_norm, sfpsisto_c, tmppolysfpsi, npw, npwx, nchip);
                }
                
                std::complex<double>* tmpcoef = batchcoef + (it - 1) % nbatch * nche_KG;
                ModuleBase::GlobalFunc::ZEROS(  psisto_c, npwx*nchip);
                ModuleBase::GlobalFunc::ZEROS(sfpsisto_c, npwx*nchip);
                for(int iche = 0 ; iche < nche_KG ; ++iche)
                {
                    std::complex<double> *tmp1 =   psisto_c;
                    std::complex<double> *tmp2 = sfpsisto_c;
                    std::complex<double> coef = tmpcoef[iche];
                    for(int ib = 0 ; ib < nchip; ++ib)
                    {
                        for(int ig = 0 ; ig < npw ; ++ig)
                        {
                            tmp1[ig] += coef * tmppolypsi[ig];
                            tmp2[ig] += conj(coef) * tmppolysfpsi[ig];
                        }
                        tmppolypsi += npwx;
                        tmppolysfpsi += npwx;
                        tmp1 += npwx;
                        tmp2 += npwx;
                    }
                }
            }
            ModuleBase::timer::tick(this->classname,"evolution_sto");
            psi::Psi<std::complex<double>> *p_exppsi = &exppsi;
#ifdef __MPI
            psi::Psi<std::complex<double>> exppsi_tot;
            if (GlobalV::NSTOGROUP > 1)
            {
                ModuleBase::timer::tick(this->classname,"bands_gather");
                exppsi_tot.resize(1,totbands,npwx);
                MPI_Allgatherv(&exppsi(0,0), ksbandper* npwx, MPI_DOUBLE_COMPLEX, 
                                &exppsi_tot(0,0), nrecv_ks, displs_ks, MPI_DOUBLE_COMPLEX, PARAPW_WORLD);
                MPI_Allgatherv(&exppsi(ksbandper,0), nchip * npwx, MPI_DOUBLE_COMPLEX, 
                                &exppsi_tot(totbands_ks,0), nrecv_sto, displs_sto, MPI_DOUBLE_COMPLEX, PARAPW_WORLD);
                p_exppsi = &exppsi_tot;
                ModuleBase::timer::tick(this->classname,"bands_gather");
            }
#endif
            const int dim_jmatrix = totbands_per*totbands_sto + nchip*totbands_ks;
            // const int dim_jmatrix = totbands_per*totbands;
            ModuleBase::ComplexMatrix j1l(ndim,dim_jmatrix), j2l(ndim,dim_jmatrix);
            ModuleBase::ComplexMatrix j1r(ndim,dim_jmatrix), j2r(ndim,dim_jmatrix);
            char transa = 'C';
            char transb = 'N';
            int totbands_per3 = ndim*totbands_per;
            int totbands3 = ndim*totbands;
            ModuleBase::timer::tick(this->classname,"matrix_multip");
            for(int id = 0 ; id < ndim ; ++id)
            {
                const int idnb = id * totbands_per;
                //<psi|sqrt(f)j_1(1-f) exp(iHt)|psi>
                zgemm_(&transa, &transb, &totbands_per, &totbands_sto, &npw, &ModuleBase::ONE, &j1sfpsi(idnb,0), &npwx, 
                        &(p_exppsi->operator()(totbands_ks,0)), &npwx, &ModuleBase::ZERO, &j1l(id,0), &totbands_per);
                if(totbands_ks > 0)
                    zgemm_(&transa, &transb, &nchip, &totbands_ks, &npw, &ModuleBase::ONE, &j1sfpsi(idnb+ksbandper, 0), &npwx, 
                        p_exppsi->get_pointer(), &npwx, &ModuleBase::ZERO, &j1l(id,totbands_per*totbands_sto), &nchip);
                // if(totbands_ks*ksbandper > 0)
                //     zgemm_(&transa, &transb, &ksbandper, &totbands_ks, &npw, &ModuleBase::ONE, &j1sfpsi(idnb,0), &npwx, 
                //         p_exppsi->get_pointer(), &npwx, &ModuleBase::ZERO, &j1l(id,totbands_per*totbands_sto + nchip*totbands_ks), &ksbandper);
                //<psi|sqrt(f)j_2(1-f) exp(iHt)|psi>
                zgemm_(&transa, &transb, &totbands_per, &totbands_sto, &npw, &ModuleBase::ONE, &j2sfpsi(idnb,0), &npwx, 
                        &(p_exppsi->operator()(totbands_ks,0)), &npwx, &ModuleBase::ZERO, &j2l(id,0), &totbands_per);
                if(totbands_ks > 0)
                    zgemm_(&transa, &transb, &nchip, &totbands_ks, &npw, &ModuleBase::ONE, &j2sfpsi(idnb+ksbandper,0), &npwx, 
                        p_exppsi->get_pointer(), &npwx, &ModuleBase::ZERO, &j2l(id,totbands_per*totbands_sto), &nchip);
                // if(totbands_ks*ksbandper > 0)
                //     zgemm_(&transa, &transb, &ksbandper, &totbands_ks, &npw, &ModuleBase::ONE, &j2sfpsi(idnb,0), &npwx, 
                //         p_exppsi->get_pointer(), &npwx, &ModuleBase::ZERO, &j2l(id,totbands_per*totbands_sto + nchip*totbands_ks), &ksbandper);
            }
            for(int id = 0; id < ndim ; ++id) // it can also use zgemm once
            {
                const int idnb = id * totbands;
                //i<exp(-iHt)sqrt(f)psi| j_1|psi> = i<psi|sqrt(f)exp(iHt) j_1|psi> = i(<psi|j_1 exp(-iHt)sqrt(f)|psi>)^+
                zgemm_(&transa, &transb, &totbands_per, &totbands_sto, &npw, &ModuleBase::IMAG_UNIT, expsfpsi.get_pointer(), &npwx,
                        &(p_j1psi->operator()(idnb+totbands_ks,0)), &npwx, &ModuleBase::ZERO, &j1r(id,0), &totbands_per);
                if(totbands_ks > 0)
                    zgemm_(&transa, &transb, &nchip, &totbands_ks, &npw, &ModuleBase::IMAG_UNIT, &expsfpsi(ksbandper,0), &npwx,
                        &(p_j1psi->operator()(idnb,0)), &npwx, &ModuleBase::ZERO, &j1r(id,totbands_per*totbands_sto), &nchip);
                // if(totbands_ks*ksbandper > 0)
                //     zgemm_(&transa, &transb, &ksbandper, &totbands_ks, &npw, &ModuleBase::IMAG_UNIT, expsfpsi.get_pointer(), &npwx,
                //         &(p_j1psi->operator()(idnb,0)), &npwx, &ModuleBase::ZERO, &j1r(id,totbands_per*totbands_sto + nchip*totbands_ks), &ksbandper);
                //i<psi|sqrt(f)exp(iHt) j_2|psi> 
                zgemm_(&transa, &transb,&totbands_per, &totbands_sto, &npw, &ModuleBase::IMAG_UNIT, expsfpsi.get_pointer(), &npwx,
                        &(p_j2psi->operator()(idnb+totbands_ks,0)), &npwx, &ModuleBase::ZERO, &j2r(id,0), &totbands_per);
                if(totbands_ks > 0)   
                    zgemm_(&transa, &transb,&nchip, &totbands_ks, &npw, &ModuleBase::IMAG_UNIT, &expsfpsi(ksbandper,0), &npwx,
                        &(p_j2psi->operator()(idnb,0)), &npwx, &ModuleBase::ZERO, &j2r(id,totbands_per*totbands_sto), &nchip);
                // if(totbands_ks*ksbandper > 0)
                //     zgemm_(&transa, &transb,&ksbandper, &totbands_ks, &npw, &ModuleBase::IMAG_UNIT, expsfpsi.get_pointer(), &npwx,
                //         &(p_j2psi->operator()(idnb,0)), &npwx, &ModuleBase::ZERO, &j2r(id,totbands_per*totbands_sto + nchip*totbands_ks), &ksbandper);
            }
            ModuleBase::timer::tick(this->classname,"matrix_multip");

#ifdef __MPI
            ModuleBase::timer::tick(this->classname,"matrix_reduce");
            MPI_Allreduce(MPI_IN_PLACE,j1l.c,ndim*dim_jmatrix,MPI_DOUBLE_COMPLEX,MPI_SUM,POOL_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,j2l.c,ndim*dim_jmatrix,MPI_DOUBLE_COMPLEX,MPI_SUM,POOL_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,j1r.c,ndim*dim_jmatrix,MPI_DOUBLE_COMPLEX,MPI_SUM,POOL_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,j2r.c,ndim*dim_jmatrix,MPI_DOUBLE_COMPLEX,MPI_SUM,POOL_WORLD);
            ModuleBase::timer::tick(this->classname,"matrix_reduce");
#endif
            int totnum = ndim*dim_jmatrix;
            int num_per = totnum / GlobalV::NPROC_IN_POOL;
            int st_per = num_per * GlobalV::RANK_IN_POOL;
            int re = totnum % GlobalV::NPROC_IN_POOL;
            if(GlobalV::RANK_IN_POOL < re)  
            {
                ++num_per;
                st_per += GlobalV::RANK_IN_POOL;
            }
            else
            {
                st_per += re;
            }
            //Re(i<psi|sqrt(f)j(1-f) exp(iHt)|psi><psi|j exp(-iHt)\sqrt(f)|psi>)
            //Im(l_ij*r_ji)=Re(i l^*_ij*r^+_ij)=Re(i l^*_i*r^+_i)
            //ddot_real = real(A^*_i * B_i)
            ModuleBase::timer::tick(this->classname,"ddot_real");
            ct11[it]
                += ModuleBase::GlobalFunc::ddot_real(num_per, j1l.c + st_per, j1r.c + st_per, false) * kv.wk[ik] / 2,
                0;
            ct12[it]
                -= ModuleBase::GlobalFunc::ddot_real(num_per, j1l.c + st_per, j2r.c + st_per, false) * kv.wk[ik] / 2,
                0;
            ct22[it]
                += ModuleBase::GlobalFunc::ddot_real(num_per, j2l.c + st_per, j2r.c + st_per, false) * kv.wk[ik] / 2,
                0;
            ModuleBase::timer::tick(this->classname,"ddot_real");
            
        }
        cout<<endl;
        if(ksbandper > 0)   delete[] en;
        delete poly_psi;
        delete poly_sfpsi;
    }
    ModuleBase::timer::tick(this->classname,"kloop");
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE,ct11,nt,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ct12,nt,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ct22,nt,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
   
    //------------------------------------------------------------------
    //                    Output
    //------------------------------------------------------------------
    if(GlobalV::MY_RANK==0)
    {
        calcondw(nt,dt,fwhmin,wcut,dw_in,ct11,ct12,ct22);
    }
    delete[] ct11;
    delete[] ct12;
    delete[] ct22;
    delete[] batchcoef;
    // delete[] cterror;
    ModuleBase::timer::tick(this->classname,"sKG");
}

void ESolver_SDFT_PW:: caldos( const int nche_dos, const double sigmain, const double emin, const double emax, const double de, const int npart)
{
    ModuleBase::TITLE(this->classname,"caldos");
    ModuleBase::timer::tick(this->classname,"caldos");
    cout<<"========================="<<endl;
    cout<<"###Calculating Dos....###"<<endl;
    cout<<"========================="<<endl;
    ModuleBase::Chebyshev<double> che(nche_dos);
    const int nk = kv.nks;
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    Stochastic_hchi& stohchi = stoiter.stohchi;
    const int npwx = wf.npwx;

    double * spolyv = nullptr;
    std::complex<double> *allorderchi = nullptr;
    if(stoiter.method == 1)
    {
        spolyv = new double [nche_dos];
        ModuleBase::GlobalFunc::ZEROS(spolyv, nche_dos);
    }
    else
    {
        spolyv = new double [nche_dos*nche_dos];
        ModuleBase::GlobalFunc::ZEROS(spolyv, nche_dos*nche_dos);
        int nchip_new = ceil((double)this->stowf.nchip_max / npart);
        allorderchi = new std::complex<double> [nchip_new * npwx * nche_dos];
    }
    ModuleBase::timer::tick(this->classname,"Tracepoly");
    cout<<"1. TracepolyA:"<<endl;
    for (int ik = 0;ik < nk;ik++)
	{
        cout<<"ik: "<<ik+1<<endl;
		if(nk > 1) 
        {
            this->p_hamilt->updateHk(ik);
        }
        stohchi.current_ik = ik;
        const int npw = kv.ngk[ik];
        const int nchipk = this->stowf.nchip[ik];
        
        std::complex<double> * pchi;
        if(GlobalV::NBANDS > 0)
            pchi = stowf.chiortho[ik].c;
        else
            pchi = stowf.chi0[ik].c;
        if(stoiter.method == 1)
        {
            che.tracepolyA(&stohchi, &Stochastic_hchi::hchi_norm, pchi, npw, npwx, nchipk);
            for(int i = 0 ; i < nche_dos ; ++i)
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
            for(int ipart = 0 ; ipart < npart ; ++ipart)
            {
                int nchipk_new = nchipk / npart;
                int start_nchipk = ipart * nchipk_new + nchipk % npart;
                if(ipart < nchipk % npart)
                {
                    nchipk_new++;
                    start_nchipk = ipart * nchipk_new;
                }
                ModuleBase::GlobalFunc::ZEROS(allorderchi, nchipk_new * npwx * nche_dos);
                std::complex<double> *tmpchi = pchi + start_nchipk * npwx;
                che.calpolyvec_complex(&stohchi, &Stochastic_hchi::hchi_norm, tmpchi, allorderchi, npw, npwx, nchipk_new);
                double* vec_all= (double *) allorderchi;
                int LDA = npwx * nchipk_new * 2;
                int M = npwx * nchipk_new * 2;
                dgemm_(&trans,&normal, &N,&N,&M,&kweight,vec_all,&LDA,vec_all,&LDA,&one,spolyv,&N);
            }
        }
    }
    if(stoiter.method == 2) delete[] allorderchi;

    ofstream ofsdos;
    int ndos = int((emax-emin) / de)+1;
    stoiter.stofunc.sigma = sigmain / ModuleBase::Ry_to_eV;
    ModuleBase::timer::tick(this->classname,"Tracepoly");

    cout<<"2. Dos:"<<endl;
    ModuleBase::timer::tick(this->classname,"DOS Loop");
    int n10 = ndos/10;
    int percent = 10;
    double *sto_dos = new double [ndos];
    double *ks_dos = new double [ndos];
    double *error = new double [ndos];
	for(int ie = 0; ie < ndos; ++ie)
	{
        double tmpks = 0;
        double tmpsto = 0;
        stoiter.stofunc.targ_e = (emin + ie * de) / ModuleBase::Ry_to_eV;
        if(stoiter.method == 1)
        {
            che.calcoef_real(&stoiter.stofunc, &Sto_Func<double>::ngauss);
		    tmpsto = BlasConnector::dot(nche_dos,che.coef_real,1,spolyv,1);
        }
        else
        {
            che.calcoef_real(&stoiter.stofunc, &Sto_Func<double>::nroot_gauss);
            tmpsto = stoiter.vTMv(che.coef_real,spolyv,nche_dos);
        }
        if(GlobalV::NBANDS > 0)
        {
            for(int ik = 0; ik < nk; ++ik)
            {
                double *en=&(this->pelec->ekb(ik, 0));
                for(int ib = 0; ib < GlobalV::NBANDS; ++ib)
                {
                    tmpks += stoiter.stofunc.gauss(en[ib]) * kv.wk[ik] / 2;
                }
            }
        }
        tmpks /= GlobalV::NPROC_IN_POOL;

        double tmperror = 0;
        if(stoiter.method == 1)
        {
            tmperror = che.coef_real[nche_dos-1] * spolyv[nche_dos-1];
        }
        else
        {
            const int norder = nche_dos;
            double last_coef = che.coef_real[norder-1];
            double last_spolyv = spolyv[norder*norder - 1];
            tmperror = last_coef *(BlasConnector::dot(norder,che.coef_real,1,spolyv+norder*(norder-1),1)
                        + BlasConnector::dot(norder,che.coef_real,1,spolyv+norder-1,norder)-last_coef*last_spolyv);
        } 

        if(ie%n10 == n10 -1) 
        {
            cout<<percent<<"%"<<" ";
            percent+=10;
        }
        sto_dos[ie] = tmpsto;
        ks_dos[ie] = tmpks;
        error[ie] = tmperror;
	}
#ifdef __MPI
	    MPI_Allreduce(MPI_IN_PLACE, ks_dos, ndos, MPI_DOUBLE, MPI_SUM , STO_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, sto_dos, ndos, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, error, ndos, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);
#endif
    if(GlobalV::MY_RANK == 0)
    {
        string dosfile = GlobalV::global_out_dir+"DOS1_smearing.dat";
        ofsdos.open(dosfile.c_str());
        double maxerror = 0;
        double sum = 0; 
        ofsdos<<setw(8)<<"## E(eV) "<<setw(20)<<"dos(eV^-1)"<<setw(20)<<"sum"<<setw(20)<<"Error(eV^-1)"<<endl;
        for(int ie = 0 ; ie < ndos ; ++ie)
        {
            double tmperror = 2.0 * abs(error[ie]);
            if(maxerror < tmperror) maxerror = tmperror;
            double dos = 2.0 * (ks_dos[ie] + sto_dos[ie]) / ModuleBase::Ry_to_eV;
            sum += dos;
	    	ofsdos <<setw(8)<< emin + ie * de <<setw(20)<< dos <<setw(20)<< sum * de <<setw(20)<< tmperror <<endl;
        }
        cout<<endl;
        cout<<"Finish DOS"<<endl;
        cout<<scientific<<"DOS max absolute Chebyshev Error: "<<maxerror<<endl;
        ofsdos.close();
    }
    delete[] sto_dos;
    delete[] ks_dos;
    delete[] error;
    delete[] spolyv;
    ModuleBase::timer::tick(this->classname,"DOS Loop");
    ModuleBase::timer::tick(this->classname,"caldos");
    return;
}

}//namespace ModuleESolver

namespace GlobalTemp
{

const ModuleBase::matrix* veff;

}