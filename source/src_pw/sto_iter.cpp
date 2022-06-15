#include "global.h"
#include "sto_iter.h"
#include "occupy.h"
#include "diago_cg.h" 
#include "../module_base/tool_title.h"
#include "../module_base/tool_quit.h"
#include "../module_base/timer.h"
#include "../src_parallel/parallel_reduce.h"

Stochastic_Iter::Stochastic_Iter()
{
    p_che = nullptr;
    spolyv = new double [1];
    change = false;
    mu0 = 0;
}

Stochastic_Iter::~Stochastic_Iter()
{
    if(p_che != nullptr) delete p_che;
    delete[] spolyv;
}

void Stochastic_Iter::init(const int dim, int* nchip_in)
{
    p_che = new ModuleBase::Chebyshev<double>(INPUT.nche_sto);
    nchip = nchip_in;
    targetne = GlobalC::CHR.nelec;
    stohchi.init();
    delete [] spolyv;
    const int norder = p_che->norder;
    spolyv = new double [norder];
    stofunc.Emin = INPUT.emin_sto;
    stofunc.Emax = INPUT.emax_sto;
}

void Stochastic_Iter::orthog(const int& ik, Stochastic_WF& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter","orthog");
    //orthogonal part
    if(GlobalV::NBANDS > 0)
    {
	    const int nchipk=stowf.nchip[ik];
	    const int npw = GlobalC::wf.npw;
	    const int npwx = GlobalC::wf.npwx;
        std::complex<double> *wfgin = stowf.chi0[ik].c, *wfgout = stowf.chiortho[ik].c;
	    for(int ig = 0 ; ig < npwx * nchipk; ++ig)
	    {
		    wfgout[ig] = wfgin[ig];
	    }

	    //orthogonal part
	    std::complex<double> *sum = new std::complex<double> [GlobalV::NBANDS * nchipk];
	    char transC='C';
	    char transN='N';
    
	    //sum(b<NBANDS, a<nchi) = < psi_b | chi_a >
	    zgemm_(&transC, &transN, &GlobalV::NBANDS, &nchipk, &npw, &ModuleBase::ONE, 
                GlobalC::wf.evc[ik].c, &npwx, wfgout, &npwx, &ModuleBase::ZERO, sum, &GlobalV::NBANDS);
	    Parallel_Reduce::reduce_complex_double_pool(sum, GlobalV::NBANDS * nchipk);
    
	    //psi -= psi * sum
	    zgemm_(&transN, &transN, &npw, &nchipk, &GlobalV::NBANDS, &ModuleBase::NEG_ONE, 
                GlobalC::wf.evc[ik].c, &npwx, sum, &GlobalV::NBANDS, &ModuleBase::ONE, wfgout, &npwx);
	    delete[] sum;
    }
}

void Stochastic_Iter::checkemm(const int& ik, int &iter, Stochastic_WF& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter","checkemm");
    if(iter > 5)
	{
        return;
	}
        
    const int norder = p_che->norder;
    std::complex<double> * pchi;
    int ntest = 1;

    if (nchip[ik] < ntest) 
	{
		ntest = nchip[ik];
	}

    for(int ichi = 0; ichi < ntest; ++ichi)
    {
        if(GlobalV::NBANDS > 0)
        {
            pchi = &stowf.chiortho[ik](ichi,0);
        }  
        else
        {
            pchi = &stowf.chi0[ik](ichi,0);
        }
        while(1)
        {
            bool converge;
            converge = p_che->checkconverge(
				&stohchi, &Stochastic_hchi::hchi_reciprocal, 
				pchi, GlobalC::wf.npw,
				stohchi.Emax, 
				stohchi.Emin, 
				5.0);

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
    if(ik == GlobalC::kv.nks-1)
    {
        stofunc.Emax = stohchi.Emax;
        stofunc.Emin = stohchi.Emin;

#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, &stofunc.Emax, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &stofunc.Emin, 1, MPI_DOUBLE, MPI_MIN , MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &change, 1, MPI_CHAR, MPI_LOR , MPI_COMM_WORLD);
#endif
        stohchi.Emin = stofunc.Emin;
        stohchi.Emax = stofunc.Emax;
        if(change)
	    {
	    	GlobalV::ofs_running<<"New Emax "<<stohchi.Emax<<" ; new Emin "<<stohchi.Emin<<std::endl;
	    }
        change = false;
    }
}

void Stochastic_Iter::itermu(int &iter) 
{
    ModuleBase::TITLE("Stochastic_Iter","itermu");
    ModuleBase::timer::tick("Stochastic_Iter","itermu");
    double dmu;
    if(iter == 1)
    {
        dmu = 2;
        th_ne = 0.1 * GlobalV::SCF_THR * GlobalC::CHR.nelec;
        // std::cout<<"th_ne "<<th_ne<<std::endl;
    }
    else
    {
        dmu = 0.1;
        th_ne = GlobalV::SCF_THR * 1e-2 * GlobalC::CHR.nelec;
    }
    this->stofunc.mu = mu0 - dmu;
    double ne1 = calne();
    double mu1 = this->stofunc.mu;

    this->stofunc.mu = mu0 + dmu;
    double ne2 = calne();
    double mu2 = this->stofunc.mu;
    double Dne = th_ne + 1;
    double ne3;
    double mu3;
    
    while(ne1 > targetne)
    {
        ne2 = ne1;
        mu2 = mu1;
        mu1 -= dmu;
        this->stofunc.mu = mu1;
        ne1 = calne();
        std::cout<<"Reset mu1 from "<<mu1+dmu<<" to "<<mu1<<std::endl;
        dmu *= 2;
    }
    while(ne2 < targetne)
    {
        ne1 = ne2;
        mu1 = mu2;
        mu2 += dmu;
        this->stofunc.mu = mu2;
        ne2 = calne();
        // cout<<"Reset mu2 from "<<mu2-dmu<<" to "<<mu2<<endl;
        dmu *= 2;
    }
    int count = 0;
    while(Dne > th_ne)
    {
        mu3 = (mu2 + mu1) / 2;
        this->stofunc.mu = mu3;
        ne3 = calne();
        if(ne3 < targetne)
        {
            ne1 = ne3;
            mu1 = mu3;
        }
        else if(ne3 > targetne)
        {
            ne2 = ne3;
            mu2 = mu3;
        }
        Dne = abs(targetne - ne3);

        count++;
        if(count > 60)
        {
            std::cout<<"Fermi energy cannot be converged. Set THNE to "<<th_ne<<std::endl;
            th_ne *= 1e1;
            if(th_ne > 1e1) ModuleBase::WARNING_QUIT("Stochastic_Iter",
                                         "Cannot converge feimi energy. Please retry with different random number");
        }
    }
    GlobalV::ofs_running<<"Converge fermi energy = "<<this->stofunc.mu<<" Ry in "<<count<<" steps."<<std::endl;
    //precision check
    double tmpre;
    tmpre = p_che->coef_real[p_che->norder-1] * spolyv[p_che->norder-1];
    MPI_Allreduce(MPI_IN_PLACE, &tmpre, 1, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);
    GlobalV::ofs_running<<"Chebyshev Precision: "<<abs(tmpre/targetne)*1e9<<"E-09"<<std::endl;
    if(tmpre/targetne > GlobalV::SCF_THR )
    {
        stringstream ss;
        ss<<tmpre/targetne;
        string fractxt,tartxt,itertxt;
        ss>>fractxt;
        ss.clear();
        ss<<GlobalV::SCF_THR;
        ss>>tartxt;
        ss.clear();
        ss<<iter+1;
        ss>>itertxt;
        string warningtxt = "Iter "+itertxt+": (Chebyshev error = "+fractxt+" > threshold = "+tartxt+" ) Please add more expansion terms for Chebychev expansion.";
        ModuleBase::WARNING("Stochastic_Chebychev", warningtxt);
    }


    GlobalC::en.ef = this->stofunc.mu = mu0 = mu3;
    
    //Set wf.wg 
    if(GlobalV::NBANDS > 0)
    {
        for(int ikk = 0; ikk < GlobalC::kv.nks; ++ikk)
        {
            double *en = GlobalC::wf.ekb[ikk];
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                GlobalC::wf.wg(ikk,iksb) = stofunc.fd(en[iksb])*GlobalC::kv.wk[ikk];
            }
        }
    }
    ModuleBase::timer::tick("Stochastic_Iter","itermu");
    return;
}

void Stochastic_Iter:: sumpolyval_k(const int& ik, Stochastic_WF& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter","sumpolyval_k");
    ModuleBase::timer::tick("Stochastic_Iter","sumpolyval");
    const int norder = p_che->norder;
    if(ik==0)   ModuleBase::GlobalFunc::ZEROS(spolyv, norder);
    std::complex<double> * pchi;
    if(GlobalV::NBANDS > 0)  pchi = stowf.chiortho[ik].c; 
    else            pchi = stowf.chi0[ik].c;
    p_che->tracepolyA(&stohchi, &Stochastic_hchi::hchi_reciprocal, pchi, GlobalC::wf.npw, GlobalC::wf.npwx, nchip[ik]);
    for(int i = 0 ; i < norder ; ++i)
    {
        spolyv[i] += p_che->polytrace[i] * GlobalC::kv.wk[ik];
    }
    ModuleBase::timer::tick("Stochastic_Iter","sumpolyval");
    return;
}


double Stochastic_Iter::calne()
{  
    ModuleBase::timer::tick("Stochastic_Iter","calne");

    p_che->calcoef_real(&stofunc,&Sto_Func<double>::nfd);
    const int norder = p_che->norder;
    double totne = 0;
    KS_ne = 0;
    double sto_ne = BlasConnector::dot(norder,p_che->coef_real,1,spolyv,1);
    if(GlobalV::NBANDS > 0)
    {
        for(int ikk = 0; ikk < GlobalC::kv.nks; ++ikk)
        {
            double *en=GlobalC::wf.ekb[ikk];
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                KS_ne += stofunc.fd(en[iksb]) * GlobalC::kv.wk[ikk];
            }
        }
    }
    KS_ne /= GlobalV::NPROC_IN_POOL;
	MPI_Allreduce(MPI_IN_PLACE, &KS_ne, 1, MPI_DOUBLE, MPI_SUM , STO_WORLD);

#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, &sto_ne, 1, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);
#endif

    totne = KS_ne + sto_ne;
    ModuleBase::timer::tick("Stochastic_Iter","calne");
    return totne;
}

void Stochastic_Iter::sum_stoband(Stochastic_WF& stowf)
{  
    ModuleBase::TITLE("Stochastic_Iter","sum_stoband");
    ModuleBase::timer::tick("Stochastic_Iter","sum_stoband");
    int nrxx = GlobalC::wfcpw->nrxx;
    int npwx = GlobalC::wf.npwx;
    const int norder = p_che->norder;

    //cal demet
    p_che->calcoef_real(&stofunc,&Sto_Func<double>::nfdlnfd);
    double stodemet = BlasConnector::dot(norder,p_che->coef_real,1,spolyv,1);

    if(GlobalV::NBANDS > 0)
    {
        for(int ikk = 0; ikk < GlobalC::kv.nks; ++ikk)
        {
            double *enb=GlobalC::wf.ekb[ikk];
            //number of electrons in KS orbitals
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                GlobalC::en.demet += stofunc.fdlnfd(enb[iksb]) * GlobalC::kv.wk[ikk];
            }
        }
    }
    GlobalC::en.demet /= GlobalV::NPROC_IN_POOL;
	MPI_Allreduce(MPI_IN_PLACE, &GlobalC::en.demet, 1, MPI_DOUBLE, MPI_SUM , STO_WORLD);

    //cal eband
    p_che->calcoef_real(&stofunc,&Sto_Func<double>::nxfd);
    double sto_eband = BlasConnector::dot(norder,p_che->coef_real,1,spolyv,1);

    //cal rho
    p_che->calcoef_real(&stofunc,&Sto_Func<double>::nroot_fd);
    
    std::complex<double> * out, *hout;
    double *sto_rho = new double [nrxx];
    //int npwall = npwx * nchip;

    double dr3 = GlobalC::ucell.omega / GlobalC::wfcpw->nxyz;
    double tmprho, tmpne;
    std::complex<double> outtem;
    double sto_ne = 0;
    ModuleBase::GlobalFunc::ZEROS(sto_rho, nrxx);

    std::complex<double> * pchi;
    std::complex<double>* porter = new std::complex<double>[nrxx];
    double out2;

    double *ksrho = nullptr;
    if(GlobalV::NBANDS > 0 && GlobalV::MY_STOGROUP==0 )
    {
        ksrho = new double [nrxx];
        ModuleBase::GlobalFunc::DCOPY(GlobalC::CHR.rho[0],ksrho,nrxx);
        ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[0],nrxx);
    }
    
    for(int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        //init k
        if(GlobalC::kv.nks > 1) GlobalC::hm.hpw.init_k(ik);
        stohchi.current_ik = ik;

        const int npw = GlobalC::kv.ngk[ik];
        double stok_eband;
        out = stowf.shchi[ik].c;
        if(GlobalV::NBANDS > 0)
            pchi = stowf.chiortho[ik].c;
        else
            pchi = stowf.chi0[ik].c;
        
        p_che->calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_reciprocal, pchi, out, npw, GlobalC::wf.npwx, nchip[ik]);

        std::complex<double> *tmpout = out;
        for(int ichi = 0; ichi < nchip[ik] ; ++ichi)
        {
            GlobalC::wfcpw->recip2real(tmpout, porter, ik);
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                GlobalC::CHR.rho[0][ir] += norm(porter[ir]) * GlobalC::kv.wk[ik];
            }
            tmpout+=npwx;
        }
    }
    delete[] porter;
   
    GlobalC::CHR.rho_mpi();
    for(int ir = 0; ir < nrxx ; ++ir)
    {
        tmprho = GlobalC::CHR.rho[0][ir] / GlobalC::ucell.omega;
        sto_rho[ir] = tmprho;
        sto_ne += tmprho;
    }
    sto_ne *= dr3;



#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE,&stodemet,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_eband,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_ne,1,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_ne,1,MPI_DOUBLE,MPI_SUM,PARAPW_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,sto_rho,nrxx,MPI_DOUBLE,MPI_SUM,PARAPW_WORLD);
#endif
    GlobalC::en.eband += sto_eband;
    GlobalC::en.demet += stodemet;
    GlobalC::en.demet *= Occupy::gaussian_parameter;

    cout.precision(12);
    GlobalV::ofs_running<<"Renormalize rho from ne = "<<sto_ne+KS_ne<<" to targetne = "<<targetne<<endl;

    double factor;
    if(abs(sto_ne) > 1e-20)
        factor = (targetne - KS_ne) / sto_ne;
    else
        factor = 1;

    if(GlobalV::MY_STOGROUP==0)
    {
        if(GlobalV::NBANDS > 0)
            ModuleBase::GlobalFunc::DCOPY(ksrho,GlobalC::CHR.rho[0],nrxx);
        else
            ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[0],nrxx);
    }
    
    
    if(GlobalV::MY_STOGROUP == 0)
    for(int is = 0 ; is < 1; ++is)
    {
        for(int ir = 0; ir < nrxx ; ++ir)
        {
            GlobalC::CHR.rho[is][ir] += sto_rho[ir] * factor;
        }
    }

    
    delete [] sto_rho;
    delete [] ksrho;
    ModuleBase::timer::tick("Stochastic_Iter","sum_stoband");
    return;
}


