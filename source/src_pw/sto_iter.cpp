#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "global.h"
#include "sto_iter.h"
#include "occupy.h"
#include "../module_base/timer.h" 

double Stochastic_Iter:: mu;
double Stochastic_Iter:: mu0;
double Stochastic_Iter:: Emin;
double Stochastic_Iter:: Emax;

Stochastic_Iter::Stochastic_Iter()
{
    mu0 = 0;
    spolyv = new double [1];
}

Stochastic_Iter::~Stochastic_Iter()
{
    delete[] spolyv;
}

void Stochastic_Iter::init(int &dim, int& chetype)
{
    nchip = GlobalC::sto_wf.nchip;
    stotype = GlobalC::sto_wf.stotype;
    //wait for init
    targetne = GlobalC::CHR.nelec;
    stoche.init( dim, chetype );
    stohchi.init();
    stohchi.get_GRA_index();
    delete [] spolyv;
    int norder = stoche.norder;
    spolyv = new double [norder];

}

void Stochastic_Iter::orthog()
{
    int nkk=1;// We temporarily use gamma k point.
    //orthogonal part
    if(GlobalV::NBANDS > 0)
    {
        for(int ik = 0; ik < nkk; ++ik)
        {
            if(stotype == "pw")
            {
                    stohchi.orthogonal_to_psi_reciprocal(GlobalC::sto_wf.chi0[ik].c,GlobalC::sto_wf.chiortho[ik].c,ik);
            }
            else
            for(int ichi = 0; ichi < nchip; ++ichi)
            {
                std::complex<double> * p0 = &GlobalC::sto_wf.chi0[ik](ichi,0);
                std::complex<double> * pchi = &GlobalC::sto_wf.chiortho[ik](ichi,0);
                stohchi.orthogonal_to_psi_real(p0,pchi,ik);
            }
        }
    }
}

void Stochastic_Iter::checkemm(int &iter)
{
    if(iter == 1)
    {
        Stochastic_hchi:: Emin = GlobalC::sto_wf.emin_sto;
        Stochastic_hchi:: Emax = GlobalC::sto_wf.emax_sto;
    }
    else if(iter > 5)
	{
        return;
	}
        
    int norder = stoche.norder;
       
    //wait for init
    
    std::complex<double> * pchi;
    int ntest = 1;

    if (nchip < ntest) 
	{
		ntest = nchip;
	}

    bool change = false;
    for(int ichi = 0; ichi < ntest; ++ichi)
    {
        if(GlobalV::NBANDS > 0)
        {
            pchi = &GlobalC::sto_wf.chiortho[0](ichi,0);
        }  
        else
        {
            pchi = &GlobalC::sto_wf.chi0[0](ichi,0);
        }
        while(1)
        {
            bool converge;
            if(stotype == "pw")
			{
                converge = stoche.checkconverge(
					stohchi.hchi_reciprocal, 
					pchi, 
					Stochastic_hchi::Emax, 
					Stochastic_hchi::Emin, 
					5);
			}
            else
			{
                converge = stoche.checkconverge(
					stohchi.hchi_real, 
					pchi, 
					Stochastic_hchi::Emax, 
					Stochastic_hchi::Emin, 
					5);
			}

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

    Emax = Stochastic_hchi:: Emax;
    Emin = Stochastic_hchi:: Emin;

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &Emax, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Emin, 1, MPI_DOUBLE, MPI_MIN , MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &change, 1, MPI_CHAR, MPI_LOR , MPI_COMM_WORLD);
#endif

    if(change)
	{
		std::cout<<"New Emax "<<Stochastic_hchi:: Emax<<" ; new Emin "<<Stochastic_hchi:: Emin<<std::endl;
	}

    Stochastic_hchi:: Emin = Emin;
    Stochastic_hchi:: Emax = Emax;
    
}

void Stochastic_Iter::itermu(int &iter) 
{
    ModuleBase::timer::tick("Stochastic_Iter","itermu");
    int nkk=1;// We temporarily use gamma k point.
    double dmu;
    if(iter == 1)
    {
        dmu = 2;
        th_ne = 0.1 * GlobalV::SCF_THR * GlobalC::CHR.nelec;
        std::cout<<"th_ne "<<th_ne<<std::endl;
    }
    else
    {
        dmu = 0.1;
        th_ne = GlobalV::SCF_THR * 1e-2 * GlobalC::CHR.nelec;
    }
    sumpolyval();
    mu = mu0 - dmu;
    double ne1 = calne();
    double mu1 = mu;

    mu = mu0 + dmu;
    double ne2 = calne();
    double mu2 = mu;
    double Dne = th_ne + 1;
    double ne3;
    double mu3;
    
    //test the domin of mu
    /*for(mu = -5; mu<5;mu+=0.2)
    {
        ne3 = calne();
        std::cout<<"mu: "<<mu<<" ; ne: "<<ne3<<std::endl;
    }
    exit(0);*/
    while(ne1 > targetne)
    {
        ne2 = ne1;
        mu2 = mu1;
        mu1 -= dmu;
        mu = mu1;
        ne1 = calne();
        std::cout<<"Reset mu1 from "<<mu1+dmu<<" to "<<mu1<<std::endl;
        dmu *= 2;
    }
    while(ne2 < targetne)
    {
        ne1 = ne2;
        mu1 = mu2;
        mu2 += dmu;
        mu = mu2;
        ne2 = calne();
        std::cout<<"Reset mu2 form "<<mu2-dmu<<" to "<<mu2<<std::endl;
        dmu *= 2;
    }
    int count = 0;
    while(Dne > th_ne)
    {
        mu3 = (mu2 + mu1) / 2;
        mu = mu3;
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
        //std::cout<<std::setw(20)<<"targetne"<<std::setw(20)<<"ne"<<std::setw(20)<<"mu"<<std::setw(20)<<"Dne"<<std::endl;
        //std::cout<<std::setw(20)<<targetne<<std::setw(20)<<ne3<<std::setw(20)<<mu3<<std::setw(20)<<Dne<<std::endl;
        count++;
        if(count > 60)
        {
            std::cout<<"Fermi energy cannot be converged. Set THNE to "<<th_ne<<std::endl;
            th_ne *= 1e4;
        }
    }
    std::cout<<"Converge fermi energy = "<<mu<<" Ry in "<<count<<" steps."<<std::endl;

    GlobalC::en.ef = mu = mu0 = mu3;
    
    //Set GlobalC::wf.wg 
    if(GlobalV::NBANDS > 0)
    {
        for(int ikk = 0; ikk < nkk; ++ikk)
        {
            double *en=GlobalC::wf.ekb[ikk];
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                GlobalC::wf.wg(ikk,iksb) = fd(en[iksb])*GlobalC::kv.wk[ikk];
            }
        }
    }
    ModuleBase::timer::tick("Stochastic_Iter","itermu");
    return;
}

void Stochastic_Iter:: sumpolyval()
{
    ModuleBase::timer::tick("Stochastic_Iter","sumpolyval");
    int norder = stoche.norder;
    //wait for init
    int nkk=1;
    
    ModuleBase::GlobalFunc::ZEROS(spolyv, norder);
    std::complex<double> * pchi;
    for(int ik = 0; ik < nkk; ++ik)
    {
        if(stotype == "pw")
        {
            if(GlobalV::NBANDS > 0)  pchi = GlobalC::sto_wf.chiortho[ik].c; 
            else            pchi = GlobalC::sto_wf.chi0[ik].c;
            stoche.calpolyval(stohchi.hchi_reciprocal, pchi, nchip);
            ModuleBase::GlobalFunc::DCOPY(stoche.polyvalue, spolyv, norder);
        }
        else
        for(int ichi = 0; ichi < nchip; ++ichi)
        {
            if(GlobalV::NBANDS > 0)  pchi = &GlobalC::sto_wf.chiortho[ik](ichi,0);
            else            pchi = &GlobalC::sto_wf.chi0[ik](ichi,0);
            stoche.calpolyval(stohchi.hchi_real , pchi);
            for(int ior = 0; ior < norder; ++ior)
            {
                spolyv[ior] += stoche.polyvalue[ior];
            }
        }
    }
    ModuleBase::timer::tick("Stochastic_Iter","sumpolyval");
    return;
}


double Stochastic_Iter::calne()
{  
    ModuleBase::timer::tick("Stochastic_Iter","calne");
    //wait for init
    int nkk = 1;

    stoche.calcoef(nfd);
    int norder = stoche.norder;
    double totne = 0;
    double sto_ne = 0;
    KS_ne = 0;
    for(int ikk = 0; ikk < nkk; ++ikk)
    {
        double stok_ne = BlasConnector::dot(norder,stoche.coef,1,spolyv,1);
        //double stok_ne = 0;
        //for(int ior = 0; ior < norder; ++ior)
        //{
        //    stok_ne += stoche.coef[ior] * spolyv[ior];
        //    //std::cout<<stoche.coef[ior]<<" "; 
        //}
        //std::cout<<std::endl;
        //std::cout<<"last term "<<stoche.coef[norder-1] * spolyv[norder-1]<<std::endl;
        if(GlobalV::NBANDS > 0)
        {
            double *en=GlobalC::wf.ekb[ikk];
            //number of electrons in KS orbitals
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                KS_ne += fd(en[iksb]) * GlobalC::kv.wk[ikk];
            }
        }
        sto_ne += stok_ne * GlobalC::kv.wk[ikk]; 
    }

#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, &sto_ne, 1, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);
#endif

    totne = KS_ne + sto_ne;
    ModuleBase::timer::tick("Stochastic_Iter","calne");
    return totne;
}

void Stochastic_Iter::sum_stoband()
{  
    ModuleBase::timer::tick("Stochastic_Iter","sum_stoband");
    int nkk=1;// We temporarily use gamma k point.
    int nrxx = stohchi.nrxx;
    int npw = GlobalC::wf.npw;
    int norder = stoche.norder;

    //cal demet
    stoche.calcoef(this->nfdlnfd);
    double stodemet = 0;

    for(int ikk = 0; ikk < nkk; ++ikk)
    {
        double stok_demet = BlasConnector::dot(norder,stoche.coef,1,spolyv,1);
        //double stok_demet=0;
        //for(int ior = 0; ior < norder; ++ior)
        //{
        //    stok_demet += stoche.coef[ior] * spolyv[ior];
        //}
        //std::cout<<"last term "<<stoche.coef[norder-1] * spolyv[norder-1]<<std::endl;
        //std::cout<<std::endl;
        if(GlobalV::NBANDS > 0)
        {
            double *enb=GlobalC::wf.ekb[ikk];
            //number of electrons in KS orbitals
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                GlobalC::en.demet += fdlnfd(enb[iksb]) * GlobalC::kv.wk[ikk];
            }
        }
        stodemet += stok_demet * GlobalC::kv.wk[ikk];
    }

    //cal eband & rho
    stoche.calcoef(this->nroot_fd);
    //stoche.calcoef(this->nfd);
    
    std::complex<double> * out, *hout;
    double *sto_rho = new double [nrxx];
    std::complex<double> *sto_rhog;
    int npwall = npw * nchip;
    if(stotype == "pw")
    {
        if(GlobalV::NBANDS == 0)
            out = new std::complex<double> [npwall];
        sto_rhog = new std::complex<double> [npw];
    }
    else
    {
        out = new std::complex<double> [nrxx];
    }

    double dr3 = GlobalC::ucell.omega / GlobalC::pw.ncxyz;
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
   
    double tmprho, tmpne;
    std::complex<double> outtem;
    double sto_ne = 0;
    double sto_eband = 0;
    ModuleBase::GlobalFunc::ZEROS(sto_rho, nrxx);

    std::complex<double> * pchi;
    std::complex<double>* porter = GlobalC::UFFT.porter;
    int* GRA_index = stohchi.GRA_index;
    double out2;

    double *ksrho;
    if(GlobalV::NBANDS > 0 && GlobalV::MY_POOL==0 && stotype == "pw")
    {
        ksrho = new double [nrxx];
        ModuleBase::GlobalFunc::DCOPY(GlobalC::CHR.rho[0],ksrho,nrxx);
        ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[0],nrxx);
    }
    
    for(int ik = 0; ik < nkk; ++ik)
    {
        double stok_eband =0;
        if(stotype == "pw")
        {
                if(GlobalV::NBANDS > 0)
                    out = pchi = GlobalC::sto_wf.chiortho[ik].c;
                else
                    pchi = GlobalC::sto_wf.chi0[ik].c;
                
                stoche.calfinalvec(stohchi.hchi_reciprocal, pchi, out, nchip);
                hout = new std::complex<double> [npwall];
                stohchi.hchi_reciprocal(out,hout,nchip);
                stok_eband = ModuleBase::GlobalFunc::ddot_real(npwall, out, hout,false) * DeltaE 
                            +  ModuleBase::GlobalFunc::ddot_real(npwall, out, out,false) * Ebar;
                sto_eband += stok_eband * GlobalC::kv.wk[ik];
                std::complex<double> *tmpout = out;
            for(int ichi = 0; ichi < nchip ; ++ichi)
            {
                ModuleBase::GlobalFunc::ZEROS( porter, GlobalC::pw.nrxx );
                for(int ig = 0; ig < npw; ++ig)
                {
                    porter[ GRA_index[ig] ] = tmpout[ig];
                }
                GlobalC::pw.FFT_wfc.FFT3D(GlobalC::UFFT.porter, 1);
                for(int ir = 0 ; ir < nrxx ; ++ir)
                {
                    GlobalC::CHR.rho[0][ir] += norm(porter[ir]);
                }
                tmpout+=npw;
            }
#ifdef __MPI
            GlobalC::CHR.rho_mpi();
#endif
            for(int ir = 0; ir < nrxx ; ++ir)
            {
                tmprho = GlobalC::CHR.rho[0][ir] * GlobalC::kv.wk[ik] / GlobalC::ucell.omega;
                sto_rho[ir] = tmprho;
                sto_ne += tmprho;
            }
            sto_ne *= dr3;
        }
        else
        {
            hout = new std::complex<double> [nrxx];
            for(int ichi = 0; ichi < nchip; ++ichi)
            {
                if(GlobalV::NBANDS > 0)
                {
                    pchi = &GlobalC::sto_wf.chiortho[ik](ichi,0);
                }  
                else
                {
                    pchi = &GlobalC::sto_wf.chi0[ik](ichi,0);
                }
                stoche.calfinalvec(stohchi.hchi_real, pchi, out);
                stohchi.hchi_real(out,hout);
                for(int ir = 0; ir < nrxx; ++ir)
                {
                    out2 = norm(out[ir]);
                    //out2 = real(conj(pchi[ir]) * out[ir]);
                    tmpne = out2 * GlobalC::kv.wk[ik];
                    stok_eband += real(conj(out[ir]) * hout[ir]) * DeltaE + Ebar * out2;
                    sto_rho[ir] += tmpne; 
                    sto_ne += tmpne;
                } 
            }
            sto_eband += stok_eband * GlobalC::kv.wk[ik];
        }
        delete [] hout;
           
    }



#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE,&stodemet,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_eband,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_ne,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,sto_rho,nrxx,MPI_DOUBLE,MPI_SUM,PARAPW_WORLD);
#endif
    GlobalC::en.eband += sto_eband;
    GlobalC::en.demet += stodemet;
    GlobalC::en.demet *= Occupy::gaussian_parameter;

    std::cout.precision(12);
    std::cout<<"Renormalize rho from ne = "<<sto_ne+KS_ne<<" to targetne = "<<targetne<<std::endl;


    double factor;
    if(stotype == "pw")
    {
        if(abs(sto_ne) > 1e-20)
            factor = (targetne - KS_ne) / sto_ne;
        else
            factor = 1;
    }
    else
    {
        if(abs(sto_ne) > 1e-20)
            factor = (targetne - KS_ne) / sto_ne / dr3;
        else
            factor = 1 / dr3;
    }
    
    if(GlobalV::MY_POOL==0 && stotype == "pw")
    {
        if(GlobalV::NBANDS > 0)
            ModuleBase::GlobalFunc::DCOPY(ksrho,GlobalC::CHR.rho[0],nrxx);
        else
            ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[0],nrxx);
    }
    
    
    if(GlobalV::MY_POOL == 0)
    for(int is = 0 ; is < 1; ++is)
    {
        for(int ir = 0; ir < nrxx ; ++ir)
        {
            GlobalC::CHR.rho[is][ir] += sto_rho[ir] * factor;
        }
    }

    
    
    


    if(stotype == "pw")
    {
        delete [] sto_rhog;
    }
    delete [] sto_rho;
    if(GlobalV::NBANDS == 0 || stotype != "pw")
        delete [] out;
    ModuleBase::timer::tick("Stochastic_Iter","sum_stoband");
    return;
}

double Stochastic_Iter:: root_fd(double e)
{
    double e_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(e_mu > 72)
        return 0;
    else
        return 1 / sqrt(1 + exp(e_mu));
}

double Stochastic_Iter:: nroot_fd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double ne_mu = (e * DeltaE + Ebar - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 72)
        return 0;
    else
        return 1 / sqrt(1 + exp(ne_mu));
}

double Stochastic_Iter:: fd(double e)
{
    double e_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(e_mu > 36)
        return 0;
    else
        return 1 / (1 + exp(e_mu));
}

double Stochastic_Iter:: nfd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double ne_mu = (e * DeltaE + Ebar - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 36)
        return 0;
    else
        return 1 / (1 + exp(ne_mu));
}

double Stochastic_Iter:: fdlnfd(double e)
{
    double e_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(e_mu > 36)
        return 0;
    else if(e_mu < -36)
        return 0;
    else
    {
        double f = 1 / (1 + exp(e_mu));
        return (f * log(f) + (1.0-f) * log(1.0-f)); 
    }
}

double Stochastic_Iter:: nfdlnfd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double ne_mu = (e * DeltaE + Ebar - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 36)
        return 0;
    else if(ne_mu < -36)
        return 0;
    else
    {
        double f = 1 / (1 + exp(ne_mu));
        return f * log(f) + (1-f) * log(1-f); 
    }
}


 void Stochastic_Iter:: test()
 {
     //=====================test============================
    /*
    std::complex<double> *in = new std::complex<double> [GlobalC::pw.nrxx];
    
    std::complex<double> *chig1 = new std::complex<double> [GlobalC::wf.npw];
    std::complex<double> *chig2 = new std::complex<double> [GlobalC::wf.npw];
    ModuleBase::GlobalFunc::ZEROS(in,GlobalC::pw.nrxx);
    ModuleBase::GlobalFunc::ZEROS(in2,GlobalC::pw.nrxx);*/

    //---------------------------------------------------
    /*//test hermit property of  hchi matrix
    Emin = -1;
    Emax = 1;
    Stochastic_hchi:: Emin = this -> Emin;
    Stochastic_hchi:: Emax = this -> Emax;
    std::complex<double> *out = new std::complex<double> [GlobalC::pw.nrxx];
    std::complex<double> *in2 = new std::complex<double> [GlobalC::pw.nrxx];
    std::cout<<"------------------------------------"<<std::endl;
    std::complex<double> cij,cji;
    double dc;
    for(int i = 0 ; i < 300 ; ++i)
    {
        if( i % 10  == 0)
            std::cout<<"We are testing "<<i+1 <<" rows..."<<std::endl;
        for(int j = i+1; j < 300 ; ++j)
        {
            in2[j] = 1;
            stohchi.hchi_real(in2,out);
            cij = out[i];
            in2[j] = 0;
            in2[i] = 1;
            stohchi.hchi_real(in2,out);
            cji = out[j];
            in2[i] = 0;
            dc = abs(conj(cij)-cji);
            if(dc > 1e-6)
            {
                std::cout<<"(i,j) = ("<<i+1<<" , "<<j+1<<") ; cij = "<<cij<<" ; cji = "<<cji<<std::endl;
            }

        }
    }
    std::cout<<"------------------------------------"<<std::endl;
    delete[] out;
    delete[] in2;*/
    //---------------------------------------------------
    
    //-------------------------------------------------------
    //compare hchi_reciprocal and h_psi
    /*Emin = -1;
    Emax = 1;
    Stochastic_hchi:: Emin = this -> Emin;
    Stochastic_hchi:: Emax = this -> Emax;
    std::complex<double> *chig1 = new std::complex<double> [GlobalC::wf.npw];
    std::complex<double> *chig2 = new std::complex<double> [GlobalC::wf.npw];
    std::complex<double> *kswf = &GlobalC::wf.evc[0](0,0);
    
    
    stohchi.hchi_reciprocal(kswf,chig1);
    
    GlobalC::hm.hpw.h_psi( kswf , chig2);
    if(GlobalV::MY_RANK==0)
    for(int i = 0; i<GlobalC::wf.npw;++i)
    {
        if(i % 100 == 0)
            std::cout<<kswf[i]<<" "<<chig1[i]<<" "<<chig2[i]<<std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(GlobalV::MY_RANK==1)
    for(int i = 0; i<GlobalC::wf.npw;++i)
    {
        std::cout.clear();
        if(i % 100 == 0)
            std::cout<<kswf[i]<<" "<<chig1[i]<<" "<<chig2[i]<<std::endl;
    }
    std::cout<<std::endl;*/
    //---------------------------------------------------

     //-------------------------------------------------------
    /*//compare hchi_real and h_psi
    Emin = -1;
    Emax = 1;
    Stochastic_hchi:: Emin = this -> Emin;
    Stochastic_hchi:: Emax = this -> Emax;
    int * GRA_index = stohchi.GRA_index;
    ModuleBase::GlobalFunc::ZEROS(in,GlobalC::pw.nrxx);
    std::complex<double> *kswf = &GlobalC::wf.evc[0](0,0);
    for ( int ig = 0 ; ig< GlobalC::wf.npw; ++ig)
    {
        in[GRA_index[ig]] = kswf [ig];
    }
    fftw_plan pp=fftw_plan_dft_3d(GlobalC::pw.nx,GlobalC::pw.ny,GlobalC::pw.nz,(fftw_complex *)in,(fftw_complex *)in2, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pp);
    stohchi.hchi_real(in2,out);
    fftw_plan pp2=fftw_plan_dft_3d(GlobalC::pw.nx,GlobalC::pw.ny,GlobalC::pw.nz,(fftw_complex *)out,(fftw_complex *)out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(pp2);
    for(int i  = 0; i <GlobalC::wf.npw;++i)
    {
        chig1[i] = out[GRA_index[i]] / GlobalC::pw.nrxx; 
    }

    GlobalC::hm.hpw.h_psi( kswf , chig2);
    for(int i = 0; i<GlobalC::wf.npw;++i)
    {
        if(i % 100 == 0)
        std::cout<<kswf[i]<<" "<<in[GRA_index[i]]<<" "<<chig1[i]<<" "<<chig2[i]<<std::endl;
    }
    std::cout<<std::endl;*/
    //---------------------------------------------------
    
    //---------------------------------------------------
    //compare eigen energy
    /*
    int * GRA_index = stohchi.GRA_index;
    std::complex<double> *chigout = new std::complex<double> [GlobalC::wf.npw];
    std::complex<double> *wave = new std::complex<double> [GlobalC::pw.nrxx];
    std::complex<double> *waveout = new std::complex<double> [GlobalC::pw.nrxx];
    Emax = 1000;
    Emin = 0;
    Stochastic_hchi:: Emin = this->Emin;
    Stochastic_hchi:: Emax = this->Emax;
    double Ebar = (Emax + Emin)/2;
    double DeltaE = (Emax - Emin)/2;
    fftw_plan pp=fftw_plan_dft_3d(GlobalC::pw.nx,GlobalC::pw.ny,GlobalC::pw.nz,(fftw_complex *)wave,(fftw_complex *)wave, FFTW_BACKWARD, FFTW_ESTIMATE);
    for(int ib = 0 ; ib < GlobalV::NBANDS ; ++ib)
    {
        std::complex<double> *kswf = &GlobalC::wf.evc[0](ib,0);
        GlobalC::hm.hpw.h_psi( kswf , chigout);
        double energy = 0;
        double norm1 =0;
        ModuleBase::GlobalFunc::ZEROS(wave,GlobalC::pw.nrxx);
        for(int ig = 0 ; ig < GlobalC::wf.npw ; ++ig)
        {
            energy += real(conj(kswf[ig]) * chigout[ig]);
            norm1 += norm (kswf[ig]);
            wave[GRA_index[ig]] = kswf[ig];
        }
        fftw_execute(pp);
        stohchi.hchi_real(wave,waveout);
        double energy2 = 0;
        double norm2 =0;
        for(int ir = 0 ; ir < GlobalC::pw.nrxx ; ++ir)
        {
            energy2 += real(conj(wave[ir]) * waveout[ir]) * DeltaE + Ebar * norm(wave[ir]);
            norm2 += norm(wave[ir]);
        }


        std::cout<<"BAND "<<ib+1<<" "<<energy/norm1<<"  "<<energy2/norm2<<std::endl;
    }
    delete []chigout;
    delete []wave;
    delete []waveout;*/
    //---------------------------------------------------

    //---------------------------------------------------
    //test ne
    /*
    int * GRA_index = stohchi.GRA_index;
    std::complex<double> *wave = new std::complex<double> [GlobalC::pw.nrxx];
    std::complex<double> *waveout = new std::complex<double> [GlobalC::pw.nrxx];
    Emax = 750;
    Emin = -100;
    Stochastic_hchi:: Emin = this->Emin;
    Stochastic_hchi:: Emax = this->Emax;
    mu = GlobalC::en.ef;
    stoche.calcoef(this->nfd);
    fftw_plan pp=fftw_plan_dft_3d(GlobalC::pw.nx,GlobalC::pw.ny,GlobalC::pw.nz,(fftw_complex *)wave,(fftw_complex *)wave, FFTW_BACKWARD, FFTW_ESTIMATE);
    for(int ib = 0 ; ib < GlobalV::NBANDS ; ++ib)
    {
        ModuleBase::GlobalFunc::ZEROS(wave,GlobalC::pw.nrxx);
        std::complex<double> *kswf = &GlobalC::wf.evc[0](ib,0);
        for(int ig = 0 ; ig < GlobalC::wf.npw ; ++ig)
        {
            wave[GRA_index[ig]] = kswf[ig];
        }
        fftw_execute(pp);
        double ne =0;
        double norm1 = 0;
        stoche.calresult(stohchi.hchi_real, GlobalC::pw.nrxx, wave, waveout);
        for(int ir = 0 ; ir < GlobalC::pw.nrxx ; ++ir)
        {   
            ne += real(conj(wave[ir]) * waveout[ir]);
            norm1 += norm(wave[ir]);
        }
        std::cout<<"Ne of Band "<<ib+1<<" is "<<ne/norm1 * GlobalC::kv.wk[0]<<std::endl;
    }
    delete []wave;
    delete []waveout;*/

    //-------------------------------------------------------------
	//test orthogonal
    /*int npw=GlobalC::wf.npw;
    char transC='C';
    char transN='N';
    int nchip = GlobalC::sto_wf.nchip;
    std::complex<double> *wave = new std::complex<double> [nchip*npw];
    for(int i = 0; i < nchip*npw; ++i)
    {
        wave[i]=GlobalC::sto_wf.chi0[0].c[i];
    }
    std::complex<double> *sum = new std::complex<double> [nchip * nchip];
	zgemm_(&transC, &transN, &nchip, &nchip, &npw, &ONE, GlobalC::sto_wf.chi0[0].c, &npw, wave, &npw, &ZERO, sum, &nchip);
	Parallel_Reduce::reduce_complex_double_pool(sum, nchip * nchip);
	if(GlobalV::MY_RANK!=0) std::cout.clear();
    double abs2 = 0;
    for(int i=0;i<nchip * nchip;++i)
    {
        if(i%nchip==int(i%nchip)) continue;
        abs2+=norm(sum[i]);
    }
	std::cout<<abs2/nchip<<std::endl;
    delete [] sum;
    delete [] wave;
	if(GlobalV::MY_RANK!=0) std::cout.setstate(ios::failbit);*/
	//-------------------------------------------------------------

    
    //=====================================================
 }
