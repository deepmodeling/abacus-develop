#include "./esolver_sdft_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw_sdft.h"
#include "module_base/timer.h"
#include "module_base/constants.h"
#include "module_base/vector3.h"
#include "module_base/complexmatrix.h"
#include "module_base/global_variable.h"
#include "module_base/global_function.h"
#include "src_pw/global.h"

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

void ESolver_SDFT_PW::sKG(const int nche_KG, const double fwhmin, const double wcut, 
                          const double dw_in, const int times)
{
    ModuleBase::timer::tick(this->classname,"sKG");
    cout<<"Calculating conductivity...."<<endl;
    
    int nw = ceil(wcut/dw_in);
    double dw =  dw_in / ModuleBase::Ry_to_eV; //converge unit in eV to Ry 
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    double dt = ModuleBase::PI/(dw*nw)/times ; //unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    int nt = ceil(sqrt(20)/sigma/dt);
    cout<<"nw: "<<nw<<" ; dw: "<<dw*ModuleBase::Ry_to_eV<<" eV"<<endl;
    cout<<"nt: "<<nt<<" ; dt: "<<dt<<" a.u.(ry^-1)"<<endl;
    assert(nw >= 1);
    assert(nt >= 1);

    ModuleBase::Chebyshev<double> che(this->nche_sto);
    ModuleBase::Chebyshev<double> chet(nche_KG);
    const int npwx = GlobalC::wf.npwx;
    const double tpiba = GlobalC::ucell.tpiba;
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    Stochastic_hchi& stohchi = stoiter.stohchi;
    const int nk = GlobalC::kv.nks;
    const int ndim = 3;

    //------------------------------
    //      Convergence test
    //------------------------------
    bool change = false;
    for (int ik = 0;ik < nk; ++ik)
	{
        int ntest = 2;
        if (this->stowf.nchip[ik] < ntest) 
	    {
	    	ntest = this->stowf.nchip[ik];
	    }

        this->phami->updateHk(ik);
        stoiter.stohchi.current_ik = ik;
        const int npw = GlobalC::kv.ngk[ik];

        complex<double> *pchi = new complex<double> [npw*4*ntest];
        //\sqrt{f}|\chi>
        for(int i = 0 ; i < ntest ; ++i)
            for(int ig = 0 ; ig < npw ; ++ig)
                pchi[i*npw+ig] = this->stowf.shchi[ik](i,ig);
        //(G+k)\sqrt{f}|\chi>
        for(int id = 0 ; id < ndim ; ++id)
        {
            for(int i = 0 ; i < ntest ; ++i)
            {
                for(int ig = 0 ; ig < npw ; ++ig)
                {
                    ModuleBase::Vector3<double> v3 = GlobalC::wfcpw->getgpluskcar(ik,i);
                    pchi[((id+1)*ntest + i)*npw + ig] = this->stowf.shchi[ik](i,ig)*v3[id] * tpiba;
                }
            }
        }
        for(int i = 0; i < ntest * 4 ; ++i)
        {
            while(1)
            {
                bool converge;
                if(this->nche_sto > nche_KG)
                {
                    converge= che.checkconverge(&stohchi, &Stochastic_hchi::hchi_reciprocal, 
	        	    	&pchi[i*npw], npw, stohchi.Emax, stohchi.Emin, 5.0);
                }
                else
                {
                    converge= chet.checkconverge(&stohchi, &Stochastic_hchi::hchi_reciprocal, 
	        	    	&pchi[i*npw], npw, stohchi.Emax, stohchi.Emin, 5.0);
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
        delete[] pchi;
    }

    //------------------------------------------------------------------
    //                    Calculate
    //------------------------------------------------------------------
    const double mu = GlobalC::en.ef;
    stoiter.stofunc.mu = mu;
    double * ct11 = new double[nt];
    double * ct12 = new double[nt];
    double * ct22 = new double[nt];
    ModuleBase::GlobalFunc::ZEROS(ct11,nt);
    ModuleBase::GlobalFunc::ZEROS(ct12,nt);
    ModuleBase::GlobalFunc::ZEROS(ct22,nt);
    double *gxyz = new double[npwx*ndim];

    stoiter.stofunc.t = dt;
    chet.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::n_sin);
    cout<<"Chebyshev precision: "<<abs(chet.coef_complex[nche_KG-1]/chet.coef_complex[0])<<endl;
    ofstream cheofs("Chebycoef");
    for(int i  = 0 ; i < nche_KG ; ++i)
    {
        cheofs<<setw(5)<<i<<setw(20)<<abs(chet.coef_complex[i]/chet.coef_complex[0])<<endl;
    }
    cheofs.close();
    
    //parallel KS orbitals
    int ksbandper,startband;
    parallelks(ksbandper,startband);
    ModuleBase::ComplexMatrix kswf(ksbandper,npwx,false);
    double *ksen = new double[ksbandper];   

    ModuleBase::timer::tick(this->classname,"kloop");
    for (int ik = 0;ik < nk;++ik)
	{
        if(nk > 1) 
        {
            this->phami->updateHk(ik);
        }
        stoiter.stohchi.current_ik = ik;
        const int npw = GlobalC::kv.ngk[ik];

        int nchip = this->stowf.nchip[ik];
        int totbands = nchip + ksbandper;
        complex<double> *leftv3 = new complex<double>[npwx*totbands*ndim];
        complex<double> *jleftv3 = new complex<double>[npwx*totbands*ndim];
        complex<double> *hleftv = new complex<double>[npwx*nchip];
        complex<double> *prightv3 = new complex<double>[npwx*totbands*ndim];
        complex<double> *jrightv3 = new complex<double>[npwx*totbands*ndim];
        complex<double> *rightv = new complex<double>[npwx*totbands];
        

        for(int ib = 0; ib < ksbandper ; ++ib)
        {
            for(int ig = 0 ; ig < npw ; ++ig)
            {
                kswf(ib,ig) = this->psi[0](ik,ib+startband,ig);
            }
            ksen[ib] = this->pelec->ekb(ik, ib+startband);
        }
        for (int ig = 0;ig < npw;ig++)
        {
            ModuleBase::Vector3<double> v3 = GlobalC::wfcpw->getgpluskcar(ik,ig);
            gxyz[ig] = v3.x * tpiba;
            gxyz[ig+npwx] = v3.y * tpiba;
            gxyz[ig+2*npwx] = v3.z * tpiba;
        }
        //-----------------------------------------------------------
        //               sto conductivity
        //-----------------------------------------------------------
        //j|left> = [(Hp+pH)/2 - mu*p]|left>
        this->phami->hPsi(this->stowf.shchi[ik].c, hleftv, nchip*npwx); //hleftv is id transitional array to get jleftv3
        for(int id = 0 ; id < ndim ; ++id)
        {
            for(int ib = 0; ib < nchip ; ++ib) // sto orbitals
		    {
		    	for (int ig = 0; ig < npw; ++ig)
		    	{
		    		leftv3[id*totbands*npwx + ib*npwx + ig] = gxyz[id*npwx + ig] * this->stowf.shchi[ik](ib,ig);
		    	}
		    }
            for(int ib = nchip ; ib < totbands ; ++ib) // KS orbital
            {
                for (int ig = 0; ig < npw; ++ig)
		    	{
		    		leftv3[id*totbands*npwx + ib*npwx + ig] = gxyz[id*npwx + ig] * kswf(ib-nchip,ig);
		    	}
            }
        }
        this->phami->hPsi(leftv3, jleftv3, totbands*ndim*npwx); //leftv3 is a transitional array to get jleftv3
        for(int id = 0 ; id < ndim ; ++id)
        {
            for (int ib = 0; ib < nchip ; ++ib) // sto orbitals
		    {
		    	for (int ig = 0; ig < npw; ++ig)
		    	{
                    int iaibig = id*totbands*npwx + ib*npwx + ig;
                    int ibig = ib *npwx + ig;
                    double ga = gxyz[id*npwx + ig];
                    jleftv3[iaibig] = (jleftv3[iaibig] + hleftv[ibig] * ga ) / 2.0 - mu * this->stowf.shchi[ik](ib,ig) * ga;
                }
            }
            for (int ib = nchip; ib < totbands ; ++ib) //KS orbitals
		    {
		    	for (int ig = 0; ig < npw; ++ig)
		    	{
                    int iaibig = id*totbands*npwx + ib*npwx + ig; 
                    double ga = gxyz[id*npwx + ig];
                    complex<double> wfikib = kswf(ib-nchip,ig);
                    jleftv3[iaibig] = (jleftv3[iaibig] + ga * ksen[ib-nchip] * wfikib ) / 2.0 - mu * wfikib * ga;
                }
            }
        }

        //p|left> Note: p^+=p, px=-id/dx  p|k>=k|k>
        for(int id = 0 ; id < ndim ; ++id)
        {
            for (int ib = 0; ib < nchip ; ++ib) // sto orbitals
		    {
		    	for (int ig = 0; ig < npw; ++ig)
		    	{
		    		leftv3[id*totbands*npwx + ib*npwx + ig] = gxyz[id*npwx + ig] * this->stowf.shchi[ik](ib,ig);
		    	}
		    }
            for(int ib = nchip; ib < totbands; ++ib) // KS orbitals
            {
                for (int ig = 0; ig < npw; ++ig)
		    	{
		    		leftv3[id*totbands*npwx + ib*npwx + ig] = gxyz[id*npwx + ig] * kswf(ib-nchip,ig);
		    	}
            }
        }

        ModuleBase::timer::tick(this->classname,"chebyshev");
        //(1-f)|left>
        che.calcoef_real(&stoiter.stofunc,&Sto_Func<double>::n_fd);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_reciprocal, leftv3, leftv3, npw, npwx, totbands*ndim);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_reciprocal, jleftv3, jleftv3, npw, npwx, totbands*ndim);
        ModuleBase::timer::tick(this->classname,"chebyshev");

                
        cout<<"ik="<<ik<<": ";
        for (int it = 1 ;it < nt ; ++it)
        {
            if(it%20==0) cout<<it<<" ";
            ModuleBase::timer::tick(this->classname,"chebyshev");
            //exp(-iH*dt/h)|left> Note: exp(iH*dt/h)^+=exp(-iHt/h)   
            chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, leftv3, leftv3, npw, npwx, totbands*ndim);
            chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, jleftv3, jleftv3, npw, npwx, totbands*ndim);


            //exp(-iH*dt/h)|right>
            if(it == 1) chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, this->stowf.shchi[ik].c, rightv, npw, npwx, nchip);
            else        chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, rightv, rightv, npw, npwx, nchip);
            ModuleBase::timer::tick(this->classname,"chebyshev");
            for(int ib = nchip ; ib < totbands ; ++ib)
            {
                for (int ig = 0; ig < npw; ++ig)
		    	{
		    		rightv[ib*npwx + ig] = exp(-ModuleBase::IMAG_UNIT * ksen[ib-nchip] * double(it) * dt) * kswf(ib-nchip,ig);
		    	}
            }

            //-i[(Hp+pH)/2 - mu*p]|right>
            this->phami->hPsi(rightv, hleftv, nchip*npwx); //hleftv is a transitional array to get jleftv3
            for(int id = 0 ; id < ndim ; ++id)
            {
                for (int ib = 0; ib < totbands ; ++ib) //KS + sto
		        {
		        	for (int ig = 0; ig < npw; ++ig)
		        	{
		        		prightv3[id*totbands*npwx + ib*npwx + ig] = gxyz[id*npwx + ig] * rightv[ib*npwx + ig];
		        	}
		        }
            }
            this->phami->hPsi(prightv3, jrightv3, totbands*ndim*npwx); //prightv3 is a transitional array to get jleftv3
            for(int id = 0 ; id < ndim ; ++id)
            {
                for (int ib = 0; ib < nchip ; ++ib) //sto orbitals
		        {
		        	for (int ig = 0; ig < npw; ++ig)
		        	{
                        int iaibig = id*totbands*npwx + ib*npwx + ig;
                        int ibig = ib *npwx + ig;
                        double ga = gxyz[id*npwx + ig];
                        jrightv3[iaibig] = -ModuleBase::IMAG_UNIT * ((jrightv3[iaibig] + hleftv[ibig] * ga ) / double(2) - mu * rightv[ib*npwx + ig] * ga);
                    }
                }
                for (int ib = nchip; ib < totbands ; ++ib) //KS orbitals
		        {
		        	for (int ig = 0; ig < npw; ++ig)
		        	{
                        int iaibig = id*totbands*npwx + ib*npwx + ig; 
                        double ga = gxyz[id*npwx + ig];
                        complex<double> wfikib = rightv[ib*npwx + ig];
                        jrightv3[iaibig] = -ModuleBase::IMAG_UNIT * ((jrightv3[iaibig] + ga * ksen[ib-nchip] * wfikib ) / double(2) - mu * wfikib * ga);
                    }
                }
            }

            //-ip|right>
            for(int id = 0 ; id < ndim ; ++id)
            {      
                for (int ib = 0; ib < totbands ; ++ib) //sto + KS
		        {
		    	    for (int ig = 0; ig < npw; ++ig)
		    	    {
		    		    prightv3[id*totbands*npwx + ib*npwx + ig] = -ModuleBase::IMAG_UNIT * gxyz[id*npwx + ig] * rightv[ib*npwx + ig];
		    	    }
		        }
            }

            //Im<left|p|right>=Re<left|-ip|right>
            for(int id = 0 ; id < ndim ; ++id)
            {
                for(int ib = 0 ; ib < nchip ; ++ib)
                {
                    ct11[it] += ModuleBase::GlobalFunc::ddot_real(npw,&leftv3[id*totbands*npwx + ib*npwx],&prightv3[id*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2,0;
                    ct12[it] -= ModuleBase::GlobalFunc::ddot_real(npw,&jleftv3[id*totbands*npwx + ib*npwx],&prightv3[id*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0;
                    ct22[it] += ModuleBase::GlobalFunc::ddot_real(npw,&jleftv3[id*totbands*npwx + ib*npwx],&jrightv3[id*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0;
                }
                for(int ib = nchip ; ib < totbands ; ++ib)
                {
                    double ei = ksen[ib-nchip];
                    double fi = stoiter.stofunc.fd(ei);
                    ct11[it] += ModuleBase::GlobalFunc::ddot_real(npw,&leftv3[id*totbands*npwx + ib*npwx],&prightv3[id*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0 * fi;
                    ct12[it] -= ModuleBase::GlobalFunc::ddot_real(npw,&jleftv3[id*totbands*npwx + ib*npwx],&prightv3[id*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0 * fi;
                    ct22[it] += ModuleBase::GlobalFunc::ddot_real(npw,&jleftv3[id*totbands*npwx + ib*npwx],&jrightv3[id*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0 * fi;
                }
            }
        }
        cout<<endl;

        delete [] leftv3;
        delete [] jleftv3;
        delete [] hleftv;
        delete [] prightv3;
        delete [] jrightv3;
        delete [] rightv;
    }
    ModuleBase::timer::tick(this->classname,"kloop");
    delete[] gxyz;
    delete[] ksen;
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
    // delete[] cterror;
    ModuleBase::timer::tick(this->classname,"sKG");
}

void ESolver_SDFT_PW::sKG_new(const int nche_KG, const double fwhmin, const double wcut, 
                          const double dw_in, const int times)
{
    ModuleBase::timer::tick(this->classname,"sKG");
    cout<<"Calculating conductivity...."<<endl;
    
    int nw = ceil(wcut/dw_in);
    double dw =  dw_in / ModuleBase::Ry_to_eV; //converge unit in eV to Ry 
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    double dt = ModuleBase::PI/(dw*nw)/times ; //unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    int nt = ceil(sqrt(20)/sigma/dt);
    cout<<"nw: "<<nw<<" ; dw: "<<dw*ModuleBase::Ry_to_eV<<" eV"<<endl;
    cout<<"nt: "<<nt<<" ; dt: "<<dt<<" a.u.(ry^-1)"<<endl;
    assert(nw >= 1);
    assert(nt >= 1);
    const int ndim = 3;

    ModuleBase::Chebyshev<double> che(this->nche_sto);
    ModuleBase::Chebyshev<double> chet(nche_KG);
    ModuleBase::Chebyshev<double> chet2(nche_KG);
    const int npwx = GlobalC::wf.npwx;
    const double tpiba = GlobalC::ucell.tpiba;
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    Stochastic_hchi& stohchi = stoiter.stohchi;
    const int nk = GlobalC::kv.nks;

    //------------------------------
    //      Convergence test
    //------------------------------
    bool change = false;
    for (int ik = 0;ik < nk; ++ik)
	{
        int ntest = 2;
        if (this->stowf.nchip[ik] < ntest) 
	    {
	    	ntest = this->stowf.nchip[ik];
	    }

        this->phami->updateHk(ik);
        stoiter.stohchi.current_ik = ik;
        const int npw = GlobalC::kv.ngk[ik];

        for(int j = 0 ; j < 2 ; ++j)
        {
        for(int i = 0; i < ntest ; ++i)
        {
            std::complex<double> *pchi;
            if(j==0)
            {
                if(GlobalV::NBANDS > 0)
                    pchi = &stowf.chiortho[ik](i,0);
                else
                    pchi = &stowf.chi0[ik](i,0);
            }
            else
            {
                if(this->nche_sto > nche_KG) break;
                pchi = &stowf.shchi[ik](i,0);
            }
            
            while(1)
            {
                bool converge;
                converge= chet.checkconverge(&stohchi, &Stochastic_hchi::hchi_reciprocal, 
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
        }

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

    //------------------------------------------------------------------
    //                    Calculate
    //------------------------------------------------------------------
    const double mu = GlobalC::en.ef;
    stoiter.stofunc.mu = mu;
    double * ct11 = new double[nt];
    double * ct12 = new double[nt];
    double * ct22 = new double[nt];
    ModuleBase::GlobalFunc::ZEROS(ct11,nt);
    ModuleBase::GlobalFunc::ZEROS(ct12,nt);
    ModuleBase::GlobalFunc::ZEROS(ct22,nt);

    stoiter.stofunc.t = dt;
    chet.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::nsin);
    chet2.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::n_sin);
    cout<<"Chebyshev precision: "<<abs(chet.coef_complex[nche_KG-1]/chet.coef_complex[0])<<endl;
    ofstream cheofs("Chebycoef");
    for(int i  = 0 ; i < nche_KG ; ++i)
    {
        cheofs<<setw(5)<<i<<setw(20)<<abs(chet.coef_complex[i]/chet.coef_complex[0])<<endl;
    }
    cheofs.close();

    
    //parallel KS orbitals
    int ksbandper,startband;
    parallelks(ksbandper,startband);
    ModuleBase::timer::tick(this->classname,"kloop");
    for (int ik = 0;ik < nk;++ik)
	{
        if(nk > 1) 
        {
            this->phami->updateHk(ik);
        }
        stoiter.stohchi.current_ik = ik;
        const int npw = GlobalC::kv.ngk[ik];

        int nchip = this->stowf.nchip[ik];
        int totbands_per = nchip + ksbandper;
        int totbands = totbands_per;
#ifdef __MPI
        MPI_Allreduce(&totbands_per,&totbands,1,MPI_INT,MPI_SUM,PARAPW_WORLD);
        const int nstogroup = GlobalV::NSTOGROUP;
        int nrecv[nstogroup];
        int displs[nstogroup];
        MPI_Allgather(&totbands_per,1,MPI_INT,nrecv,1,MPI_INT,PARAPW_WORLD);
        displs[0] = 0;
        for(int i = 1; i < nstogroup ; ++i)
        {
            displs[i] = displs[i-1] + nrecv[i-1];
        }
        for(int i = 0 ; i < nstogroup ; ++i)
        {
            nrecv[i] *= npwx;
            displs[i] *= npwx;
        }
#endif
        

        //-----------------------------------------------------------
        //               sto conductivity
        //-----------------------------------------------------------
        //before loop

        //|psi>
        psi::Psi<std::complex<double>> psi0(1,totbands_per,npwx); //|psi>
        psi::Psi<std::complex<double>> sfpsi0(1,totbands_per,npwx); //sqrt(f)|psi>
        psi::Psi<std::complex<double>> hpsi0(1,totbands_per,npwx); //h|psi>
        psi::Psi<std::complex<double>> hsfpsi0(1,totbands_per,npwx); //h*sqrt(f)|psi>
        //j|psi> j1=p  j2=(Hp+pH)/2 - mu*p
        psi::Psi<std::complex<double>> j1psi(ndim,totbands_per,npwx);
        psi::Psi<std::complex<double>> j2psi(ndim,totbands_per,npwx);
        //(1-f)*j*sqrt(f)|psi>
        psi::Psi<std::complex<double>> j1sfpsi(ndim,totbands_per,npwx);
        psi::Psi<std::complex<double>> j2sfpsi(ndim,totbands_per,npwx);
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
                if(GlobalV::NBANDS > 0)
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
        for(int id = 0 ; id < ndim ; ++id)
        {
            for(int ib = 0 ; ib < totbands_per ; ++ib )
            {
                for(int ig = 0; ig < npw ; ++ig)
                {
                    std::complex<double> psiig = psi0(ib,ig);
                    std::complex<double> sfpsiig = sfpsi0(ib,ig);
                    double gplusk_a = GlobalC::wfcpw->getgpluskcar(ik,ig)[id];
                    j1psi(id,ib,ig) = tpiba * gplusk_a * psiig;
                    j1sfpsi(id,ib,ig) = tpiba * gplusk_a * sfpsiig;
                }
            }
        }
        this->phami->hPsi(psi0.get_pointer(), hpsi0.get_pointer(), totbands_per*npwx);
        this->phami->hPsi(sfpsi0.get_pointer(), hsfpsi0.get_pointer(), totbands_per*npwx);
        this->phami->hPsi(j1psi.get_pointer(), j2psi.get_pointer(), ndim*totbands_per*npwx);
        this->phami->hPsi(j1sfpsi.get_pointer(), j2sfpsi.get_pointer(), ndim*totbands_per*npwx);
        GlobalV::ofs_running<<"11111111"<<endl;
        GlobalV::ofs_running<<j2psi(0,0,3)<<" "<<j2psi(1,0,3)<<" "<<j2psi(2,0,3)<<endl;
        GlobalV::ofs_running<<j2sfpsi(0,0,3)<<" "<<j2sfpsi(1,0,3)<<" "<<j2sfpsi(2,0,3)<<endl;
        for(int id = 0 ; id < ndim ; ++id)
        {
            for(int ib = 0 ; ib < totbands_per ; ++ib )
            {
                for(int ig = 0; ig < npw ; ++ig)
                {
                    std::complex<double> hpsiig = hpsi0(ib,ig);
                    std::complex<double> hsfpsiig = hsfpsi0(ib,ig);
                    double gplusk_a = GlobalC::wfcpw->getgpluskcar(ik,ig)[id];
                    j2psi(id,ib,ig) = (j2psi(id,ib,ig) + tpiba * gplusk_a * hpsiig)/2.0 - mu * j1psi(id,ib,ig);
                    j2sfpsi(id,ib,ig) = (j2sfpsi(id,ib,ig) + tpiba * gplusk_a * hsfpsiig)/2.0 - mu * j1sfpsi(id,ib,ig);
                }
            }
        }
        GlobalV::ofs_running<<"222222"<<endl;
        GlobalV::ofs_running<<j2psi(0,0,3)<<" "<<j2psi(1,0,3)<<" "<<j2psi(2,0,3)<<endl;
        GlobalV::ofs_running<<j2sfpsi(0,0,3)<<" "<<j2sfpsi(1,0,3)<<" "<<j2sfpsi(2,0,3)<<endl;

        //(1-f)
        che.calcoef_real(&stoiter.stofunc,&Sto_Func<double>::n_fd);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_reciprocal, j1sfpsi.get_pointer(), j1sfpsi.get_pointer(), npw, npwx, totbands_per*ndim);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_reciprocal, j2sfpsi.get_pointer(), j2sfpsi.get_pointer(), npw, npwx, totbands_per*ndim);
        GlobalV::ofs_running<<"33333333"<<endl;
        GlobalV::ofs_running<<j2psi(0,0,3)<<" "<<j2psi(1,0,3)<<" "<<j2psi(2,0,3)<<endl;
        GlobalV::ofs_running<<j2sfpsi(0,0,3)<<" "<<j2sfpsi(1,0,3)<<" "<<j2sfpsi(2,0,3)<<endl;
        psi::Psi<std::complex<double>> *p_j1psi = &j1psi; 
        psi::Psi<std::complex<double>> *p_j2psi = &j2psi; 
#ifdef __MPI
        psi::Psi<std::complex<double>> j1psi_tot,j2psi_tot;
        if (GlobalV::NSTOGROUP > 1)
        {
            j1psi_tot.resize(ndim,totbands,npwx);
            j2psi_tot.resize(ndim,totbands,npwx);
            for(int id = 0 ; id < ndim ; ++id)
            {
                MPI_Allgatherv(&j1psi(id,0,0), totbands_per * npwx, mpicomplex, 
                            &j1psi_tot(id,0,0), nrecv, displs, mpicomplex, PARAPW_WORLD);
                MPI_Allgatherv(&j2psi(id,0,0), totbands_per * npwx, mpicomplex, 
                            &j2psi_tot(id,0,0), nrecv, displs, mpicomplex, PARAPW_WORLD);
            }
            p_j1psi = &j1psi_tot;
            p_j2psi = &j2psi_tot;
        }
#endif


        //loop of t
        psi::Psi<std::complex<double>> exppsi(1,totbands_per,npwx);
        psi::Psi<std::complex<double>> expsfpsi(1,totbands_per,npwx);
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
            
            //exp(iHdt)|chi>
            chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, &exppsi(ksbandper,0), &exppsi(ksbandper,0), npw, npwx, nchip);
            //exp(-iHdt)|shchi>
            chet2.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, &expsfpsi(ksbandper,0), &expsfpsi(ksbandper,0), npw, npwx, nchip);
            psi::Psi<std::complex<double>> *p_exppsi = &exppsi;
#ifdef __MPI
            psi::Psi<std::complex<double>> exppsi_tot;
            if (GlobalV::NSTOGROUP > 1)
            {
                exppsi_tot.resize(1,totbands,npwx);
                MPI_Allgatherv(&exppsi(0,0), totbands_per * npwx, mpicomplex, 
                                &exppsi_tot(0,0), nrecv, displs, mpicomplex, PARAPW_WORLD);
                p_exppsi = &exppsi_tot;
            }
#endif
            ModuleBase::ComplexMatrix j1l(ndim,totbands_per*totbands), j2l(ndim,totbands_per*totbands);
            ModuleBase::ComplexMatrix j1r(ndim,totbands_per*totbands), j2r(ndim,totbands_per*totbands);
            char transa = 'C';
            char transb = 'N';
            int totbands_per3 = ndim*totbands_per;
            int totbands3 = ndim*totbands;
            for(int id = 0 ; id < ndim ; ++id)
            {
                //<psi|sqrt(f)j_1(1-f) exp(iHt)|psi>
                zgemm_(&transa, &transb,&totbands_per, &totbands, &npw, &ModuleBase::ONE, &j1sfpsi(id,0,0), &npwx, 
                        p_exppsi->get_pointer(), &npwx, &ModuleBase::ZERO, &j1l(id,0), &totbands_per);
                //<psi|sqrt(f)j_2(1-f) exp(iHt)|psi>
                zgemm_(&transa, &transb,&totbands_per, &totbands, &npw, &ModuleBase::ONE, &j2sfpsi(id,0,0), &npwx, 
                        p_exppsi->get_pointer(), &npwx, &ModuleBase::ZERO, &j2l(id,0), &totbands_per);
            }
            for(int id = 0; id < ndim ; ++id) // it can also use zgemm once
            {
                //i<exp(-iHt)sqrt(f)psi| j_1|psi> = i<psi|sqrt(f)exp(iHt) j_1|psi> = i(<psi|j_1 exp(-iHt)sqrt(f)|psi>)^+
                zgemm_(&transa, &transb,&totbands_per, &totbands, &npw, &ModuleBase::IMAG_UNIT, expsfpsi.get_pointer(), &npwx,
                        &(p_j1psi->operator()(id,0,0)), &npwx, &ModuleBase::ZERO, &j1r(id,0), &totbands_per);
                //i<psi|sqrt(f)exp(iHt) j_2|psi> 
                zgemm_(&transa, &transb,&totbands_per, &totbands, &npw, &ModuleBase::IMAG_UNIT, expsfpsi.get_pointer(), &npwx,
                        &(p_j2psi->operator()(id,0,0)), &npwx, &ModuleBase::ZERO, &j2r(id,0), &totbands_per);
            }
            GlobalV::ofs_running<<"11111111"<<endl;
            GlobalV::ofs_running<<j2l(0,23)<<" "<<j2l(1,23)<<" "<<j2l(2,23)<<endl;
            GlobalV::ofs_running<<j2r(0,23)<<" "<<j2r(1,23)<<" "<<j2r(2,23)<<endl;

#ifdef __MPI
            MPI_Allreduce(MPI_IN_PLACE,j1l.c,ndim*totbands_per*totbands*2,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,j2l.c,ndim*totbands_per*totbands*2,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,j1r.c,ndim*totbands_per*totbands*2,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,j2r.c,ndim*totbands_per*totbands*2,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#endif
            GlobalV::ofs_running<<"2222222"<<endl;
            GlobalV::ofs_running<<j2l(0,23)<<" "<<j2l(1,23)<<" "<<j2l(2,23)<<endl;
            GlobalV::ofs_running<<j2r(0,23)<<" "<<j2r(1,23)<<" "<<j2r(2,23)<<endl;
            //Re(i<psi|sqrt(f)j(1-f) exp(iHt)|psi><psi|j exp(-iHt)\sqrt(f)|psi>)
            if(GlobalV::RANK_IN_POOL==0)
            {
                //Im(l_ij*r_ji)=Re(i l^*_ij*r^+_ij)=Re(i l^*_i*r^+_i)
                //ddot_real = real(A^*_i * B_i)
                ct11[it] += ModuleBase::GlobalFunc::ddot_real(totbands*totbands_per*ndim,j1l.c,j1r.c,false) * GlobalC::kv.wk[ik] / 2,0;
                ct12[it] -= ModuleBase::GlobalFunc::ddot_real(totbands*totbands_per*ndim,j1l.c,j2r.c,false) * GlobalC::kv.wk[ik] / 2,0;
                ct22[it] += ModuleBase::GlobalFunc::ddot_real(totbands*totbands_per*ndim,j2l.c,j2r.c,false) * GlobalC::kv.wk[ik] / 2,0;
            }
        }
        cout<<endl;
        if(ksbandper > 0)   delete[] en;

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
    // delete[] cterror;
    ModuleBase::timer::tick(this->classname,"sKG");
}

}