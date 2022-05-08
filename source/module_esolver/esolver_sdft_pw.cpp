#include "./esolver_sdft_pw.h"

namespace ModuleESolver
{

ESolver_SDFT_PW::ESolver_SDFT_PW()
{
    classname = "ESolver_SDFT_PW";
    basisname = "PW";
}

void ESolver_SDFT_PW::Init(Input &inp, UnitCell_pseudo &cell)
{
    // ESolver_KS_PW::Init(inp, cell);
    // STO_WF.alloc();
	// stoiter.alloc( wf.npwx );
}

void ESolver_SDFT_PW::beforescf(const int istep)
{
    ESolver_KS_PW::beforescf(istep);
    // STO_WF.init();
	// stoiter.init( wf.npwx );
	// if(NITER==0)
	// {
	// 	int iter = 1;
	// 	ETHR = 0.1*DRHO2/ std::max(1.0, ucell.nelec);
	// 	double *h_diag = new double[wf.npwx * NPOL];
	// 	for (int ik = 0;ik < kv.nks;ik++)
	// 	{
	// 		if(NBANDS > 0 && MY_STOGROUP == 0)
	// 		{
	// 			this->c_bands_k(ik,h_diag,istep+1);
	// 		}
	// 		else
	// 		{
	// 			hm.hpw.init_k(ik); 
	// 		}
	// 		if(NBANDS > 0)
	// 		{
	// 			MPI_Bcast(wf.evc[ik].c, wf.npwx*NBANDS*2, MPI_DOUBLE , 0, PARAPW_WORLD);
	// 			MPI_Bcast(wf.ekb[ik], NBANDS, MPI_DOUBLE, 0, PARAPW_WORLD);
	// 		}
	// 		stoiter.stoche.ndmin = wf.npw;
	// 		stoiter.orthog(ik);
	// 		stoiter.checkemm(ik,iter);
	// 	}
	// 	for (int ik = 0;ik < kv.nks;ik++)
	// 	{
	// 		//init k
	// 		if(kv.nks > 1) hm.hpw.init_k(ik);
	// 		stoiter.stoche.ndmin = wf.npw;

	// 		stoiter.sumpolyval_k(ik);
	// 	}
	// 	delete [] h_diag;
	// 	stoiter.itermu(iter);
	// 	stoiter.stoche.calcoef(Stochastic_Iter::nroot_fd);
	// 	for(int ik = 0; ik < kv.nks; ++ik)
    // 	{
    // 	    //init k
    // 	    if(kv.nks > 1) hm.hpw.init_k(ik);
    // 		stoiter.stoche.ndmin = wf.npw;
	// 		complex<double> * out, *pchi;
	// 		out = STO_WF.shchi[ik].c;
    //     	if(NBANDS > 0)
    //     	    pchi = STO_WF.chiortho[ik].c;
    //     	else
    //     	    pchi = STO_WF.chi0[ik].c;
	
    //     	stoiter.stoche.calfinalvec(stoiter.stohchi.hchi_reciprocal, pchi, out, stoiter.nchip[ik]);
	// 	}

	// }
}

void ESolver_SDFT_PW::eachiterfinish(int iter, bool conv_elec)
{

}
void ESolver_SDFT_PW::afterscf(const int iter, bool conv_elec)
{

}

void ESolver_SDFT_PW::hamilt2density(int istep, int iter, double ethr)
{
    
}

void ESolver_SDFT_PW::cal_Energy(energy &en)
{

}

void ESolver_SDFT_PW::cal_Force(ModuleBase::matrix &force)
{
	
}
void ESolver_SDFT_PW::cal_Stress(ModuleBase::matrix &stress)
{

}


}