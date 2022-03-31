//=========================================================
//AUTHOR : Peize Lin & jingan
//DATE   : 2019-01-16
//UPDATE : 2019-06-28 
//=========================================================

#include "../module_base/blas_connector.h"
#include "../module_base/scalapack_connector.h"
#include "../module_base/timer.h"
#include "global_fp.h"
#include "../src_pw/global.h"

#include "../src_external/src_test/test_function.h"
#include "../src_external/src_test/src_global/complexmatrix-test.h"

#include "src_lcao/local_orbital_charge.h"

void Local_Orbital_Charge::init_dm_2d()
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->dm_gamma.resize(GlobalV::NSPIN);
	}
	else
	{
		this->dm_k.resize(GlobalC::kv.nks);
	}
}


//for gamma_only
void Local_Orbital_Charge::cal_dm(const ModuleBase::matrix& wg,
    std::vector<ModuleBase::matrix>& wfc_gamma,
    std::vector<ModuleBase::matrix>& dm_gamma)
{
	ModuleBase::TITLE("Local_Orbital_Charge", "cal_dm");
	
	#ifdef TEST_DIAG
    static int istep=0;
    std::ofstream ofs("wfc_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
    ofs<<wfc_gamma<<std::endl;
	#endif
	
	// dm = wfc.T * wg * wfc.conj()
	// dm[is](iw1,iw2) = \sum_{ib} wfc[is](ib,iw1).T * wg(is,ib) * wfc[is](ib,iw2).conj()
	assert(wg.nc<=GlobalV::NLOCAL);
    assert(wg.nr==GlobalV::NSPIN);
    for(int is=0; is!=GlobalV::NSPIN; ++is)
    {
        std::vector<double> wg_local(this->ParaV->ncol,0.0);
        for(int ib_global=0; ib_global!=wg.nc; ++ib_global)
        {
            const int ib_local = this->ParaV->trace_loc_col[ib_global];
            if(ib_local>=0)
            {
                wg_local[ib_local] = wg(is,ib_global);
            }
        }
        
        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        ModuleBase::matrix wg_wfc(wfc_gamma[is]);
        for(int ir=0; ir!=wg_wfc.nr; ++ir)
        {
            BlasConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
        const double one_float=1.0, zero_float=0.0;
        const int one_int=1;
        const char N_char='N', T_char='T';
        dm_gamma[is].create( this->ParaV->ncol, this->ParaV->nrow );
    #ifdef __MPI
        pdgemm_(
            &N_char, &T_char, 
            &GlobalV::NLOCAL, &GlobalV::NLOCAL, &wg.nc,
            &one_float,
            wg_wfc.c, &one_int, &one_int, this->ParaV->desc,
            wfc_gamma[is].c, &one_int, &one_int, this->ParaV->desc,
            &zero_float,
            dm_gamma[is].c, &one_int, &one_int, this->ParaV->desc);
    #else
        const int lda=GlobalV::NLOCAL;
        dgemm_(
            &N_char, &T_char, 
            &GlobalV::NLOCAL, &GlobalV::NLOCAL, &wg.nc,
            &one_float,
            wg_wfc.c, &lda,
            wfc_gamma[is].c, &lda,
            &zero_float,
            dm_gamma[is].c, &lda);
    #endif    
    }
    #ifdef TEST_DIAG
    static int istep=0;
    std::ofstream ofs("dm_" + ModuleBase::GlobalFunc::TO_STRING(istep) + "_" + ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
    ofs << dm_gamma << std::endl;
	#endif

	return;
}

//for multi-k
void Local_Orbital_Charge::cal_dm(const ModuleBase::matrix& wg,    // wg(ik,ib), cal all dm 
    std::vector<ModuleBase::ComplexMatrix>& wfc_k,
    std::vector<ModuleBase::ComplexMatrix> &dm_k)
{
	ModuleBase::TITLE("Local_Orbital_Charge", "cal_dm");
	
	#ifdef TEST_DIAG
    static int istep=0;
    std::ofstream ofs("wfc_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
    ofs<<wfc_k<<std::endl;
	#endif
	//
	// dm = wfc.T * wg * wfc.conj()
	// dm[ik](iw1,iw2) = \sum_{ib} wfc[ik](ib,iw1).T * wg(ik,ib) * wfc[ik](ib,iw2).conj()
	assert(wg.nc<=GlobalV::NLOCAL);
    assert(wg.nr==GlobalC::kv.nks);
    for(int ik=0; ik!=GlobalC::kv.nks; ++ik)
    {
        std::vector<double> wg_local(this->ParaV->ncol,0.0);
        for(int ib_global=0; ib_global!=wg.nc; ++ib_global)
        {
            const int ib_local = this->ParaV->trace_loc_col[ib_global];
            if(ib_local>=0)
            {
                wg_local[ib_local] = wg(ik,ib_global);
            }
        }

        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw).conj();
        ModuleBase::ComplexMatrix wg_wfc = conj(wfc_k[ik]);
        for(int ir=0; ir!=wg_wfc.nr; ++ir)
        {
            BlasConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );
        }
        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
        const double one_float[2]={1.0, 0.0}, zero_float[2]={0.0, 0.0};
        const int one_int=1;
        const char N_char='N', T_char='T';
        dm_k[ik].create( this->ParaV->ncol, this->ParaV->nrow );
#ifdef __MPI
        pzgemm_(
            &N_char, &T_char,
            &GlobalV::NLOCAL, &GlobalV::NLOCAL, &wg.nc,
            &one_float[0],
            wg_wfc.c, &one_int, &one_int, this->ParaV->desc_wfc,
            wfc_k[ik].c, &one_int, &one_int, this->ParaV->desc_wfc,
            &zero_float[0],
            dm_k[ik].c, &one_int, &one_int, this->ParaV->desc);
    #else
        const int lda=GlobalV::NLOCAL;
        const complex<double> one_complex={1.0,0.0}, zero_complex={0.0,0.0};
        zgemm_(
            &N_char, &T_char, 
            &GlobalV::NLOCAL, &GlobalV::NLOCAL, &wg.nc,
            &one_complex,
            wg_wfc.c, &lda,
            wfc_k[ik].c, &lda,
            &zero_complex,
            dm_k[ik].c, &lda);
    #endif
    }
	
	#ifdef TEST_DIAG
    static int istep=0;
    std::ofstream ofs("dm_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
    ofs<<dm_k<<std::endl;
	#endif

	return;
}

//must cal cal_dm first
void Local_Orbital_Charge::cal_dm_R(
    std::vector<ModuleBase::ComplexMatrix> &dm_k,
    Record_adj& ra,    //ra.for_2d();
    double** dm2d)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "cal_dm_R");
    assert(dm_k[0].nr > 0 && dm_k[0].nc > 0); //must call cal_dm first

    for (int ik = 0;ik < GlobalC::kv.nks;++ik)
    {
        // allocate memory and pointer for each ispin
        int ispin = 0;
        if (GlobalV::NSPIN == 2)
        {
            ispin = GlobalC::kv.isk[ik];
        }
        for (int T1 = 0;T1 < GlobalC::ucell.ntype;++T1)
        {
            for (int I1 = 0;I1 < GlobalC::ucell.atoms[T1].na;++I1)
            {
                const int iat = GlobalC::ucell.itia2iat(T1, I1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                //irr: number of adjacent orbital pairs int this proc
                const int irrstart = this->ParaV->nlocstart[iat];

                int count = 0;
                for (int cb = 0;cb < ra.na_each[iat];++cb)
                {
                    const int T2 = ra.info[iat][cb][3];
                    const int I2 = ra.info[iat][cb][4];
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                    //-----------------
                    // exp[i * R * k]
                    //-----------------
                    const std::complex<double> phase =
                        exp(ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (
                            GlobalC::kv.kvec_d[ik].x * ra.info[iat][cb][0] +
                            GlobalC::kv.kvec_d[ik].y * ra.info[iat][cb][1] +
                            GlobalC::kv.kvec_d[ik].z * ra.info[iat][cb][2]
                            ));
                    for (int iw1 = 0;iw1 < GlobalC::ucell.atoms[T1].nw;++iw1)
                    {
                        int iw1_all = start1 + iw1;
                        int mu = this->ParaV->trace_loc_row[iw1_all];
                        if (mu < 0)continue;
                        for (int iw2 = 0;iw2 < GlobalC::ucell.atoms[T2].nw;++iw2)
                        {
                            int iw2_all = start2 + iw2;
                            int nu = this->ParaV->trace_loc_col[iw2_all];
                            if (nu < 0)continue;
                            //Caution: output of pzgemm_ : col first in **each** proc itself !!
                            dm2d[ispin][irrstart + count] += (dm_k[ik](nu, mu) * phase).real();
                            ++count;
                        }//iw2
                    }//iw1
                }//TI2(cb)
                assert(count == this->ParaV->nlocdim[iat]);
            }//I1
        }//T1
    }//ik
    ModuleBase::timer::tick("Local_Orbital_Charge", "cal_dm_R");
    return;
}