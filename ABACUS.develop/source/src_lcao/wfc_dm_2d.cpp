//=========================================================
//AUTHOR : Peize Lin & jingan
//DATE   : 2019-01-16
//UPDATE : 2019-06-28 
//=========================================================

#include "wfc_dm_2d.h"
#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"

#include "src_external/src_test/test_function.h"
#include "src_external/src_test/src_global/complexmatrix-test.h"

void Wfc_Dm_2d::init()
{
	TITLE("Wfc_Dm_2d", "init");
	if(GAMMA_ONLY_LOCAL)
	{
		wfc_gamma.resize(NSPIN);
		dm_gamma.resize(NSPIN);
	}
	else
	{
		wfc_k.resize(kv.nks);
		dm_k.resize(kv.nks);
	}
}

void Wfc_Dm_2d::cal_dm(const matrix &wg)
{
	TITLE("Wfc_Dm_2d", "cal_dm");
	
	#ifdef TEST_DIAG
	{
		static int istep=0;
		ofstream ofs("wfc_"+TO_STRING(istep++)+"_"+TO_STRING(MY_RANK));
		if(GAMMA_ONLY_LOCAL)
		{
			ofs<<wfc_gamma<<endl;
		}
		else
		{
			ofs<<wfc_k<<endl;
		}
	}
	#endif
	
	// dm = wfc.T * wg * wfc.conj()
	// dm[ik](iw1,iw2) = \sum_{ib} wfc[ik](ib,iw1).T * wg(ik,ib) * wfc[ik](ib,iw2).conj()
	assert(wg.nc<=NLOCAL);
	if(GAMMA_ONLY_LOCAL)
	{
		assert(wg.nr==NSPIN);
		for(int is=0; is!=NSPIN; ++is)
		{
			std::vector<double> wg_local(ParaO.ncol,0.0);
			for(int ib_global=0; ib_global!=wg.nc; ++ib_global)
			{
				const int ib_local = ParaO.trace_loc_col[ib_global];
				if(ib_local>=0)
				{
					wg_local[ib_local] = wg(is,ib_global);
				}
			}
			
			// wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
			matrix wg_wfc(wfc_gamma[is]);
			for(int ir=0; ir!=wg_wfc.nr; ++ir)
			{
				LapackConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );
			}

			// C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
			const double one_float=1.0, zero_float=0.0;
			const int one_int=1;
			const char N_char='N', T_char='T';
			dm_gamma[is].create( wfc_gamma[is].nr, wfc_gamma[is].nc );
			pdgemm_(
				&N_char, &T_char,
				&NLOCAL, &NLOCAL, &wg.nc,
				&one_float,
				wg_wfc.c, &one_int, &one_int, ParaO.desc,
				wfc_gamma[is].c, &one_int, &one_int, ParaO.desc,
				&zero_float,
				dm_gamma[is].c, &one_int, &one_int, ParaO.desc);
		}
	}
	else
	{
		assert(wg.nr==kv.nks);
		for(int ik=0; ik!=kv.nks; ++ik)
		{
			std::vector<double> wg_local(ParaO.ncol,0.0);
			for(int ib_global=0; ib_global!=wg.nc; ++ib_global)
			{
				const int ib_local = ParaO.trace_loc_col[ib_global];
				if(ib_local>=0)
				{
					wg_local[ib_local] = wg(ik,ib_global);
				}
			}

			// wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw).conj();
			ComplexMatrix wg_wfc = conj(wfc_k[ik]);
			for(int ir=0; ir!=wg_wfc.nr; ++ir)
			{
				LapackConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );
			}

			// C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
			const double one_float=1.0, zero_float=0.0;
			const int one_int=1;
			const char N_char='N', T_char='T';
			dm_k[ik].create( wfc_k[ik].nr, wfc_k[ik].nc );
			pzgemm_(
				&N_char, &T_char,
				&NLOCAL, &NLOCAL, &wg.nc,
				&one_float,
				wg_wfc.c, &one_int, &one_int, ParaO.desc,
				wfc_k[ik].c, &one_int, &one_int, ParaO.desc,
				&zero_float,
				dm_k[ik].c, &one_int, &one_int, ParaO.desc);
		}
	}
	
	#ifdef TEST_DIAG
	{
		static int istep=0;
		ofstream ofs("dm_"+TO_STRING(istep)+"_"+TO_STRING(MY_RANK));
		if(GAMMA_ONLY_LOCAL)
		{
			ofs<<dm_gamma<<endl;
		}
		else
		{
			ofs<<dm_k<<endl;
		}
	}
	#endif

	return;
}
