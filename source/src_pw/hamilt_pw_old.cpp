#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"
#include "global.h"
#include "hamilt_pw.h"
#include "../module_base/blas_connector.h"
#include "../src_io/optical.h" // only get judgement to calculate optical matrix or not.
#include "myfunc.h"
#include "../module_base/timer.h"

int Hamilt_PW::moved = 0;

Hamilt_PW::Hamilt_PW()
{
    hpsi = new std::complex<double>[1];
    spsi = new std::complex<double>[1];
}

Hamilt_PW::~Hamilt_PW()
{
    delete[] hpsi;
    delete[] spsi;
}


void Hamilt_PW::allocate(
	const int &npwx,
	const int &npol,
	const int &nkb,
	const int &nrxx)
{
    ModuleBase::TITLE("Hamilt_PW","allocate");

	assert(npwx > 0);
	assert(npol > 0);
	assert(nkb >=0);

    delete[] hpsi;
    delete[] spsi;

    this->hpsi = new std::complex<double> [npwx * npol];
    this->spsi = new std::complex<double> [npwx * npol];

    ModuleBase::GlobalFunc::ZEROS(this->hpsi, npwx * npol);
    ModuleBase::GlobalFunc::ZEROS(this->spsi, npwx * npol);

    return;
}


void Hamilt_PW::init_k(const int ik)
{
    ModuleBase::TITLE("Hamilt_PW","init_k");

	// mohan add 2010-09-30
	// (1) Which spin to use.
	if(GlobalV::NSPIN==2)
	{
		GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
	}

	this->current_ik = ik;

	// (2) Take the local potential.
	for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	{
		GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);//mohan add 2007-11-12
	}

	// (3) Calculate nonlocal pseudopotential vkb
	//if (GlobalC::ppcell.nkb > 0 && !LINEAR_SCALING) xiaohui modify 2013-09-02
	if(GlobalC::ppcell.nkb > 0 && (GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")) //xiaohui add 2013-09-02. Attention...
	{
		GlobalC::ppcell.getvnl(ik, GlobalC::ppcell.vkb);
	}

	// (4) The number of wave functions.
	GlobalC::wf.npw = GlobalC::kv.ngk[ik];

	// (5) The index of plane waves.
	// (7) ik
	GlobalV::CURRENT_K = ik;

    return;
}


//----------------------------------------------------------------------
// Hamiltonian diagonalization in the subspace spanned
// by nstart states psi (atomic or random wavefunctions).
// Produces on output n_band eigenvectors (n_band <= nstart) in evc.
//----------------------------------------------------------------------
void Hamilt_PW::diagH_subspace(
    const int ik,
    const int nstart,
    const int n_band,
    const ModuleBase::ComplexMatrix &psi,
    ModuleBase::ComplexMatrix &evc,
    double *en)
{
    ModuleBase::TITLE("Hamilt_PW","diagH_subspace");
    ModuleBase::timer::tick("Hamilt_PW","diagH_subspace");

	assert(nstart!=0);
	assert(n_band!=0);

    ModuleBase::ComplexMatrix hc(nstart, nstart);
    ModuleBase::ComplexMatrix sc(nstart, nstart);
    ModuleBase::ComplexMatrix hvec(nstart,n_band);

	int dmin=0;
	int dmax=0;
	const int npw = GlobalC::kv.ngk[ik];

	if(GlobalV::NSPIN != 4)
	{
		dmin= npw;
		dmax = GlobalC::wf.npwx;
	}
	else
	{
		dmin = GlobalC::wf.npwx*GlobalV::NPOL;
		dmax = GlobalC::wf.npwx*GlobalV::NPOL;
	}

	std::complex<double> *aux=new std::complex<double> [dmax*nstart];
	std::complex<double> *paux = aux;
	std::complex<double> *ppsi = psi.c;

	this->h_psi(psi.c, aux, nstart);

	char trans1 = 'C';
	char trans2 = 'N';
	zgemm_(&trans1,&trans2,&nstart,&nstart,&dmin,&ModuleBase::ONE,psi.c,&dmax,aux,&dmax,&ModuleBase::ZERO,hc.c,&nstart);
	hc=transpose(hc,false);

	zgemm_(&trans1,&trans2,&nstart,&nstart,&dmin,&ModuleBase::ONE,psi.c,&dmax,psi.c,&dmax,&ModuleBase::ZERO,sc.c,&nstart);
	sc=transpose(sc,false);

	delete []aux;

	// Peize Lin add 2019-03-09
#ifdef __LCAO
#ifdef __MPI
	if(GlobalV::BASIS_TYPE=="lcao_in_pw")
	{
		auto add_Hexx = [&](const double alpha)
		{
			for (int m=0; m<nstart; ++m)
			{
				for (int n=0; n<nstart; ++n)
				{
					hc(m,n) += alpha * GlobalC::exx_lip.get_exx_matrix()[ik][m][n];
				}
			}
		};
		if(XC_Functional::get_func_type()==4)
		{
			if ( Exx_Global::Hybrid_Type::HF   == GlobalC::exx_lcao.info.hybrid_type ) // HF
			{
				add_Hexx(1);
			}
			else if (Exx_Global::Hybrid_Type::PBE0 == GlobalC::exx_lcao.info.hybrid_type || 
					Exx_Global::Hybrid_Type::HSE  == GlobalC::exx_lcao.info.hybrid_type) // PBE0 or HSE
			{
				add_Hexx(GlobalC::exx_global.info.hybrid_alpha);
			}
		}
	}
#endif
#endif

	if(GlobalV::NPROC_IN_POOL>1)
	{
		Parallel_Reduce::reduce_complex_double_pool( hc.c, nstart*nstart );
		Parallel_Reduce::reduce_complex_double_pool( sc.c, nstart*nstart );
	}

	// after generation of H and S matrix, diag them
    GlobalC::hm.diagH_LAPACK(nstart, n_band, hc, sc, nstart, en, hvec);


	// Peize Lin add 2019-03-09
#ifdef __LCAO
#ifdef __MPI
	if("lcao_in_pw"==GlobalV::BASIS_TYPE)
	{
		switch(GlobalC::exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				GlobalC::exx_lip.k_pack->hvec_array[ik] = hvec;
				break;
		}
	}
#endif
#endif

    //=======================
    //diagonize the H-matrix
    //=======================

// for tests
/*
		std::cout << std::setprecision(3);
		out.printV3(GlobalV::ofs_running,GlobalC::kv.kvec_c[ik]);
		out.printcm_norm("sc",sc,1.0e-4);
		out.printcm_norm("hvec",hvec,1.0e-4);
		out.printcm_norm("hc",hc,1.0e-4);
		std::cout << std::endl;
*/

	std::cout << std::setprecision(5);

//--------------------------
// KEEP THIS BLOCK FOR TESTS
//--------------------------
/*
	std::cout << "  hc matrix" << std::endl;
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			double a = hc(i,j).real();
			if(abs(a) < 1.0e-5) a = 0;
			std::cout << std::setw(6) << a;
		}
		std::cout << std::endl;
	}

	std::cout << "  sc matrix" << std::endl;
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{
			double a = sc(i,j).real();
			if(abs(a) < 1.0e-5) a = 0;
			std::cout << std::setw(6) << a;
		}
		std::cout << std::endl;
	}

	std::cout << "\n Band Energy" << std::endl;
	for(int i=0; i<GlobalV::NBANDS; i++)
	{
		std::cout << " e[" << i+1 << "]=" << en[i] * ModuleBase::Ry_to_eV << std::endl;
	}
*/
//--------------------------
// KEEP THIS BLOCK FOR TESTS
//--------------------------


	if((GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") && GlobalV::CALCULATION=="nscf" && !Optical::opt_epsilon2)
	{
		GlobalV::ofs_running << " Not do zgemm to get evc." << std::endl;
	}
	else if((GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
		&& ( GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax")) //pengfei 2014-10-13
	{
		// because psi and evc are different here,
		// I think if psi and evc are the same,
		// there may be problems, mohan 2011-01-01
		char transa = 'N';
		char transb = 'T';
		zgemm_( &transa,
				&transb,
				&dmax, // m: row of A,C
				&n_band, // n: col of B,C
				&nstart, // k: col of A, row of B
				&ModuleBase::ONE, // alpha
				psi.c, // A
				&dmax, // LDA: if(N) max(1,m) if(T) max(1,k)
				hvec.c, // B
				&n_band, // LDB: if(N) max(1,k) if(T) max(1,n)
				&ModuleBase::ZERO,  // belta
				evc.c, // C
				&dmax ); // LDC: if(N) max(1, m)
	}
	else
	{
		// As the evc and psi may refer to the same matrix, we first
		// create a temporary matrix to story the result. (by wangjp)
		// qianrui improve this part 2021-3-13
		char transa = 'N';
		char transb = 'T';
		ModuleBase::ComplexMatrix evctmp(n_band, dmin,false);
		zgemm_(&transa,&transb,&dmin,&n_band,&nstart,&ModuleBase::ONE,psi.c,&dmax,hvec.c,&n_band,&ModuleBase::ZERO,evctmp.c,&dmin);
		for(int ib=0; ib<n_band; ib++)
		{
			for(int ig=0; ig<dmin; ig++)
			{
				evc(ib,ig) = evctmp(ib,ig);
			}
		}
	}
    //out.printr1_d("en",en,n_band);

//	std::cout << "\n bands" << std::endl;
//	for(int ib=0; ib<n_band; ib++)
//	{
//		std::cout << " ib=" << ib << " " << en[ib] * ModuleBase::Ry_to_eV << std::endl;
//	}

    //out.printcm_norm("hvec",hvec,1.0e-8);

    ModuleBase::timer::tick("Hamilt_PW","diagH_subspace");
    return;
}


void Hamilt_PW::h_1psi( const int npw_in, const std::complex < double> *psi,
                        std::complex<double> *hpsi, std::complex < double> *spsi)
{
    this->h_psi(psi, hpsi);

    for (int i=0;i<npw_in;i++)
    {
        spsi[i] = psi[i];
    }
    return;
}


void Hamilt_PW::s_1psi
(
    const int dim,
    const std::complex<double> *psi,
    std::complex<double> *spsi
)
{
    for (int i=0; i<dim; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}


void Hamilt_PW::h_psi(const std::complex<double> *psi_in, std::complex<double> *hpsi, const int m)
{
    ModuleBase::timer::tick("Hamilt_PW","h_psi");
    int i = 0;
    int j = 0;
    int ig= 0;
	const int ik = this->current_ik;
	const double tpiba2 = GlobalC::ucell.tpiba2;
	const int npwx = GlobalC::wf.npwx;
	const int npw = GlobalC::wf.npw;
	const int nrxx = GlobalC::rhopw->nrxx;

	//if(GlobalV::NSPIN!=4) ModuleBase::GlobalFunc::ZEROS(hpsi, npw);
	//else ModuleBase::GlobalFunc::ZEROS(hpsi, GlobalC::wf.npwx * GlobalV::NPOL);//added by zhengdy-soc
	const int dmax = npwx * GlobalV::NPOL;

	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	std::complex<double> *tmhpsi;
	const std::complex<double> *tmpsi_in;
 	if(GlobalV::T_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		for(int ib = 0 ; ib < m; ++ib)
		{
			for(ig = 0;ig < npw; ++ig)
			{
				tmhpsi[ig] = GlobalC::wfcpw->getgk2(ik,ig) * tpiba2 * tmpsi_in[ig];
			}
			if(GlobalV::NSPIN==4){
				for(ig=npw; ig < npwx; ++ig)
				{
					tmhpsi[ig] = 0;
				}
				tmhpsi +=npwx;
				tmpsi_in += npwx;
				for (ig = 0;ig < npw ;++ig)
				{
					tmhpsi[ig] = GlobalC::wfcpw->getgk2(ik,ig) * tpiba2 * tmpsi_in[ig];
				}
				for(ig=npw; ig < npwx; ++ig)
				{
					tmhpsi[ig] =0;
				}
			}
			tmhpsi += npwx;
			tmpsi_in += npwx;
		}
	}

	//------------------------------------
	//(2) the local potential.
	//-----------------------------------
	ModuleBase::timer::tick("Hamilt_PW","vloc");
	std::complex<double>* porter = new std::complex<double>[GlobalC::wfcpw->nmaxgr];
	if(GlobalV::VL_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		for(int ib = 0 ; ib < m; ++ib)
		{
			if(GlobalV::NSPIN!=4){
				GlobalC::wfcpw->recip2real(tmpsi_in, porter, ik);
				for (int ir=0; ir< nrxx; ir++)
				{
					porter[ir] *=  GlobalC::pot.vr_eff1[ir];
				}
				GlobalC::wfcpw->real2recip(porter, tmhpsi, ik, true);
			}
			else
			{
				std::complex<double>* porter1 = new std::complex<double>[nrxx];
				// fft to real space and doing things.
				GlobalC::wfcpw->recip2real(tmpsi_in, porter, ik);
				GlobalC::wfcpw->recip2real(tmpsi_in+npwx, porter1, ik);
				std::complex<double> sup,sdown;
				for (int ir=0; ir< nrxx; ir++)
				{
					sup = porter[ir] * (GlobalC::pot.vr_eff(0,ir) + GlobalC::pot.vr_eff(3,ir)) +
						porter1[ir] * (GlobalC::pot.vr_eff(1,ir) - std::complex<double>(0.0,1.0) * GlobalC::pot.vr_eff(2,ir));
					sdown = porter1[ir] * (GlobalC::pot.vr_eff(0,ir) - GlobalC::pot.vr_eff(3,ir)) +
					porter[ir] * (GlobalC::pot.vr_eff(1,ir) + std::complex<double>(0.0,1.0) * GlobalC::pot.vr_eff(2,ir));
					porter[ir] = sup;
					porter1[ir] = sdown;
				}
				// (3) fft back to G space.
				GlobalC::wfcpw->real2recip(porter, tmhpsi, ik, true);
				GlobalC::wfcpw->real2recip(porter1, tmhpsi+npwx, ik, true);

				delete[] porter1;
			}
			tmhpsi += dmax;
			tmpsi_in += dmax;
		}
	}
	ModuleBase::timer::tick("Hamilt_PW","vloc");

	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	ModuleBase::timer::tick("Hamilt_PW","vnl");
	if(GlobalV::VNL_IN_H)
	{
		if ( GlobalC::ppcell.nkb > 0)
		{
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			//qianrui optimize 2021-3-31
			int nkb=GlobalC::ppcell.nkb;
			ModuleBase::ComplexMatrix becp(GlobalV::NPOL * m, nkb, false);
			char transa = 'C';
			char transb = 'N';
			if(m==1 && GlobalV::NPOL==1)
			{
				int inc = 1;
				zgemv_(&transa, &npw, &nkb, &ModuleBase::ONE, GlobalC::ppcell.vkb.c, &GlobalC::wf.npwx, psi_in, &inc, &ModuleBase::ZERO, becp.c, &inc);
			}
			else
			{
				int npm = GlobalV::NPOL * m;
				zgemm_(&transa,&transb,&nkb,&npm,&npw,&ModuleBase::ONE,GlobalC::ppcell.vkb.c,&GlobalC::wf.npwx,psi_in,&GlobalC::wf.npwx,&ModuleBase::ZERO,becp.c,&nkb);
				//add_nonlocal_pp is moddified, thus tranpose not needed here.
				//if(GlobalV::NONCOLIN)
				//{
				//	ModuleBase::ComplexMatrix partbecp(GlobalV::NPOL, nkb ,false);
				//	for(int ib = 0; ib < m; ++ib)
				//	{
//
				//		for ( i = 0;i < GlobalV::NPOL;i++)
				//			for (j = 0;j < nkb;j++)
				//				partbecp(i, j) = tmbecp[i*nkb+j];
				//		for (j = 0; j < nkb; j++)
				//			for (i = 0;i < GlobalV::NPOL;i++)
				//				tmbecp[j*GlobalV::NPOL+i] = partbecp(i, j);
				//		tmbecp += GlobalV::NPOL * nkb;
				//	}
				//}
			}

			Parallel_Reduce::reduce_complex_double_pool( becp.c, nkb * GlobalV::NPOL * m);

			this->add_nonlocal_pp(hpsi, becp.c, m);
			//======================================================================
			/*std::complex<double> *becp = new std::complex<double>[ GlobalC::ppcell.nkb * GlobalV::NPOL ];
			ModuleBase::GlobalFunc::ZEROS(becp,GlobalC::ppcell.nkb * GlobalV::NPOL);
			for (i=0;i< GlobalC::ppcell.nkb;i++)
			{
				const std::complex<double>* p = &GlobalC::ppcell.vkb(i,0);
				const std::complex<double>* const p_end = p + npw;
				const std::complex<double>* psip = psi_in;
				for (;p<p_end;++p,++psip)
				{
					if(!GlobalV::NONCOLIN) becp[i] += psip[0]* conj( p[0] );
					else{
						becp[i*2] += psip[0]* conj( p[0] );
						becp[i*2+1] += psip[GlobalC::wf.npwx]* conj( p[0] );
					}
				}
			}
			Parallel_Reduce::reduce_complex_double_pool( becp, GlobalC::ppcell.nkb * GlobalV::NPOL);
			this->add_nonlocal_pp(hpsi, becp);
			delete[] becp;*/
			//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		}
	}
	ModuleBase::timer::tick("Hamilt_PW","vnl");
	//------------------------------------
	// (4) the metaGGA part
	//------------------------------------
	ModuleBase::timer::tick("Hamilt_PW","meta");
	if(XC_Functional::get_func_type() == 3)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		for(int ib = 0; ib < m; ++ib)
		{
			for(int j=0; j<3; j++)
			{
				for (int ig = 0;ig < GlobalC::kv.ngk[GlobalV::CURRENT_K] ; ig++)
				{
					double fact = GlobalC::wfcpw->getgpluskcar(ik,ig)[j] * GlobalC::ucell.tpiba;
					porter[ig] = tmpsi_in[ig] * complex<double>(0.0,fact);
				}
				
				GlobalC::wfcpw->recip2real(porter,porter,ik);

				for (int ir = 0; ir < nrxx; ir++)
				{
					porter[ir] *= GlobalC::pot.vofk(GlobalV::CURRENT_SPIN,ir);
				}
				GlobalC::wfcpw->real2recip(porter,porter,ik);

				for (int ig = 0;ig < npw ; ig++)
				{
					double fact = GlobalC::wfcpw->getgpluskcar(ik,ig)[j] * GlobalC::ucell.tpiba;
					tmhpsi[ig] -= complex<double>(0.0,fact) * porter[ig];
				}
			}//x,y,z directions
		}
	}
	delete[] porter;
	ModuleBase::timer::tick("Hamilt_PW","meta");
    ModuleBase::timer::tick("Hamilt_PW","h_psi");
    return;
}

//--------------------------------------------------------------------------
// this function sum up each non-local pseudopotential located on each atom,
//--------------------------------------------------------------------------
void Hamilt_PW::add_nonlocal_pp(
	std::complex<double> *hpsi_in,
	const std::complex<double> *becp,
	const int m)
{
    ModuleBase::timer::tick("Hamilt_PW","add_nonlocal_pp");

	// number of projectors
	int nkb = GlobalC::ppcell.nkb;

	std::complex<double> *ps  = new std::complex<double> [nkb * GlobalV::NPOL * m];
    ModuleBase::GlobalFunc::ZEROS(ps, GlobalV::NPOL * m * nkb);

    int sum = 0;
    int iat = 0;
    if(GlobalV::NSPIN!=4)
	{
		for (int it=0; it<GlobalC::ucell.ntype; it++)
		{
			const int nproj = GlobalC::ucell.atoms[it].nh;
			for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				// each atom has nproj, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip=0; ip<nproj; ip++)
				{
					for (int ip2=0; ip2<nproj; ip2++)
					{
						for(int ib = 0; ib < m ; ++ib)
						{
							ps[(sum + ip2) * m + ib] +=
							GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip2)
							* becp[ib * nkb + sum + ip];
						}//end ib
					}// end ih
				}//end jh
				sum += nproj;
				++iat;
			} //end na
		} //end nt
	}
	else
	{
		for (int it=0; it<GlobalC::ucell.ntype; it++)
		{
			int psind=0;
			int becpind=0;
			std::complex<double> becp1=std::complex<double>(0.0,0.0);
			std::complex<double> becp2=std::complex<double>(0.0,0.0);

			const int nproj = GlobalC::ucell.atoms[it].nh;
			for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				// each atom has nproj, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip=0; ip<nproj; ip++)
				{
					for (int ip2=0; ip2<nproj; ip2++)
					{
						for(int ib = 0; ib < m ; ++ib)
						{
							psind = (sum+ip2) * 2 * m + ib * 2;
							becpind = ib*nkb*2 + sum + ip;
							becp1 =  becp[becpind];
							becp2 =  becp[becpind + nkb];
							ps[psind] += GlobalC::ppcell.deeq_nc(0, iat, ip2, ip) * becp1
								+GlobalC::ppcell.deeq_nc(1, iat, ip2, ip) * becp2;
							ps[psind +1] += GlobalC::ppcell.deeq_nc(2, iat, ip2, ip) * becp1
								+GlobalC::ppcell.deeq_nc(3, iat, ip2, ip) * becp2;
						}//end ib
					}// end ih
				}//end jh
				sum += nproj;
				++iat;
			} //end na
		} //end nt
	}

	/*
    for (int ig=0;ig<GlobalC::wf.npw;ig++)
    {
        for (int i=0;i< GlobalC::ppcell.nkb;i++)
        {
            hpsi_in[ig]+=ps[i]*GlobalC::ppcell.vkb(i,ig);
        }
    }
	*/


	// use simple method.
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//qianrui optimize 2021-3-31
	char transa = 'N';
	char transb = 'T';
	if(GlobalV::NPOL==1 && m==1)
	{
		int inc = 1;
		zgemv_(&transa,
			&GlobalC::wf.npw,
			&GlobalC::ppcell.nkb,
			&ModuleBase::ONE,
			GlobalC::ppcell.vkb.c,
			&GlobalC::wf.npwx,
			ps,
			&inc,
			&ModuleBase::ONE,
			hpsi_in,
			&inc);
	}
	else
	{
		int npm = GlobalV::NPOL*m;
		zgemm_(&transa,
			&transb,
			&GlobalC::wf.npw,
			&npm,
			&GlobalC::ppcell.nkb,
			&ModuleBase::ONE,
			GlobalC::ppcell.vkb.c,
			&GlobalC::wf.npwx,
			ps,
			&npm,
			&ModuleBase::ONE,
			hpsi_in,
			&GlobalC::wf.npwx);
	}

	//======================================================================
	/*if(!GlobalV::NONCOLIN)
	for(int i=0; i<GlobalC::ppcell.nkb; i++)
	{
		std::complex<double>* p = &GlobalC::ppcell.vkb(i,0);
		std::complex<double>* p_end = p + GlobalC::wf.npw;
		std::complex<double>* hp = hpsi_in;
		std::complex<double>* psp = &ps[i];
		for (;p<p_end;++p,++hp)
		{
			hp[0] += psp[0] * p[0];
		}
	}
	else
	for(int i=0; i<GlobalC::ppcell.nkb; i++)
	{
		std::complex<double>* p = &GlobalC::ppcell.vkb(i,0);
		std::complex<double>* p_end = p + GlobalC::wf.npw;
		std::complex<double>* hp = hpsi_in;
		std::complex<double>* hp1 = hpsi_in + GlobalC::wf.npwx;
		std::complex<double>* psp = &ps[i*2];
		for (;p<p_end;p++,++hp,++hp1)
		{
			hp[0] += psp[0] * (p[0]);
			hp1[0] += psp[1] * (p[0]);
		}
	}*/
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	delete[] ps;
    ModuleBase::timer::tick("Hamilt_PW","add_nonlocal_pp");
    return;
}

void Hamilt_PW::diag_zheev(const int &npw_in, ModuleBase::ComplexMatrix &psi, const int &nband, double *em, double *err)
{
    ModuleBase::TITLE("Hamilt_PW","diag_zheev");
    assert(nband < npw_in) ;

    // if flag =0, this means success.
    // RedM means Reduced matrix, because the dimension of plane wave
    // is much larger than nbands.

    ModuleBase::ComplexMatrix RedM(nband, nband);
    std::complex<double> * eta =  new std::complex<double>[npw_in] ;
    std::complex<double> * hpsi1 =  new std::complex<double>[npw_in] ;
    std::complex<double> * spsi1 =  new std::complex<double>[npw_in] ;

	ModuleBase::GlobalFunc::ZEROS(eta, npw_in);
	ModuleBase::GlobalFunc::ZEROS(hpsi1, npw_in);
	ModuleBase::GlobalFunc::ZEROS(spsi1, npw_in);

    double * tmpen = new double[nband] ;

    assert(eta  !=  0) ;
    assert(hpsi1  !=    0) ;
    assert(spsi1    !=  0) ;
    assert(tmpen != 0) ;

    // <j|H|i>, where j < = i is calculated.
    // first calculate eta =|i>, and hpsi = H|i>

    std::complex<double> tmp ;
    std::complex<double> tmp1 ;

    // calculate tmpen[i]

    for (int i = 0; i < nband ; i++)
    {
        dcopy(psi, i, eta);
        h_1psi(npw_in, eta, hpsi1, spsi1) ;

        tmp = ModuleBase::ZERO ;
        tmp1 = ModuleBase::ZERO ;

        for (int ig = 0; ig < npw_in; ig++)
        {
            tmp += conj(eta[ig]) * hpsi1[ig] ;
            tmp1 += conj(eta[ig]) * eta[ig] ;
        }

		tmp = tmp * conj(tmp1) / (norm(tmp1));

		tmpen[i] = tmp.real() ;

        for (int j = 0; j <= i; j++)
        {
            // calculate H[i,j] = <i|eta> where |eta> = H|j>
            // RedM(j, i) = <j|H|i>

            tmp = ModuleBase::ZERO ;

            for (int ig = 0; ig < npw_in; ig++)
            {
                tmp += conj(psi(j, ig)) * hpsi1[ig];
            }

            RedM(j, i) = tmp ;

            if (i != j) RedM(i, j) = conj(tmp) ;
        }
    }

    // This is for calling zheev.
    char  jobz = 'V' ;  // eigenvalues and eigenvectors

    char  uplo = 'U' ;  // the upper is stored.

    const int lwork = 2 * nband;

    std::complex<double> * work = new std::complex<double>[lwork]() ;

    double * rwork = new double[3*nband-2] ;

    int info = 0 ;

    // diag by calling zheev
    // basis are psi_0, psi_1, psi_2...

    LapackConnector::zheev(jobz, uplo, nband, RedM, nband, em, work, lwork, rwork, &info);

    delete[] eta  ;

    delete[] hpsi1 ;

    delete[] spsi1 ;

    //  change back to plane wave basis.

//  std::cout << " calling zheev is performed " << std::endl ;
//  std::cout << " zheev output info... " << std::endl;

    // delete the allocated data array.

    delete[] work ;

    delete[] rwork ;

    // std::cout the infomation from zheev...

    // std::cout the infomation from zheev...
    if (info == 0)
    {
        std::cout << "  successful exit of zheev " << std::endl ;
    }
    else if (info < 0)
    {
        std::cout << " the i-th argument had an illegal value. info =  " << info << std::endl ;
    }
    else
    {
        std::cout << "the algorithm failed to converge. info = " << info << std::endl ;
    }

    ModuleBase::ComplexMatrix kpsi(nband, npw_in);

    // kpsi = c_1 \psi_1 + c_2 \psi_2 + \cdots + c_N \psi_N
    // For the mth wavefunction, |psi(m)> = U(1, m) \psi_1 + ... U(k, m) \psi_k + ...
    // So, |psi(m)>_j means the jth element of |psi(m)>
    //  We have |psi(m)>_j =  U(1, m) \psi(1, j) + U(2, m) \psi(2, j) + ...
    // = \sum_{k=1}^{nband} U(k, m) \psi(k, j)

    // Store the wavefunction in kpsi

    for (int m = 0; m < nband; m++)
    {
        for (int j = 0; j < npw_in; j++)
        {
            tmp = ModuleBase::ZERO ;

            for (int  k = 0; k < nband; k++)
            {
                tmp += RedM(k, m) * psi(k, j);
            }

            kpsi(m, j) = tmp ;
        }
    }


    // update the wavefunction of psi
    for (int m = 0; m < nband ; m++)
    {
        for (int ig = 0; ig < npw_in; ig++)
        {
            psi(m, ig) = kpsi(m, ig);
        }
    }

    // calculate error of the last results.
    // err = || H\psi -E\psi||

    std::cout << " callilng cal_err " << std::endl ;

    cal_err(npw_in, psi, nband, em, err);


    // out put the results.
    std::cout<<std::setw(6)<<"Bands"
        <<std::setw(12)<<"energy(ev)"
        <<std::setw(12)<<"err"
        <<std::setw(25)<<"||H * psi - E * psi ||\n";

    for (int m = 0; m < 5; m++)
    {
        std::cout << std::setw(6) << m
             << std::setw(12) << em[m] * ModuleBase::Ry_to_eV
             << std::setw(12) << err[m]
             << std::setw(25) << tmpen[m] * ModuleBase::Ry_to_eV << std::endl ;
    }

    std::cout << " end of diag_zheev " << std::endl ;

    return;
}

void Hamilt_PW::cal_err
(
    const int &npw_in,
    ModuleBase::ComplexMatrix &psi,
    const int &nband,
    double *em,
    double *err
)
{
//	ModuleBase::TITLE("Hamilt_PW", "cal_err");
//	std::cout << "\n npw_in = " << npw_in << std::endl;

    assert(nband < npw_in);
    ModuleBase::timer::tick("Hamilt_PW", "cal_err") ;

    std::complex<double> *psitmp =  new std::complex<double>[npw_in]();
    std::complex<double> *hpsitmp =  new std::complex<double>[npw_in]();
    std::complex<double> *spsitmp =  new std::complex<double>[npw_in]();

    std::complex<double> tmp1 ;

    for (int m = 0; m < nband; m++)
    {
//		std::cout << "\n m = " << m << std::endl;
        dcopy(psi, m, psitmp) ;
        h_1psi(npw_in, psitmp, hpsitmp, spsitmp);

        std::complex<double> tmp = ModuleBase::ZERO;

        for (int ig=0;  ig<npw_in; ig++)
        {
            tmp1 =  hpsitmp[ig] - em[m] * psitmp[ig];
            tmp +=  conj(tmp1) * tmp1;
        }

        // err[m] = ||H\psitmp - \lambda_m \psitmp||
        err[m] = sqrt( tmp.real() ) ;
    }

    delete[] psitmp;
    delete[] hpsitmp;
    delete[] spsitmp;

//	std::cout << " calculate error of the wavefunctions " << std::endl ;
    ModuleBase::timer::tick("Hamilt_PW", "cal_err") ;
    return;
}