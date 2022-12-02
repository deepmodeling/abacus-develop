#include "wavefunc.h"
#include "global.h"
#include "../src_lcao/wavefunc_in_pw.h"
#include "../src_io/winput.h"
#include "../src_io/chi0_hilbert.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hamilt/hamilt_pw.h"

wavefunc::wavefunc()
{
	out_wfc_pw = 0;
}

wavefunc::~wavefunc()
{
	if(GlobalV::test_deconstructor)
	{
		std::cout << " ~wavefunc()" << std::endl;
	}
	if(this->irindex == nullptr) 
	{
		delete[] this->irindex;		
		this->irindex=nullptr;
	}
}

psi::Psi<std::complex<double>>* wavefunc::allocate(const int nks)
{	
	ModuleBase::TITLE("wavefunc","allocate");

	this->npwx = GlobalC::wfcpw->npwk_max;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"npwx",npwx);

	assert(npwx > 0);
	assert(nks > 0);

	// allocate for kinetic energy

	// if use spin orbital, do not double nks but double allocate evc and wanf2.
	int prefactor = 1;
	if(GlobalV::NSPIN==4) prefactor = GlobalV::NPOL;//added by zhengdy-soc

	const int nks2 = nks;

	psi::Psi<std::complex<double>>* psi_out = nullptr;
	if(GlobalV::CALCULATION=="nscf" && GlobalC::wf.mem_saver==1)
	{
		//initial psi rather than evc
		psi_out = new psi::Psi<std::complex<double>>(1, GlobalV::NBANDS, npwx * GlobalV::NPOL, GlobalC::kv.ngk.data());
		if(GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			wanf2[0].create(GlobalV::NLOCAL, npwx * GlobalV::NPOL);
			std::cout << " Memory for wanf2 (MB): " << 
				ModuleBase::Memory::record("wavefunc","wanf2",GlobalV::NLOCAL*(prefactor*npwx),"complexmatrix") << std::endl;
		}
		std::cout << " MEMORY FOR PSI (MB)  : " << 
			ModuleBase::Memory::record("wavefunc","psi",GlobalV::NBANDS*(prefactor*npwx),"complexmatrix") << std::endl;
	}
	else if(GlobalV::BASIS_TYPE!="pw")
	{
		this->evc = new ModuleBase::ComplexMatrix [nks2];
		this->wanf2 = new ModuleBase::ComplexMatrix [nks2];
		for (int ik = 0; ik < nks2; ik++)
		{
			this->evc[ik].create(GlobalV::NBANDS, npwx * GlobalV::NPOL);//added by zhengdy-soc
			//Mohan add 2010-1-10
			if((GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") || winput::out_spillage==2)
			{
				this->wanf2[ik].create(GlobalV::NLOCAL, npwx * GlobalV::NPOL);//added by zhengdy-soc
			}
		};
		std::cout << " MEMORY FOR PSI (MB)  : " << 
		ModuleBase::Memory::record("wavefunc","evc",nks2*GlobalV::NBANDS*(prefactor*npwx),"complexmatrix") << std::endl;
	}
	else
	{
		//initial psi rather than evc
		psi_out = new psi::Psi<std::complex<double>>(nks2, GlobalV::NBANDS, npwx * GlobalV::NPOL, GlobalC::kv.ngk.data());
		std::cout << " MEMORY FOR PSI (MB)  : " << 
		ModuleBase::Memory::record("wavefunc","psi",nks2*GlobalV::NBANDS*(prefactor*npwx),"complexmatrix") << std::endl;
	}
	return psi_out;

	//showMemStats();
}

//===================================================================
// This routine computes an estimate of the start_ wavefunctions
// from superposition of atomic wavefunctions or random wave functions.
//===================================================================
void wavefunc::wfcinit(psi::Psi<std::complex<double>>* psi_in)
{
    ModuleBase::TITLE("wavefunc","wfcinit");
    ModuleBase::timer::tick("wavefunc","wfcinit");

    this->wfcinit_k(psi_in);

    GlobalC::en.demet = 0.0;

    ModuleBase::timer::tick("wavefunc","wfcinit");
    return;
}

int wavefunc::get_starting_nw(void)const
{
    if (init_wfc == "file")
    {
		throw std::runtime_error("wavefunc::get_starting_nw. start_ wfc from file: not implemented yet! "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__)); 	// Peize Lin change 2019-05-01
        //ModuleBase::WARNING_QUIT("wfcinit_k","\n start_ wfc from file: not implemented yet!");
        //**********************************************************************
        // ... read the wavefunction into memory (if it is not done in c_bands)
        //**********************************************************************
    }
    else if (init_wfc.substr(0,6) == "atomic")
    {
        if (GlobalC::ucell.natomwfc >= GlobalV::NBANDS)
        {
            if(GlobalV::test_wf)GlobalV::ofs_running << " Start wave functions are all pseudo atomic wave functions." << std::endl;
        }
        else
        {
            if(GlobalV::test_wf)GlobalV::ofs_running << " Start wave functions are atomic + "
            << GlobalV::NBANDS - GlobalC::ucell.natomwfc
            << " random wave functions." << std::endl;
        }
        return max(GlobalC::ucell.natomwfc,  GlobalV::NBANDS);
    }
    else if (init_wfc == "random")
    {
        if(GlobalV::test_wf)GlobalV::ofs_running << " Start wave functions are all random." << std::endl;
        return GlobalV::NBANDS;
    }
    else
    {
		throw std::runtime_error("wavefunc::get_starting_nw. Don't know what to do! Please Check source code! "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__)); 	// Peize Lin change 2019-05-01
        //ModuleBase::WARNING_QUIT("get_starting_nw","Don't know what to do! Please Check source code!");
    }
}



#ifdef __LCAO
//We are not goint to support lcao_in_paw until
//the obsolete GlobalC::hm is replaced by the 
//refactored moeules (psi, hamilt, etc.)
/*
void wavefunc::LCAO_in_pw_k(const int &ik, ModuleBase::ComplexMatrix &wvf)
{
	ModuleBase::TITLE("wavefunc","LCAO_in_pw_k");
	ModuleBase::timer::tick("wavefunc","LCAO_in_pw_k");

	assert(GlobalV::BASIS_TYPE=="lcao_in_pw");

	static bool ltable = false;

	if(!ltable)
	{
		this->table_local.create(GlobalC::ucell.ntype, GlobalC::ucell.nmax_total, GlobalV::NQX);

		// GlobalC::ORB.orbital_file: file name of the numerical atomic orbitals (NAOs)
		// table_local: generate one-dimensional table for NAOs
		Wavefunc_in_pw::make_table_q(GlobalC::ORB.orbital_file, this->table_local);
		ltable = true;
	}

	Wavefunc_in_pw::produce_local_basis_in_pw(ik, wvf, this->table_local);

	//-------------------------------------------------------------
	// (2) diago to get ElecState::ekb, then the weights can be calculated.
	//-------------------------------------------------------------
    GlobalC::hm.hpw.allocate(this->npwx, GlobalV::NPOL, GlobalC::ppcell.nkb, GlobalC::wfcpw->nrxx);
	GlobalC::hm.hpw.init_k(ik);

	//GlobalC::hm.diagH_subspace(ik ,GlobalV::NLOCAL, GlobalV::NBANDS, wvf, wvf, ekb[ik]);
//	for(int ib=0; ib<GlobalV::NBANDS; ib++)
//	{
//		std::cout << " ib=" << ib << " e=" << ekb[ik][ib] << std::endl;
//	}

//	ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"CONSTRUCT_LOCAL_BASIS_IN_PW");

	ModuleBase::timer::tick("wavefunc","LCAO_in_pw_k");
	return;
}


void wavefunc::LCAO_in_pw_k_q(const int &ik, ModuleBase::ComplexMatrix &wvf, ModuleBase::Vector3<double> q)   // pengfei  2016-11-23
{
	ModuleBase::TITLE("wavefunc","LCAO_in_pw_k_q");
	ModuleBase::timer::tick("wavefunc","LCAO_in_pw_k_q");
	//assert(LOCAL_BASIS==4); xiaohui modify 2013-09-01
	assert(GlobalV::BASIS_TYPE=="lcao_in_pw"); //xiaohui add 2013-09-01. Attention! How about "BASIS_TYPE=="lcao""???

	Wavefunc_in_pw::produce_local_basis_q_in_pw(ik, wvf, this->table_local, q);

	ModuleBase::timer::tick("wavefunc","LCAO_in_pw_k_q");
	return;
}
*/
#endif

void wavefunc::diago_PAO_in_pw_k2(const int &ik, psi::Psi<std::complex<double>> &wvf, hamilt::Hamilt<double>* phm_in)
{
	ModuleBase::TITLE("wavefunc","diago_PAO_in_pw_k2");
	// (6) Prepare for atmoic orbitals or random orbitals
	const int starting_nw = this->get_starting_nw();
	if(starting_nw == 0) return;
	assert(starting_nw > 0);

	const int nbasis = wvf.get_nbasis();
	const int nbands = wvf.get_nbands();
	const int current_nbasis = GlobalC::kv.ngk[ik];

	//special case here! use Psi(k-1) for the initialization of Psi(k)
	//this method should be tested.
	/*if(GlobalV::CALCULATION == "nscf" && GlobalC::ucell.natomwfc == 0 && ik>0)
	{
		//this is memsaver case
		if(wvf.get_nk() == 1)
		{
			return;
		}
		else
		{
			ModuleBase::GlobalFunc::COPYARRAY(&wvf(ik-1, 0, 0), &wvf(ik, 0, 0), wvf.get_nbasis()* wvf.get_nbands());
			return;
		}
	}
	*/

	ModuleBase::ComplexMatrix wfcatom(starting_nw, nbasis);//added by zhengdy-soc
	if(GlobalV::test_wf)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "starting_nw", starting_nw);
	if(init_wfc.substr(0,6)=="atomic")
	{
		this->atomic_wfc(ik, current_nbasis, GlobalC::ucell.lmax_ppwf, wfcatom, GlobalC::ppcell.tab_at, GlobalV::NQX, GlobalV::DQ);
		if( init_wfc == "atomic+random" && starting_nw == GlobalC::ucell.natomwfc )//added by qianrui 2021-5-16
		{
			this->atomicrandom(wfcatom,0,starting_nw,ik, GlobalC::wfcpw);
		}

		//====================================================
		// If not enough atomic wfc are available, complete
		// with random wfcs
		//====================================================
		this->random(wfcatom, GlobalC::ucell.natomwfc, nbands, ik, GlobalC::wfcpw);
	}
	else if(init_wfc=="random")
	{
		this->random(wfcatom,0,nbands,ik, GlobalC::wfcpw);
	}

	// (7) Diago with cg method.
	std::vector<double> etatom(starting_nw, 0.0);
	//if(GlobalV::DIAGO_TYPE == "cg") xiaohui modify 2013-09-02
	if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
	{
		if(phm_in!= nullptr)
		{
			// hsolver::DiagoIterAssist<double>::diagH_subspace_init(phm_in,
      //                          wfcatom,
      //                          wvf,
      //                          etatom.data());
			hsolver::DiagoIterAssist<double>::diagH_subspace_init(phm_in,
                         wfcatom.c,
				   							 wfcatom.nr,
				   							 wfcatom.nc,
                         wvf,
                         etatom.data());
			return;
		}
		else
		{
			ModuleBase::WARNING_QUIT("wavefunc","Psi does not exist!");
			//this diagonalization method is obsoleted now
			//GlobalC::hm.diagH_subspace(ik ,starting_nw, nbands, wfcatom, wfcatom, etatom.data());
		}
	}

	assert(nbands <= wfcatom.nr);
	for (int ib=0; ib<nbands; ib++)
	{
		for (int ig=0; ig<nbasis; ig++)
		{
			wvf(ib, ig) = wfcatom(ib, ig);
		}
	}
}

void wavefunc::diago_PAO_in_pw_k2_device(const psi::DEVICE_CPU* ctx, const int &ik, psi::Psi<std::complex<double>, psi::DEVICE_CPU> &wvf, hamilt::Hamilt<double, psi::DEVICE_CPU>* phm_in)
{
	this->diago_PAO_in_pw_k2(ik, wvf, phm_in);
}

#if ((defined __CUDA) || (defined __ROCM))
void wavefunc::diago_PAO_in_pw_k2_device(const psi::DEVICE_GPU* ctx, const int &ik, psi::Psi<std::complex<double>, psi::DEVICE_GPU> &wvf, hamilt::Hamilt<double, psi::DEVICE_GPU>* phm_in)
{
	ModuleBase::TITLE("wavefunc","diago_PAO_in_pw_k2");
	// (6) Prepare for atmoic orbitals or random orbitals
	const int starting_nw = this->get_starting_nw();
	if(starting_nw == 0) return;
	assert(starting_nw > 0);

	const int nbasis = wvf.get_nbasis();
	const int nbands = wvf.get_nbands();
	const int current_nbasis = GlobalC::kv.ngk[ik];

	ModuleBase::ComplexMatrix wfcatom(starting_nw, nbasis);//added by zhengdy-soc
	if(GlobalV::test_wf)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "starting_nw", starting_nw);
	if(init_wfc.substr(0,6)=="atomic")
	{
		this->atomic_wfc(ik, current_nbasis, GlobalC::ucell.lmax_ppwf, wfcatom, GlobalC::ppcell.tab_at, GlobalV::NQX, GlobalV::DQ);
		if( init_wfc == "atomic+random" && starting_nw == GlobalC::ucell.natomwfc )//added by qianrui 2021-5-16
		{
			this->atomicrandom(wfcatom,0,starting_nw,ik, GlobalC::wfcpw);
		}

		//====================================================
		// If not enough atomic wfc are available, complete
		// with random wfcs
		//====================================================
		this->random(wfcatom, GlobalC::ucell.natomwfc, nbands, ik, GlobalC::wfcpw);
	}
	else if(init_wfc=="random")
	{
		this->random(wfcatom,0,nbands,ik, GlobalC::wfcpw);
	}

	// (7) Diago with cg method.
	std::vector<double> etatom(starting_nw, 0.0);
	//if(GlobalV::DIAGO_TYPE == "cg") xiaohui modify 2013-09-02
	if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
	{
		if(phm_in!= nullptr)
		{
			// hsolver::DiagoIterAssist<double>::diagH_subspace_init(phm_in,
      //                          wfcatom,
      //                          wvf,
      //                          etatom.data());
			std::complex<double> *d_wfcatom = nullptr;
			psi::DEVICE_CPU * cpu_ctx = {};
			psi::DEVICE_GPU * gpu_ctx = {};
			resmem_complex_op()(gpu_ctx, d_wfcatom, wfcatom.nr * wfcatom.nc);
			syncmem_complex_h2d_op()(gpu_ctx, cpu_ctx, d_wfcatom, wfcatom.c, wfcatom.nr * wfcatom.nc);
			hsolver::DiagoIterAssist<double, psi::DEVICE_GPU>::diagH_subspace_init(
												 phm_in,
                         d_wfcatom,
				   							 wfcatom.nr,
				   							 wfcatom.nc,
                         wvf,
                         etatom.data());
			delmem_complex_op()(gpu_ctx, d_wfcatom);
			return;
		}
		else
		{
			//this diagonalization method is obsoleted now
			//GlobalC::hm.diagH_subspace(ik ,starting_nw, nbands, wfcatom, wfcatom, etatom.data());
		}
	}

	assert(nbands <= wfcatom.nr);
	for (int ib=0; ib<nbands; ib++)
	{
		for (int ig=0; ig<nbasis; ig++)
		{
			wvf(ib, ig) = wfcatom(ib, ig);
		}
	}

}
#endif

void wavefunc::wfcinit_k(psi::Psi<std::complex<double>>* psi_in)
{
	ModuleBase::TITLE("wavefunc","wfcinit_k");

	if(GlobalV::BASIS_TYPE=="pw") 
	{
		this->irindex = new int [GlobalC::wfcpw->fftnxy];
		GlobalC::wfcpw->getfftixy2is(this->irindex);
    #if defined(__CUDA) || defined(__UT_USE_CUDA)
    GlobalC::wfcpw->get_ig2ixyz_k();
    #endif
	}
	if(GlobalV::CALCULATION=="nscf")
	{
		return;
	}

	//---------------------------------------------------
	//  calculte the overlap <i,0 | e^{i(q+G)r} | j,R>
	//---------------------------------------------------
#ifdef __LCAO
//We are not goint to support lcao_in_paw until
//the obsolete GlobalC::hm is replaced by the 
//refactored moeules (psi, hamilt, etc.)
/*
	if((!GlobalC::chi0_hilbert.epsilon) && GlobalC::chi0_hilbert.kmesh_interpolation )    // pengfei  2016-11-23
	{
		GlobalC::chi0_hilbert.Parallel_G();    // for parallel: make sure in each core, G(all_gcars(GlobalC::sf.gcars))  are the same

		// iw1->i, iw2->j, R store the positions of the neighbor unitcells that |i,0> and |j,R> have overlaps
		R = new ModuleBase::Vector3<int>** [GlobalV::NLOCAL];
		for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
		{
			R[iw1] = new ModuleBase::Vector3<int>* [GlobalV::NLOCAL];
			for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
			{
				R[iw1][iw2] = new ModuleBase::Vector3<int>[GlobalC::chi0_hilbert.lcao_box[0]*GlobalC::chi0_hilbert.lcao_box[1]*GlobalC::chi0_hilbert.lcao_box[2]];
			}
		}

		Rmax = new int* [GlobalV::NLOCAL];
		for(int iw=0; iw<GlobalV::NLOCAL; iw++)
		{
			Rmax[iw] = new int[GlobalV::NLOCAL];
		}

		int NR; // The Max number of the overlaps for each iw1,iw2;
		NR = get_R(GlobalC::chi0_hilbert.lcao_box[0],GlobalC::chi0_hilbert.lcao_box[1],GlobalC::chi0_hilbert.lcao_box[2]);

		// store the overlap relationship to "nearest.dat"
		std::stringstream ss;
		ss << GlobalV::global_out_dir <<"nearest.dat";
		std::ofstream ofs(ss.str().c_str());
		ofs << NR << std::endl;
		std::cout <<"NR = "<<NR<<std::endl;    // Max
		for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
		{
			for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
			{
				ofs<<iw1<<"   "<<iw2<<"    "<<Rmax[iw1][iw2]<<std::endl;   // iw1, iw2, and how many overlaps between iw1 and iw2
				for(int i=0; i<Rmax[iw1][iw2]; i++)
				{
					ofs<<R[iw1][iw2][i].x<<"  "<<R[iw1][iw2][i].y<<"  "<<R[iw1][iw2][i].z<<std::endl;   // positions
				}
			}
		}
		ofs.close();

		int NG = GlobalC::chi0_hilbert.dim;  // chi0's dimension

		std::complex<double> ***wanf2_q;    // <j,0 | k+G+q>

		wanf2_q = new std::complex<double> **[GlobalC::kv.nks];
		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			wanf2_q[ik] = new std::complex<double> *[GlobalV::NLOCAL];
			for(int iw=0; iw<GlobalV::NLOCAL; iw++)
			{
				wanf2_q[ik][iw] = new std::complex<double>[npwx];
			}
		}

		std::complex<double> overlap_aux[GlobalV::NLOCAL][GlobalV::NLOCAL][NG][NR];     // <i,0 | e^{i(q+G)r} | j,R>
		std::complex<double> overlap[GlobalV::NLOCAL][GlobalV::NLOCAL][NG][NR];
		double overlap_aux_R[GlobalV::NLOCAL][GlobalV::NLOCAL][NG][NR];       //real part
		double overlap_R[GlobalV::NLOCAL][GlobalV::NLOCAL][NG][NR];
		double overlap_aux_I[GlobalV::NLOCAL][GlobalV::NLOCAL][NG][NR];       //imag part
		double overlap_I[GlobalV::NLOCAL][GlobalV::NLOCAL][NG][NR];
		
		ModuleBase::ComplexMatrix Mat;
		Mat.create(GlobalV::NLOCAL,npwx);
		ModuleBase::Vector3<double> qg; // q+G
		ModuleBase::Vector3<double> gkqg, Rcar[GlobalV::NLOCAL][GlobalV::NLOCAL][NR];  // k+G+qg, Rcartesian
		
		for(int iw1=0;iw1<GlobalV::NLOCAL; iw1++)
		{
			for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
			{
				for(int i=0; i<NR; i++)
				{
					Rcar[iw1][iw2][i].x = 0.0; Rcar[iw1][iw2][i].y = 0.0;Rcar[iw1][iw2][i].z = 0.0;
				}
			}
		}

		for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
		{
			for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
			{
				for(int i=0; i<Rmax[iw1][iw2]; i++)
				{
					Rcar[iw1][iw2][i].x = GlobalC::ucell.latvec.e11 * R[iw1][iw2][i].x
					+ GlobalC::ucell.latvec.e21 * R[iw1][iw2][i].y + GlobalC::ucell.latvec.e31 * R[iw1][iw2][i].z;
					Rcar[iw1][iw2][i].y = GlobalC::ucell.latvec.e12 * R[iw1][iw2][i].x
					+ GlobalC::ucell.latvec.e22 * R[iw1][iw2][i].y + GlobalC::ucell.latvec.e32 * R[iw1][iw2][i].z;
					Rcar[iw1][iw2][i].z = GlobalC::ucell.latvec.e13 * R[iw1][iw2][i].x
					+ GlobalC::ucell.latvec.e23 * R[iw1][iw2][i].y + GlobalC::ucell.latvec.e33 * R[iw1][iw2][i].z;
				}
			}
		}


		double arg;
		std::complex<double> phase;
		for(int iq=0; iq<GlobalC::chi0_hilbert.nq; iq++)
		{
			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					for(int ig=0; ig<NG; ig++)
					{
						ModuleBase::GlobalFunc::ZEROS( overlap_aux[iw1][iw2][ig], NR);
					}
				}
			}
		}

		// main
		for(int iq=0; iq<GlobalC::chi0_hilbert.nq; iq++)    // loop over iq
		{
			for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
			{
				for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
				{
					for(int ig=0; ig<NG; ig++)
					{
						ModuleBase::GlobalFunc::ZEROS( overlap_aux[iw1][iw2][ig], NR);
					}
				}
			}

			for(int g=0; g<NG; g++)    // loop over ig
			{
				qg.x = GlobalC::chi0_hilbert.qcar[iq][0] + GlobalC::chi0_hilbert.all_gcar[g].x;
				qg.y = GlobalC::chi0_hilbert.qcar[iq][1] + GlobalC::chi0_hilbert.all_gcar[g].y;
				qg.z = GlobalC::chi0_hilbert.qcar[iq][2] + GlobalC::chi0_hilbert.all_gcar[g].z;
				std::cout <<"qg = "<<qg.x<<" "<<qg.y<<" "<<qg.z<<std::endl;
				for(int ik=0; ik<GlobalC::kv.nks; ik++)
				{
					this->LCAO_in_pw_k_q(ik, Mat, qg);
					for(int iw=0; iw<GlobalV::NLOCAL; iw++)
					{
						for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
						{
							wanf2_q[ik][iw][ig] = Mat(iw,ig);
						}
					}
				}
				for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)    // loop over iw1
				{
					for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)  // loop over iw2
					{
						for(int ik=0;ik<GlobalC::kv.nks;ik++)         // loop over ik
						{
							for(int ig=0;ig<GlobalC::kv.ngk[ik];ig++)    // loop over ig
							{
								gkqg = GlobalC::wfcpw->getgpluskcar(ik,ig) + qg;
								for(int ir=0; ir<Rmax[iw1][iw2]; ir++)   // Rmax
								{
									arg = gkqg * Rcar[iw1][iw2][ir] * ModuleBase::TWO_PI;
									phase = std::complex<double>( cos(arg),  -sin(arg) );
									overlap_aux[iw1][iw2][g][ir] += conj(GlobalC::wf.wanf2[ik](iw1,ig))
									* wanf2_q[ik][iw2][ig] * phase/static_cast<double>(GlobalC::kv.nks);
									// Peize Lin add static_cast 2018-07-14
								}
							}
						}
					}
				}
			}

			for(int g=0; g<NG; g++)
			{
				for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
				{
					for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
					{
						for(int ir=0; ir<NR; ir++)
						{
							overlap_aux_R[iw1][iw2][g][ir] = overlap_aux[iw1][iw2][g][ir].real();
							overlap_aux_I[iw1][iw2][g][ir] = overlap_aux[iw1][iw2][g][ir].imag();
						}
					}
				}
			}

#ifdef __MPI
			MPI_Allreduce(overlap_aux_R,overlap_R,GlobalV::NLOCAL * GlobalV::NLOCAL * NG * NR,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
			MPI_Allreduce(overlap_aux_I,overlap_I,GlobalV::NLOCAL * GlobalV::NLOCAL * NG * NR,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
#endif

			for(int g=0; g<NG; g++)
			{
				for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
				{
					for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
					{
						for(int ir=0; ir<NR; ir++)
						{
							overlap[iw1][iw2][g][ir] = std::complex<double>( overlap_R[iw1][iw2][g][ir],
							overlap_I[iw1][iw2][g][ir]);
						}
					}
				}
			}

			//------------------------------
			// store the overlap in q_(iq)
			//------------------------------
			std::stringstream ss1;
			ss1 << GlobalV::global_out_dir <<"q_"<<iq;
			std::ofstream ofs1(ss1.str().c_str());
			ofs1<<NG<<std::endl;
			for(int g=0; g<NG; g++)
			{
				for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
				{
					for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
					{
						for(int ir=0; ir<Rmax[iw1][iw2]; ir++)
						{
							ofs1<<overlap[iw1][iw2][g][ir]<<"  ";
						}
						ofs1<<std::endl;
					}
				}
				// mohan update 2021-02-24
				ofs1<<std::endl; ofs1<<std::endl;
			}
			ofs1.close();
		}

		for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
		{
			for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
			{
				delete[] R[iw1][iw2];
			}
			delete[] R[iw1];
		}
		delete[] R;

		for(int iw=0; iw<GlobalV::NLOCAL; iw++)
		{
			delete[] Rmax[iw];
		}
		delete[] Rmax;

		for(int ik=0; ik<GlobalC::kv.nks; ik++)
		{
			for(int iw=0; iw<GlobalV::NLOCAL; iw++)
			{
				delete[] wanf2_q[ik][iw];
			}
			delete[] wanf2_q[ik];
		}
		delete[] wanf2_q;

	}
*/
#endif
	return;
}

//--------------------------------------------
// get the nearest unitcell positions
// that exist overlaps between two orbitals
// iw1 and iw2
//--------------------------------------------
int wavefunc::get_R(int ix, int iy, int iz)   // pengfei 2016-11-23
{
	int count;
	ModuleBase::Vector3<double> r,r1,r2;

	for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
	{
		for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
		{
			int it1 = iw2it(iw1); int ia1 = iw2ia(iw1);
			int it2 = iw2it(iw2); int ia2 = iw2ia(iw2);
			//std::cout <<"iw1= "<<iw1<<" iw2= "<<iw2<<" it1= "<<it1<<" ia1= "<<ia1<<" it2= "<<it2<<" ia2= "<<ia2<<std::endl;
			count = 0;

			for(int nx=-int(ix/2);nx<=int(ix/2);nx++)
			{
				for(int ny=-int(iy/2);ny<=int(iy/2);ny++)
				{
					for(int nz=-int(iz/2);nz<=int(iz/2);nz++)
					{
						//std::cout <<"count = "<<count<<std::endl;
						//std::cout<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
						r1.x = GlobalC::ucell.atoms[it1].tau[ia1].x * GlobalC::ucell.lat0;
						r1.y = GlobalC::ucell.atoms[it1].tau[ia1].y * GlobalC::ucell.lat0;
						r1.z = GlobalC::ucell.atoms[it1].tau[ia1].z * GlobalC::ucell.lat0;
						r2.x = (GlobalC::ucell.atoms[it2].tau[ia2].x
						+ GlobalC::ucell.latvec.e11 * nx + GlobalC::ucell.latvec.e21 * ny + GlobalC::ucell.latvec.e31 * nz) * GlobalC::ucell.lat0;
						r2.y = (GlobalC::ucell.atoms[it2].tau[ia2].y
						+ GlobalC::ucell.latvec.e12 * nx + GlobalC::ucell.latvec.e22 * ny + GlobalC::ucell.latvec.e32 * nz) * GlobalC::ucell.lat0;
						r2.z = (GlobalC::ucell.atoms[it2].tau[ia2].z
						+ GlobalC::ucell.latvec.e13 * nx + GlobalC::ucell.latvec.e23 * ny + GlobalC::ucell.latvec.e33 * nz) * GlobalC::ucell.lat0;
						r = r2 - r1;
						double distance = sqrt(r*r);

						if(distance < (GlobalC::ucell.atoms[it1].Rcut + GlobalC::ucell.atoms[it2].Rcut))
						{
							R[iw1][iw2][count].x = nx;
							R[iw1][iw2][count].y = ny;
							R[iw1][iw2][count].z = nz;
							count++;
						}
					}
				}
			}
			Rmax[iw1][iw2] = count;
		}
	}

	int NR = 0;
	for(int iw1=0;iw1<GlobalV::NLOCAL;iw1++)
	{
		for(int iw2=0;iw2<GlobalV::NLOCAL;iw2++)
		{
			if(Rmax[iw1][iw2] > NR)
			{
				NR = Rmax[iw1][iw2];
			}
		}
	}

	return NR;
}


int wavefunc::iw2it( int iw)    // pengfei 2016-11-23
{
    int ic, type;
    ic =0;
    for(int it =0; it<GlobalC::ucell.ntype; it++)
	{
        for(int ia = 0; ia<GlobalC::ucell.atoms[it].na; ia++)
        {
            for(int L=0; L<GlobalC::ucell.atoms[it].nwl+1; L++)
			{
                for(int N=0; N<GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i=0; i<(2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           type = it;
                        }
                        ic++;
					}
                }
			}
        }
	}
    return type;
}

int wavefunc::iw2ia( int iw)    // pengfei 2016-11-23
{
	int ic, na;
	ic =0;
	for(int it =0; it<GlobalC::ucell.ntype; it++)
	{
		for(int ia = 0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			for(int L=0; L<GlobalC::ucell.atoms[it].nwl+1; L++)
				for(int N=0; N<GlobalC::ucell.atoms[it].l_nchi[L]; N++)
				{
					for(int i=0; i<(2*L+1); i++)
					{
						if(ic == iw)
						{
							na = ia;
						}
						ic++;
					}
				}
		}
	}
	return na;
}

//LiuXh add a new function here,
//20180515
void wavefunc::init_after_vc(const int nks)
{
    ModuleBase::TITLE("wavefunc","init");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"npwx",this->npwx);

    assert(this->npwx > 0);
    assert(nks > 0);
    assert(GlobalV::NBANDS > 0);

    const int nks2 = nks;
	const int nbasis = this->npwx * GlobalV::NPOL;

	if((GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") || winput::out_spillage==2)
	{
		if(wanf2 != nullptr)delete[] wanf2;
		this->wanf2 = new ModuleBase::ComplexMatrix [nks2];
		for (int ik = 0; ik < nks2; ik++)
		{
			this->wanf2[ik].create(GlobalV::NLOCAL, nbasis);
		}
	}

	std::cout << " MEMORY FOR PSI (MB)  : " <<
	ModuleBase::Memory::record("wavefunc","psi",nks*GlobalV::NBANDS*nbasis,"complexmatrix") << std::endl;

    if(GlobalV::test_wf)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"psi allocation","Done");
        if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"wanf2 allocation","Done");
        }
    }

    return;
}