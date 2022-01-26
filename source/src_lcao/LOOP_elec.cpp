#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "input_update.h"
#include "../src_io/chi0_hilbert.h"
#include "LCAO_evolve.h"
#include "dftu.h"
//
#include "../module_neighbor/sltk_atom_arrange.h"
#include "LCAO_nnr.h"
#include "../src_io/istate_charge.h"
#include "../src_io/istate_envelope.h"
#include "ELEC_scf.h"
#include "ELEC_nscf.h"
#include "ELEC_cbands_gamma.h"
#include "ELEC_cbands_k.h"
#include "ELEC_evolve.h"
//
#include "../module_base/lapack_connector.h"
#include "../module_base/blas_connector.h"
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "../src_pw/vdwd2.h"
#include "../src_pw/vdwd3.h"
#include "../module_base/timer.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif

LOOP_elec::LOOP_elec(){
	copy_wfc_flag = true ;
}

void LOOP_elec::solve_elec_stru(const int &istep)
{
    ModuleBase::TITLE("LOOP_elec","solve_elec_stru"); 
    ModuleBase::timer::tick("LOOP_elec","solve_elec_stru"); 

	// prepare HS matrices, prepare grid integral
	this->set_matrix_grid();

	// need to add loop for k 
	// convert 2d format work to WFC_K_backup
	int zero_int=0;
	if (INPUT.tddft==1 ) 
	{
		if (!copy_wfc_flag)
		{
			GlobalC::LOWF.allocate_k_backup(GlobalC::GridT);
			GlobalC::LOWF.set_trace_aug_backup(GlobalC::GridT);
			this->store_WFC_1(zero_int);
		}
		else 
		{
			copy_wfc_flag = false ;
			this->allocate_work(zero_int);
		}
	}

	// density matrix extrapolation and prepare S,T,VNL matrices 
	this->before_solver(istep);
	// do self-interaction calculations / nscf/ tddft, etc. 
	this->solver(istep);

	// need to add loop for k 
	// convert WFC_K to 2d format work
	if (INPUT.tddft==1)
	{ 
		this->store_WFC_2(zero_int); 
	}

    ModuleBase::timer::tick("LOOP_elec","solve_elec_stru"); 
	return;
}


void LOOP_elec::set_matrix_grid(void)
{
    ModuleBase::TITLE("LOOP_elec","set_matrix_grid"); 
    ModuleBase::timer::tick("LOOP_elec","set_matrix_grid"); 

	// (1) Find adjacent atoms for each atom.
	GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
		GlobalV::ofs_running,
		GlobalV::OUT_LEVEL,
		GlobalC::ORB.get_rcutmax_Phi(), 
		GlobalC::ucell.infoNL.get_rcutmax_Beta(), 
		GlobalV::GAMMA_ONLY_LOCAL);

	atom_arrange::search(
		GlobalV::SEARCH_PBC,
		GlobalV::ofs_running,
		GlobalC::GridD, 
		GlobalC::ucell, 
		GlobalV::SEARCH_RADIUS, 
		GlobalV::test_atom_input);

	//ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"SEARCH ADJACENT ATOMS");

	// (3) Periodic condition search for each grid.
	GlobalC::GridT.set_pbc_grid(
			GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz,
			GlobalC::pw.bx, GlobalC::pw.by, GlobalC::pw.bz,
			GlobalC::pw.nbx, GlobalC::pw.nby, GlobalC::pw.nbz,
			GlobalC::pw.nbxx, GlobalC::pw.nbzp_start, GlobalC::pw.nbzp);

	// (2) If k point is used here, allocate HlocR after atom_arrange.
	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		// For each atom, calculate the adjacent atoms in different cells
		// and allocate the space for H(R) and S(R).
		GlobalC::LNNR.cal_nnr();
		GlobalC::LM.allocate_HS_R(GlobalC::LNNR.nnr);
#ifdef __DEEPKS
		GlobalC::ld.allocate_V_deltaR(GlobalC::LNNR.nnr);
#endif

		// need to first calculae lgd.
		// using GlobalC::GridT.init.
		GlobalC::LNNR.cal_nnrg(GlobalC::GridT);
	}

    ModuleBase::timer::tick("LOOP_elec","set_matrix_grid"); 
	return;
}


void LOOP_elec::before_solver(const int &istep)
{
    ModuleBase::TITLE("LOOP_elec","before_solver"); 
    ModuleBase::timer::tick("LOOP_elec","before_solver"); 

	// set the augmented orbitals index.
	// after ParaO and GridT, 
	// this information is used to calculate
	// the force.
	GlobalC::LOWF.set_trace_aug(GlobalC::GridT); //LiuXh modify 2021-09-06, clear memory, WFC_GAMMA_aug not used now

	// init density kernel and wave functions.
	GlobalC::LOC.allocate_dm_wfc(GlobalC::GridT);

	//======================================
	// do the charge extrapolation before the density matrix is regenerated.
	// mohan add 2011-04-08
	// because once atoms are moving out of this processor,
	// the density matrix will not std::map the new atomic configuration,
	//======================================
	// THIS IS A BUG, BECAUSE THE INDEX GlobalC::GridT.trace_lo
	// HAS BEEN REGENERATED, SO WE NEED TO
	// REALLOCATE DENSITY MATRIX FIRST, THEN READ IN DENSITY MATRIX,
	// AND USE DENSITY MATRIX TO DO RHO GlobalV::CALCULATION.-- mohan 2013-03-31
	//======================================
	if(GlobalC::pot.extra_pot=="dm" && istep>1)//xiaohui modify 2015-02-01
	{
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::pw.nrxx);
			std::stringstream ssd;
			ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM" ;
			// reading density matrix,
			GlobalC::LOC.read_dm(is, ssd.str() );
		}

		// calculate the charge density
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			GlobalC::UHM.GG.cal_rho(GlobalC::LOC.DM);
		}
		else
		{
			GlobalC::UHM.GK.cal_rho_k();
		}

		// renormalize the charge density
		GlobalC::CHR.renormalize_rho();

		// initialize the potential
		GlobalC::pot.init_pot( istep-1, GlobalC::pw.strucFac );
	}


	// (9) compute S, T, Vnl, Vna matrix.
    GlobalC::UHM.set_lcao_matrices();

#ifdef __DEEPKS
    //for each ionic step, the overlap <psi|alpha> must be rebuilt
    //since it depends on ionic positions
    if (GlobalV::out_descriptor)
    {
		//build and save <psi(0)|alpha(R)> at beginning
        GlobalC::ld.build_psialpha(GlobalV::FORCE,
			GlobalC::ucell,
			GlobalC::ORB,
			GlobalC::GridD,
			GlobalC::ParaO,
			GlobalC::UOT);
    }
#endif

    ModuleBase::timer::tick("LOOP_elec","before_solver"); 
	return;
}

void LOOP_elec::solver(const int &istep)
{
    ModuleBase::TITLE("LOOP_elec","solver"); 
    ModuleBase::timer::tick("LOOP_elec","solver"); 

	// self consistent calculations for electronic ground state
	if (GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md"
			|| GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax") //pengfei 2014-10-13
	{
		//Peize Lin add 2016-12-03
		switch(GlobalC::exx_lcao.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				GlobalC::exx_lcao.cal_exx_ions();
				break;
			case Exx_Global::Hybrid_Type::No:
			case Exx_Global::Hybrid_Type::Generate_Matrix:
				break;
			default:
				throw std::invalid_argument(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}

		// No exx
		if( Exx_Global::Hybrid_Type::No==GlobalC::exx_global.info.hybrid_type  )
		{
			ELEC_scf es;
			es.scf(istep-1);
		}
		else if( Exx_Global::Hybrid_Type::Generate_Matrix == GlobalC::exx_global.info.hybrid_type )
		{
			Exx_Opt_Orb exx_opt_orb;
			exx_opt_orb.generate_matrix();
		}
		else    // Peize Lin add 2016-12-03
		{
			ELEC_scf es;
			es.scf(istep-1);
			if( GlobalC::exx_global.info.separate_loop )
			{
				for( size_t hybrid_step=0; hybrid_step!=GlobalC::exx_global.info.hybrid_step; ++hybrid_step )
				{
					GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);
					GlobalC::exx_lcao.cal_exx_elec();
					
					ELEC_scf es;
					es.scf(istep-1);
					if(ELEC_scf::iter==1)     // exx converge
					{
						break;
					}
				}
			}
			else
			{
				GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);

				ELEC_scf es;
				es.scf(istep-1);
				
			}
		}
	}
	else if (GlobalV::CALCULATION=="nscf")
	{
		ELEC_nscf::nscf(GlobalC::UHM);
	}
	else if (GlobalV::CALCULATION=="istate")
	{
		IState_Charge ISC;
		ISC.begin();
	}
	else if (GlobalV::CALCULATION=="ienvelope")
	{
		IState_Envelope IEP;
		IEP.begin();
	}
	else
	{
		ModuleBase::WARNING_QUIT("LOOP_elec::solver","CALCULATION type not supported");
	}

    ModuleBase::timer::tick("LOOP_elec","solver"); 
	return;
}

inline int globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localIndex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localIndex%nblk;
    return gIndex;
}

inline int work2WFC(
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	std::complex<double>* work,
	std::complex<double>** WFC,
	std::complex<double>** WFCAUG)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                WFC[igcol][mu_local]=work[j*naroc[0]+i];
            }
	        int mu_aug=GlobalC::LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                WFCAUG[igcol][mu_aug]=work[j*naroc[0]+i];
            }
        }
    }
    return 0;
}

inline int WFC2work(
	int naroc[2],
	int nb,
	int dim0,
	int dim1,
	int iprow,
	int ipcol,
	std::complex<double>* work,
	std::complex<double>** WFC,
	std::complex<double>** WFCAUG)
{
    for(int j=0; j<naroc[1]; ++j)
    {
        int igcol=globalIndex(j, nb, dim1, ipcol);
        if(igcol>=GlobalV::NBANDS) continue;
        for(int i=0; i<naroc[0]; ++i)
        {
            int igrow=globalIndex(i, nb, dim0, iprow);
	        int mu_local=GlobalC::GridT.trace_lo[igrow];
            if(mu_local>=0)
            {
                work[j*naroc[0]+i]=WFC[igcol][mu_local];
            }
	        int mu_aug=GlobalC::LOWF.trace_aug[igrow];
            if(mu_aug>=0)
            {
                work[j*naroc[0]+i]=WFCAUG[igcol][mu_aug];
            }
        }
    }
    return 0;
}

void LOOP_elec::allocate_work(const int &ik)
{
	int ncol,nrow;
	const int inc = 1;
	int info;
	ncol = GlobalC::ParaO.ncol;
	nrow = GlobalC::ParaO.nrow;
	int nprocs, myid;
	long maxnloc,nloc; // maximum number of elements in local matrix
	nloc=GlobalC::ParaO.nloc;
	MPI_Comm comm_2D_0;
	MPI_Comm comm_2D;
	comm_2D_0 = MPI_COMM_WORLD;
	int dim0,dim1;
	dim0 = (int)sqrt((double)GlobalV::DSIZE); 
	while (GlobalV::DSIZE%dim0!=0)
	{
		dim0 = dim0 - 1;
	}
	assert(dim0 > 0);
	dim1=GlobalV::DSIZE/dim0;
	int dim[2];
	dim[0]=dim0;
	dim[1]=dim1;
	int period[2]={1,1};
	int reorder=0;
	MPI_Cart_create(comm_2D_0,2,dim,period,reorder,&comm_2D);
    MPI_Status status;
    MPI_Comm_size(comm_2D, &nprocs);
    MPI_Comm_rank(comm_2D, &myid);
    MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
    MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
	int naroc[2]; // maximum number of row or column
	int nb;
	nb = GlobalC::ParaO.nb;
	this->work1=new std::complex<double>[maxnloc];
	ModuleBase::GlobalFunc::ZEROS(work1,maxnloc);
	this->work2=new std::complex<double>[maxnloc];
	ModuleBase::GlobalFunc::ZEROS(work2,maxnloc);
	// need to add deletion of work in the last step
}

void LOOP_elec::store_WFC_1(const int &ik)
{
	int ncol,nrow;
	const int inc = 1;
	int info;
	ncol = GlobalC::ParaO.ncol;
	nrow = GlobalC::ParaO.nrow;
	int nprocs, myid;
	long maxnloc,nloc; // maximum number of elements in local matrix
	nloc=GlobalC::ParaO.nloc;
	MPI_Comm comm_2D_0;
	MPI_Comm comm_2D;
	comm_2D_0 = MPI_COMM_WORLD;
	int dim0,dim1;
	dim0 = (int)sqrt((double)GlobalV::DSIZE); 
	while (GlobalV::DSIZE%dim0!=0)
	{
		dim0 = dim0 - 1;
	}
	assert(dim0 > 0);
	dim1=GlobalV::DSIZE/dim0;
	int dim[2];
	dim[0]=dim0;
	dim[1]=dim1;
	int period[2]={1,1};
	int reorder=0;
	MPI_Cart_create(comm_2D_0,2,dim,period,reorder,&comm_2D);
    MPI_Status status;
    MPI_Comm_size(comm_2D, &nprocs);
    MPI_Comm_rank(comm_2D, &myid);
    MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
    MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
	int naroc[2]; // maximum number of row or column
	int nb;
	nb = GlobalC::ParaO.nb;
	for(int iprow=0; iprow<dim0; ++iprow)
    {
        for(int ipcol=0; ipcol<dim1; ++ipcol)
        {
            const int coord[2]={iprow, ipcol};
            int src_rank;
            MPI_Cart_rank(comm_2D, coord, &src_rank);
            if(myid==src_rank)
        	{
				zcopy_(&nloc, work1, &inc, work2, &inc);
                naroc[0]=nrow;
                naroc[1]=ncol;
			}				
			info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, comm_2D);
            info=MPI_Bcast(work2, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, comm_2D);
			info=work2WFC(naroc, nb,dim0, dim1, iprow, ipcol, work2, GlobalC::LOWF.WFC_K_backup[ik], GlobalC::LOWF.WFC_K_aug_backup[ik]);
		}
	}
}

void LOOP_elec::store_WFC_2(const int &ik)
{
	int ncol,nrow;
	const int inc = 1;
	int info;
	ncol = GlobalC::ParaO.ncol;
	nrow = GlobalC::ParaO.nrow;
	int nprocs, myid;
	long maxnloc,nloc; // maximum number of elements in local matrix
	nloc=GlobalC::ParaO.nloc;
	MPI_Comm comm_2D_0;
	MPI_Comm comm_2D;
	comm_2D_0 = MPI_COMM_WORLD;
	int dim0,dim1;
	dim0 = (int)sqrt((double)GlobalV::DSIZE); 
	while (GlobalV::DSIZE%dim0!=0)
	{
		dim0 = dim0 - 1;
	}
	assert(dim0 > 0);
	dim1=GlobalV::DSIZE/dim0;
	int dim[2];
	dim[0]=dim0;
	dim[1]=dim1;
	int period[2]={1,1};
	int reorder=0;
	MPI_Cart_create(comm_2D_0,2,dim,period,reorder,&comm_2D);
    MPI_Status status;
    MPI_Comm_size(comm_2D, &nprocs);
    MPI_Comm_rank(comm_2D, &myid);
    MPI_Reduce(&nloc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, comm_2D);
    MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, comm_2D);
	int naroc[2]; // maximum number of row or column
	int nb;
	nb = GlobalC::ParaO.nb;
	for(int iprow=0; iprow<dim0; ++iprow)
    {
        for(int ipcol=0; ipcol<dim1; ++ipcol)
        {
            const int coord[2]={iprow, ipcol};
            int src_rank;
            MPI_Cart_rank(comm_2D, coord, &src_rank);
            if(myid==src_rank)
        	{
                naroc[0]=nrow;
                naroc[1]=ncol;
				info=WFC2work(naroc, nb, dim0, dim1, iprow, ipcol, work1, GlobalC::LOWF.WFC_K[ik], GlobalC::LOWF.WFC_K_aug[ik]);
			}
		}
	}
}
