#include "run_lcao.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_io/optical.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "module_neighbor/sltk_atom_arrange.h"
#include "src_lcao/LOOP_cell.h"
#include "src_io/print_info.h"
#include "module_symmetry/symmetry.h"
#include "src_lcao/run_md_lcao.h"

Run_lcao::Run_lcao(){}
Run_lcao::~Run_lcao(){}


void Run_lcao::lcao_line(void)
{
    ModuleBase::TITLE("Run_lcao","lcao_line");
	ModuleBase::timer::tick("Run_lcao","lcao_line");

    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put GlobalV::NLOCAL and GlobalV::NBANDS as input parameters
#ifdef __LCAO
    GlobalC::ucell.setup_cell( GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalC::out, GlobalV::global_atom_card, GlobalV::ofs_running);
#else
    GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir, GlobalC::out, GlobalV::global_atom_card, GlobalV::ofs_running);
#endif
	if(INPUT.test_just_neighbor)
	{
		//test_search_neighbor();
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
			GlobalV::test_atom_input,
			INPUT.test_just_neighbor);
	}
	// setup GlobalV::NBANDS
	// Yu Liu add 2021-07-03
	GlobalC::CHR.cal_nelec();

	// mohan add 2010-09-06
	// Yu Liu move here 2021-06-27
	// because the number of element type
	// will easily be ignored, so here
	// I warn the user again for each type.
	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		GlobalC::xcf.which_dft(GlobalC::ucell.atoms[it].dft);
	}

    //GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir , GlobalV::global_atom_card , GlobalV::ofs_running, GlobalV::NLOCAL, GlobalV::NBANDS);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(GlobalC::ucell, GlobalC::out, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set(GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, GlobalC::ucell.G, GlobalC::ucell.latvec );
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(GlobalC::ucell, GlobalC::kv, GlobalC::xcf);

    // * reading the localized orbitals/projectors
	// * construct the interpolation tables.

	GlobalC::LOWF.orb_con.set_orb_tables(
		GlobalV::ofs_running,
		GlobalC::UOT,
		GlobalC::ORB,
		GlobalC::ucell.ntype,
		GlobalC::ucell.lmax,
		INPUT.lcao_ecut,
		INPUT.lcao_dk,
		INPUT.lcao_dr,
		INPUT.lcao_rmax,
		GlobalC::ucell.lat0,
		INPUT.out_descriptor,
		INPUT.out_r_matrix,
		Exx_Abfs::Lmax,
		GlobalV::FORCE,
		GlobalV::MY_RANK,
		GlobalC::ucell.infoNL.nprojmax,
		GlobalC::ucell.infoNL.nproj,
		GlobalC::ucell.infoNL.Beta);

	// * allocate H and S matrices according to computational resources
	// * set the 'trace' between local H/S and global H/S
	GlobalC::LM.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, GlobalC::ParaO);

//--------------------------------------
// cell relaxation should begin here
//--------------------------------------

    // Initalize the plane wave basis set
    GlobalC::pw.gen_pw(GlobalV::ofs_running, GlobalC::ucell, GlobalC::kv);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT PLANEWAVE");
    std::cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx <<" * " << GlobalC::pw.ny <<" * "<< GlobalC::pw.nz << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx <<" * " << GlobalC::pw.nby <<" * "<< GlobalC::pw.nbz << std::endl;

    // the symmetry of a variety of systems.
    if(GlobalV::CALCULATION == "test")
    {
        Cal_Test::test_memory();
        ModuleBase::QUIT();
    }

    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz, GlobalC::pw.nczp,
        GlobalC::pw.nrxx, GlobalC::pw.nbz, GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04


	// Inititlize the charge density.
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::pw.nrxx, GlobalC::pw.ngmc);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT CHARGE");

	// Initializee the potential.
    GlobalC::pot.allocate(GlobalC::pw.nrxx);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT POTENTIAL");


	// Peize Lin add 2018-11-30
	if(GlobalV::CALCULATION=="nscf")
	{
		switch(GlobalC::exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				GlobalC::exx_global.info.set_xcfunc(GlobalC::xcf);
				break;
		}
	}

	if(GlobalV::CALCULATION=="md")
	{
		Run_MD_LCAO run_md_lcao;
		run_md_lcao.opt_cell();
	}
	else // cell relaxations
	{
		LOOP_cell lc;
		lc.opt_cell();

		GlobalC::en.perform_dos();
	}

	ModuleBase::timer::tick("Run_lcao","lcao_line");
    return;
}
