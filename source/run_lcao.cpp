#include "run_lcao.h"
#include "src_pw/global.h"
#include "input.h"
#include "src_io/optical.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "module_neighbor/sltk_atom_arrange.h"
#include "src_lcao/LOOP_cell.h"
#include "src_io/print_info.h"
#include "src_lcao/run_md_lcao.h"

Run_lcao::Run_lcao(){}
Run_lcao::~Run_lcao(){}


void Run_lcao::lcao_line(ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("Run_lcao","lcao_line");
    ModuleBase::timer::tick("Run_lcao", "lcao_line");
    
    //-----------------------init Cell--------------------------
    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put GlobalV::NLOCAL and GlobalV::NBANDS as input parameters
#ifdef __LCAO
    GlobalC::ucell.setup_cell( GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#else
    GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
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


    // the symmetry of a variety of systems.
    if (GlobalV::CALCULATION == "test")
    {
        Cal_Test::test_memory();
        ModuleBase::QUIT();
    }
    //-----------------------init Cell--------------------------


    //------------------------------------------------------------
    //---------------------Init ESolver-------------------------
    p_esolver->Init(INPUT, GlobalC::ucell);
    //------------------------------------------------------------

    //------------------init Basis_lcao----------------------
    // Init Basis should be put outside of Ensolver.
    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    ORB_control orb_con(
        GlobalV::GAMMA_ONLY_LOCAL,
        GlobalV::NLOCAL, GlobalV::NBANDS,
        GlobalV::NSPIN, GlobalV::DSIZE,
        GlobalV::NB2D, GlobalV::DCOLOR,
        GlobalV::DRANK, GlobalV::MY_RANK,
        GlobalV::CALCULATION, GlobalV::KS_SOLVER);
    
    Init_Basis_lcao(orb_con, INPUT, GlobalC::ucell);
    //------------------init Basis_lcao----------------------


    //---------------------------MD/Relax------------------
    if (GlobalV::CALCULATION == "md")
	{
		Run_MD_LCAO run_md_lcao(orb_con.ParaV);
		run_md_lcao.opt_cell(orb_con, p_esolver);
	}
	else // cell relaxations
	{
        LOOP_cell lc(orb_con.ParaV);
        //keep wfc_gamma or wfc_k remaining
        lc.opt_cell(orb_con, p_esolver);
    }
    //---------------------------MD/Relax------------------

	ModuleBase::timer::tick("Run_lcao","lcao_line");
    return;
}

void Run_lcao::Init_Basis_lcao(ORB_control& orb_con, Input& inp, UnitCell_pseudo& ucell)
{
    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    orb_con.read_orb_first(
        GlobalV::ofs_running,
        GlobalC::ORB,
        ucell.ntype,
        ucell.lmax,
        inp.lcao_ecut,
        inp.lcao_dk,
        inp.lcao_dr,
        inp.lcao_rmax,
        GlobalV::deepks_setorb,
        inp.out_mat_r,
        GlobalV::CAL_FORCE,
        GlobalV::MY_RANK);

    ucell.infoNL.setupNonlocal(
        ucell.ntype,
		ucell.atoms,
		GlobalV::ofs_running,
        GlobalC::ORB);

#ifdef __MPI   
	orb_con.set_orb_tables(
		GlobalV::ofs_running,
		GlobalC::UOT,
		GlobalC::ORB,
		ucell.lat0,
		GlobalV::deepks_setorb,
		Exx_Abfs::Lmax,
		ucell.infoNL.nprojmax,
		ucell.infoNL.nproj,
        ucell.infoNL.Beta);
#else
	int Lmax=0;
	orb_con.set_orb_tables(
		GlobalV::ofs_running,
		GlobalC::UOT,
		GlobalC::ORB,
		ucell.lat0,
		GlobalV::deepks_setorb,
		Lmax,
		ucell.infoNL.nprojmax,
	    ucell.infoNL.nproj,
        ucell.infoNL.Beta);
#endif

    if (orb_con.setup_2d)
        orb_con.setup_2d_division(GlobalV::ofs_running, GlobalV::ofs_warning);
}