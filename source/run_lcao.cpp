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


    //---------------------------MD/Relax------------------
    if (GlobalV::CALCULATION == "md")
	{
		Run_MD_LCAO run_md_lcao;
		run_md_lcao.opt_cell(p_esolver);
	}
	else // cell relaxations
	{
        LOOP_cell lc;
        //keep wfc_gamma or wfc_k remaining
        lc.opt_cell(p_esolver);
    }
    //---------------------------MD/Relax------------------

	ModuleBase::timer::tick("Run_lcao","lcao_line");
    return;
}
