#include "driver_classic.h"
#include "run_md_classic.h"
#include "../input.h"
#include "../module_base/tool_title.h"
#include "../module_base/timer.h"
#include "../module_base/global_variable.h"
#include "../module_base/memory.h"
#include "../src_io/print_info.h"

Driver_classic::Driver_classic(){}
Driver_classic::~Driver_classic(){}

void Driver_classic::init()
{
	ModuleBase::TITLE("Driver_classic", "init");

	time_t time_start = std::time(NULL);
	ModuleBase::timer::start();

	// (1) read the input parameters.
	Driver_classic::reading();

	// (2) welcome to the classic MD!
	Driver_classic::classic_world();

	// (3) output information
	time_t	time_finish= std::time(NULL);
	Print_Info::print_time(time_start, time_finish);

	// (4) close all of the running logs 
	INPUT.close_log();

	return;
}

void Driver_classic::reading(void)
{
	ModuleBase::timer::tick("Driver_classic", "reading");

	// (1) read INPUT 
	INPUT.Init( GlobalV::global_in_card );

    // (2) Print the parameters into INPUT file.
    std::stringstream ss1;
    ss1 << GlobalV::global_out_dir << GlobalV::global_in_card;
    INPUT.Print( ss1.str() );

	ModuleBase::timer::tick("Driver_classic","reading");
	return;
}

void Driver_classic::convert(UnitCell_pseudo &ucell_c)
{
    ModuleBase::TITLE("Driver_classic","convert");
	ModuleBase::timer::tick("Driver_classic","convert");

    if(INPUT.stru_file!="") GlobalV::stru_file = INPUT.stru_file;
    GlobalV::CALCULATION = INPUT.calculation;
	GlobalV::BASIS_TYPE = INPUT.basis_type;
    GlobalV::OUT_LEVEL = INPUT.out_level;
    GlobalV::SEARCH_RADIUS = max(INPUT.search_radius,(INPUT.mdp.lj_rcut+2)*ModuleBase::ANGSTROM_AU);
	GlobalV::SEARCH_PBC = INPUT.search_pbc;
    GlobalV::MD_NSTEP = INPUT.mdp.md_nstep;

	GlobalV::CAL_FORCE = INPUT.cal_force;
    GlobalV::FORCE_THR = INPUT.force_thr;

	GlobalV::CAL_STRESS = INPUT.cal_stress;
    GlobalV::STRESS_THR = INPUT.stress_thr;
    GlobalV::PRESS1 = INPUT.press1;
    GlobalV::PRESS2 = INPUT.press2;
    GlobalV::PRESS3 = INPUT.press3;

    ucell_c.latName = INPUT.latname; 
	ucell_c.ntype = INPUT.ntype;
	ucell_c.init_vel = INPUT.init_vel;

}

void Driver_classic::classic_world(void)
{
	ModuleBase::TITLE("Driver_classic", "classic_world");

	Run_MD_CLASSIC run_md_classic;

	Driver_classic::convert(run_md_classic.ucell_c);

	if(GlobalV::CALCULATION != "md")
	{
		ModuleBase::WARNING_QUIT("Driver_classic::classic_world","CALCULATION must be md!");
	}

    run_md_classic.classic_md_line();

	ModuleBase::timer::finish( GlobalV::ofs_running );

	ModuleBase::Memory::print_all( GlobalV::ofs_running ) ;

	return;
}