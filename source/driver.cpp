#include "driver.h"
#include "input.h"
#include "input_conv.h"
#include "run_pw.h"
#include "src_pw/global.h"
#ifdef __LCAO
#include "run_lcao.h"
#include "src_lcao/global_fp.h"
#endif
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "src_io/print_info.h"
#include "module_base/timer.h"
#include "module_base/memory.h"
#include "module_esolver/esolver.h"

Driver::Driver(){}


Driver::~Driver(){}


void Driver::init()
{
	ModuleBase::TITLE("Driver","init");

	time_t time_start = std::time(NULL);
	ModuleBase::timer::start();

	// (1) read the input parameters.
	this->reading();

	// (2) welcome to the atomic world!
	this->atomic_world();


	// (3) output information
	time_t	time_finish= std::time(NULL);
	Print_Info::print_time(time_start, time_finish);

	// (4) close all of the running logs 

	INPUT.close_log();

	return;
}


void Driver::reading(void)
{
	ModuleBase::timer::tick("Driver","reading");

	// (1) read INPUT 
	INPUT.Init( GlobalV::global_in_card );

	// (2) copy the variables from INPUT to each class
	Input_Conv::Convert();

	// (3) define the 'DIAGONALIZATION' world in MPI
	Parallel_Global::split_diag_world(GlobalV::DIAGO_PROC);
	Parallel_Global::split_grid_world(GlobalV::DIAGO_PROC);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DRANK",GlobalV::DRANK+1);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DSIZE",GlobalV::DSIZE);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DCOLOR",GlobalV::DCOLOR+1);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"GRANK",GlobalV::GRANK+1);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"GSIZE",GlobalV::GSIZE);

#ifdef __MPI
    // (4)  divide the GlobalV::NPROC processors into GlobalV::KPAR for k-points parallelization.
    GlobalC::Pkpoints.init_pools();
#endif

    // (5) Read in parameters about wannier functions.
    winput::Init( GlobalV::global_wannier_card );

    // (6) Print the parameters into INPUT file.
    std::stringstream ss1;
    ss1 << GlobalV::global_out_dir << GlobalV::global_in_card;
    INPUT.Print( ss1.str() );
    //ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"READING CARDS");

	ModuleBase::timer::tick("Driver","reading");
	return;
}

void Driver::atomic_world(void)
{
	ModuleBase::TITLE("Driver","atomic_world");
	//--------------------------------------------------
	// choose basis sets:
	// pw: plane wave basis set
	// lcao_in_pw: LCAO expaned by plane wave basis set
	// lcao: linear combination of atomic orbitals
	//--------------------------------------------------
	string use_ensol;
	ModuleESolver::ESolver *p_esolver;
	if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
	{
		use_ensol = "ksdft_pw"; 
		//We set it temporarily
		//Finally, we have ksdft_pw, ksdft_lcao, sdft_pw, ofdft, lj, eam, etc.
		ModuleESolver::init_esolver(p_esolver, use_ensol);
		Run_pw::plane_wave_line(p_esolver);
		ModuleESolver::clean_esolver(p_esolver);
	}
#ifdef __LCAO
	else if(GlobalV::BASIS_TYPE=="lcao")
    {
        use_ensol = "ksdft_lcao";
        ModuleESolver::init_esolver(p_esolver, use_ensol);
        Run_lcao::lcao_line(p_esolver);
        ModuleESolver::clean_esolver(p_esolver);
    }
#endif

	ModuleBase::timer::finish( GlobalV::ofs_running );

	ModuleBase::Memory::print_all( GlobalV::ofs_running ) ;

	return;
}
