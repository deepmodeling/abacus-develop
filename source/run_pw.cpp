#include "run_pw.h"
#include "src_pw/global.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "src_io/optical.h"
#include "src_io/numerical_basis.h"
#include "src_io/numerical_descriptor.h"
#include "src_io/print_info.h"
#include "src_ions/Cell_PW.h"
#include "src_pw/run_md_pw.h"

Run_pw::Run_pw(){}
Run_pw::~Run_pw(){}

void Run_pw::plane_wave_line(ModuleESolver::ESolver *p_esolver)
{
    ModuleBase::TITLE("Run_pw","plane_wave_line");
	ModuleBase::timer::tick("Run_pw","plane_wave_line");

    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put GlobalV::NLOCAL and GlobalV::NBANDS as input parameters
#ifdef __LCAO
    GlobalC::ucell.setup_cell( GlobalC::ORB, GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#else
    GlobalC::ucell.setup_cell( GlobalV::global_pseudo_dir, GlobalV::stru_file, GlobalV::ofs_running);
#endif

    // mohan add 2010-10-10, just to test the symmetry of a variety
    // of systems.
    if(GlobalV::CALCULATION == "test")
    {
        Cal_Test::test_memory();
        ModuleBase::QUIT();
    }

    //------------------------------------------------------------
    //---------------------Init ESolver-------------------------
    //------------------------------------------------------------
    p_esolver->Init(INPUT, GlobalC::ucell);

    

    if(GlobalV::CALCULATION == "md")
    {
        Run_MD_PW run_md_pw;
        run_md_pw.md_cells_pw(p_esolver);
    }
    else
    {
        Cell_PW cpws;
        cpws.opt_cells_pw(p_esolver);
    }

    // cout<<"cpws SUCCESS"<<endl;


    // caoyu add 2020-11-24, mohan updat 2021-01-03
    if(GlobalV::BASIS_TYPE=="pw" && GlobalV::deepks_out_labels)
    {
        Numerical_Descriptor nc;
        nc.output_descriptor(GlobalC::wf.evc, INPUT.deepks_descriptor_lmax);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"GENERATE DESCRIPTOR FOR DEEPKS");
    }


    if(GlobalV::BASIS_TYPE=="pw" && winput::out_spillage) //xiaohui add 2013-09-01
    {
        //std::cout << "\n Output Spillage Information : " << std::endl;
        // calculate spillage value.
#ifdef __LCAO
        if ( winput::out_spillage == 3)
        {
            GlobalV::BASIS_TYPE="pw"; 
            std::cout << " NLOCAL = " << GlobalV::NLOCAL << std::endl;

            for (int ik=0; ik<GlobalC::kv.nks; ik++)
            {
                GlobalC::wf.wanf2[ik].create(GlobalV::NLOCAL, GlobalC::wf.npwx);
				if(GlobalV::BASIS_TYPE=="pw")
                {
					std::cout << " ik=" << ik + 1 << std::endl;

                    GlobalV::BASIS_TYPE="lcao_in_pw";
					GlobalC::wf.LCAO_in_pw_k(ik, GlobalC::wf.wanf2[ik]);
                    GlobalV::BASIS_TYPE="pw";
                }
            }

            //Spillage sp;
            //sp.get_both(GlobalV::NBANDS, GlobalV::NLOCAL, GlobalC::wf.wanf2, GlobalC::wf.evc);
        }
#endif

        // output overlap
        if ( winput::out_spillage <= 2 )
        {
            Numerical_Basis numerical_basis;
            numerical_basis.output_overlap(GlobalC::wf.evc);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"BASIS OVERLAP (Q and S) GENERATION.");
        }
    }


	if(Optical::opt_epsilon2 && (GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw"))
	{
		Optical opt;
		opt.cal_epsilon2(GlobalV::NBANDS);
	}

	// compute density of states
	GlobalC::en.perform_dos_pw();

	ModuleBase::timer::tick("Run_pw","plane_wave_line");
    return;
}
