#include "run_pw.h"
#include "src_pw/global.h"
#include "src_pw/energy.h"
#include "input.h"
#include "src_io/optical.h"
#include "src_io/cal_test.h"
#include "src_io/winput.h"
#include "src_io/numerical_basis.h"
#include "src_io/numerical_descriptor.h"
#include "src_io/print_info.h"
#include "module_symmetry/symmetry.h"
#include "src_ions/Cell_PW.h"

Run_pw::Run_pw(){}
Run_pw::~Run_pw(){}

void Run_pw::plane_wave_line(void)
{
    TITLE("Run_pw","plane_wave_line");
	timer::tick("Run_pw","plane_wave_line");

    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put NLOCAL and NBANDS as input parameters
    ucell.setup_cell( global_pseudo_dir, out, global_atom_card, ofs_running);
    //ucell.setup_cell( global_pseudo_dir , global_atom_card , ofs_running, NLOCAL, NBANDS);

    // setup NBANDS 
	// Yu Liu add 2021-07-03
	CHR.cal_nelec();

    // mohan add 2010-09-06
	// Yu Liu move here 2021-06-27
	// because the number of element type
	// will easily be ignored, so here
	// I warn the user again for each type.
	for(int it=0; it<ucell.ntype; it++)
	{
		xcf.which_dft(ucell.atoms[it].dft);
    }

    DONE(ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (Symmetry::symm_flag)
    {
        symm.analy_sys();
        DONE(ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    kv.set( symm, global_kpoint_card, NSPIN, ucell.G, ucell.latvec );
    DONE(ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info PI;
    PI.setup_parameters();

    // Initalize the plane wave basis set
    pw.gen_pw(ofs_running, ucell, kv);
    DONE(ofs_running,"INIT PLANEWAVE");
    cout << " UNIFORM GRID DIM     : " << pw.nx <<" * " << pw.ny <<" * "<< pw.nz << endl;
    cout << " UNIFORM GRID DIM(BIG): " << pw.nbx <<" * " << pw.nby <<" * "<< pw.nbz << endl;

    // mohan add 2010-10-10, just to test the symmetry of a variety
    // of systems.
    if(CALCULATION == "test")
    {
        Cal_Test::test_memory();
        QUIT();
    }

    // mohan add 2010-09-13
    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    Pgrid.init(pw.ncx, pw.ncy, pw.ncz, pw.nczp,
        pw.nrxx, pw.nbz, pw.bz); // mohan add 2010-07-22, update 2011-05-04

//----------------------------------------------------------
// 1 read in initial data:
//   a lattice structure:atom_species,atom_positions,lattice vector
//   b k_points
//   c pseudopotential
// 2 setup planeware basis, FFT,structure factor, ...
// 3 initialize local and nonlocal pseudopotential in G_space
// 4 initialize charge desity and warefunctios in G_space
//----------------------------------------------------------

    //=====================================
    // init charge/potential/wave functions
    //=====================================
    CHR.allocate(NSPIN, pw.nrxx, pw.ngmc);
    pot.allocate(pw.nrxx);

    Cell_PW cpws;
    cpws.opt_cells_pw();


    // caoyu add 2020-11-24, mohan updat 2021-01-03
    if(BASIS_TYPE=="pw" && INPUT.out_descriptor==1)
    {
        Numerical_Descriptor nc;
        nc.output_descriptor(wf.evc, INPUT.lmax_descriptor);
        DONE(ofs_running,"GENERATE DESCRIPTOR FOR DEEPKS");
    }


    if(BASIS_TYPE=="pw" && winput::out_spillage) //xiaohui add 2013-09-01
    {
        //cout << "\n Output Spillage Information : " << endl;
        // calculate spillage value.
#ifdef __LCAO
        if ( winput::out_spillage == 3)
        {
            BASIS_TYPE="pw"; 
            cout << " NLOCAL = " << NLOCAL << endl;

            for (int ik=0; ik<kv.nks; ik++)
            {
                wf.wanf2[ik].create(NLOCAL, wf.npwx);
				if(BASIS_TYPE=="pw")
                {
					cout << " ik=" << ik + 1 << endl;

                    BASIS_TYPE="lcao_in_pw";
					wf.LCAO_in_pw_k(ik, wf.wanf2[ik]);
                    BASIS_TYPE="pw";
                }
            }

            //Spillage sp;
            //sp.get_both(NBANDS, NLOCAL, wf.wanf2, wf.evc);
        }
#endif

        // output overlap
        if ( winput::out_spillage <= 2 )
        {
            Numerical_Basis numerical_basis;
            numerical_basis.output_overlap(wf.evc);
            DONE(ofs_running,"BASIS OVERLAP (Q and S) GENERATION.");
        }
    }


	if(Optical::opt_epsilon2 && (BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw"))
	{
		Optical opt;
		opt.cal_epsilon2(NBANDS);
	}

	// compute density of states
	en.perform_dos_pw();

	timer::tick("Run_pw","plane_wave_line");
    return;
}
