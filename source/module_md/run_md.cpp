#include "run_md.h"
#include "run_md_pw.h"
#ifdef __LCAO
#include "run_md_lcao.h"
#endif
#include "run_md_classic.h"
#include "../src_pw/global.h"
#include "../input.h"
#include "../src_io/cal_test.h"
#include "../src_io/print_info.h"
#include "../module_symmetry/symmetry.h"

Run_md::Run_md(){}
Run_md::~Run_md(){}


void Run_md::md_line(void)
{
	TITLE("Run_md","md_line");
	timer::tick("Run_md","md_line");

	if(INPUT.mdp.md_potential==0)
	{
		Run_md::ai_md_line();
	}
	else
	{
		Run_md::classic_md_line();
	}

	timer::tick("Run_md","md_line");
    return;
}


void Run_md::ai_md_line(void)
{
    // Setup the unitcell.
    // improvement: a) separating the first reading of the atom_card and subsequent
    // cell relaxation. b) put NLOCAL and NBANDS as input parameters
    ucell.setup_cell( global_pseudo_dir, out, global_atom_card, ofs_running);

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

    //ucell.setup_cell( global_pseudo_dir , global_atom_card , ofs_running, NLOCAL, NBANDS);
    DONE(ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (Symmetry::symm_flag)
    {
        symm.analy_sys(ucell, out);
        DONE(ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    kv.set(symm, global_kpoint_card, NSPIN, ucell.G, ucell.latvec );
    DONE(ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info PI;
    PI.setup_parameters();

#ifdef __LCAO
    // * reading the localized orbitals/projectors 
	// * construct the interpolation tables.
	if(BASIS_TYPE=="lcao")
	{
		LOWF.orb_con.set_orb_tables(
		ofs_running,
		UOT, 
		ORB,
		ucell.ntype,
		ucell.lmax,
		INPUT.lcao_ecut,
		INPUT.lcao_dk,
		INPUT.lcao_dr,
		INPUT.lcao_rmax, 
		ucell.lat0, 
		INPUT.out_descriptor,
		INPUT.out_r_matrix,
		Exx_Abfs::Lmax,
		FORCE,
		MY_RANK);

		// * allocate H and S matrices according to computational resources
		// * set the 'trace' between local H/S and global H/S
		LM.divide_HS_in_frag(GAMMA_ONLY_LOCAL, ParaO);
	}
#endif

    // Initalize the plane wave basis set
    pw.gen_pw(ofs_running, ucell, kv);
    DONE(ofs_running,"INIT PLANEWAVE");
    cout << " UNIFORM GRID DIM     : " << pw.nx <<" * " << pw.ny <<" * "<< pw.nz << endl;
    cout << " UNIFORM GRID DIM(BIG): " << pw.nbx <<" * " << pw.nby <<" * "<< pw.nbz << endl;

    // the symmetry of a variety of systems.
    if(CALCULATION == "test")
    {
        Cal_Test::test_memory();
        QUIT();
    }

    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    Pgrid.init(pw.ncx, pw.ncy, pw.ncz, pw.nczp,
        pw.nrxx, pw.nbz, pw.bz); // mohan add 2010-07-22, update 2011-05-04

	// Inititlize the charge density.
    CHR.allocate(NSPIN, pw.nrxx, pw.ngmc);
    DONE(ofs_running,"INIT CHARGE");

	// Initializee the potential.
    pot.allocate(pw.nrxx);
    DONE(ofs_running,"INIT POTENTIAL");

	if(BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw")
	{
		Run_MD_PW run_md_pw;
        run_md_pw.md_cells_pw();
	}
#ifdef __LCAO
	else if(BASIS_TYPE=="lcao")
	{
		Run_MD_LCAO run_md_lcao;
		run_md_lcao.opt_cell();
	}
#endif

    return;
}

void Run_md::classic_md_line(void)
{
	// Setup the unitcell.
    ucell.setup_cell( global_pseudo_dir, out, global_atom_card, ofs_running);
	DONE(ofs_running, "SETUP UNITCELL");

	Print_Info PI;
    PI.setup_parameters();

	Run_MD_CLASSIC run_md_classic;
	run_md_classic.md_cells_classic();

	return;
}