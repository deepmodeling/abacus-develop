//#include "../src_pw/tools.h"
#include "ORB_control.h"
//#include "src_global/sltk_atom_arrange.h"
#include "ORB_gen_tables.h"
#include "build_st_pw.h"
#include "../src_pdiag/pdiag_double.h"

//#include "../src_pw/global.h"

ORB_control::ORB_control()
{}

ORB_control::~ORB_control()
{}

void ORB_control::set_orb_tables(
	ORB_gen_tables &OGT, 
	LCAO_Orbitals &orb,
	const int &ntype, // mohan add 2021-04-26
	const int &lmax, // mohan add 2021-04-26 
	const double &lcao_ecut_in, // mohan add 2021-04-16
	const double &lcao_dk_in, // mohan add 2021-04-16
	const double &lcao_dr_in, // mohan add 2021-04-16
	const double &lcao_rmax_in, // mohan add 2021-04-16
	const double &lat0,
	const int &out_descriptor,
	const int &out_r_matrix,
	const int &Lmax_exx,
	const int &my_rank) // mohan add 2021-04-26
{
    TITLE("ORB_control","set_orb_tables");
	timer::tick("ORB_control","set_orb_tables",'B');
    
	//=============================================================================
    // (1) FUNCTION : use 'info' to generate 'Numerical Orbital'
    // (1) RESULT : We have 'Numerical Orbital' for calculate S-table and T-table.
	//=============================================================================

	// mohan add 2021-04-16
	assert(ntype>0);
	assert(lmax>=0);
	assert(lcao_ecut_in>0.0);
	assert(lcao_dk_in>0.0);
	assert(lcao_dr_in>0.0);
	assert(lcao_rmax_in>0.0);

	// mohan add 2021-04-16
	orb.ecutwfc = lcao_ecut_in;
	orb.dk = lcao_dk_in;
	orb.dR = lcao_dr_in;
	orb.Rmax = lcao_rmax_in;
	
    orb.Read_Orbitals(
		ntype, 
		lmax, 
		out_descriptor, 
		out_r_matrix, 
		my_rank);

	if(CALCULATION=="test")
	{
		timer::tick("ORB_control","set_orb_tables",'B');
		return;
	}

    //=============================================================================
    // (2) FUNCTION : Generate Gaunt_Coefficients and S-table using OGT.init
    // 	   Must have 'Numerical Orbital' infomation
    // (2) RESULT : we have tabulated S table for use.
    //=============================================================================
    const int job0 = 3;
    // job0 :
    // 1: generate overlap table
    // 2: generate kinetic table
    // 3: generate overlap & kinetic table
    OGT.gen_tables(job0, orb, Lmax_exx, out_descriptor);
    // init lat0, in order to interpolated value from this table.

	assert(lat0>0.0);
    OGT.set_unit(lat0);


	timer::tick("ORB_control","set_orb_tables",'B');
    return;
}

void ORB_control::clear_after_ions(
	ORB_gen_tables &OGT, 
	LCAO_Orbitals &orb,
	const int &out_descriptor)
{
    TITLE("ORB_control","clear_after_ions");
    OGT.MOT.Destroy_Table(orb);
    OGT.tbeta.Destroy_Table_Beta(orb);
    
	//caoyu add 2021-03-18
    if (out_descriptor>0) 
	{
        OGT.talpha.Destroy_Table_Alpha();
    }
    return;
}
