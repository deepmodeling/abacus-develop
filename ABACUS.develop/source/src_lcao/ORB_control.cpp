#include "../src_pw/tools.h"
#include "../src_pw/global.h"
#include "ORB_control.h"
#include "src_global/sltk_atom_arrange.h"
#include "ORB_gen_tables.h"
#include "build_st_pw.h"

ORB_control::ORB_control()
{}

ORB_control::~ORB_control()
{
	if(test_deconstructor)
	{
		cout << " ~ORB_control()" << endl;
	}
}

void ORB_control::set_orb_tables(void)
{
    TITLE("ORB_control","set_orb_tables");
	timer::tick("ORB_control","set_orb_tables",'B');
    
	//=============================================================================
    // (1) FUNCTION : use 'info' to generate 'Numerical Orbital'
    // (1) RESULT : We have 'Numerical Orbital' for calculate S-table and T-table.
	//=============================================================================
    ORB.Read_Orbitals();

	if(CALCULATION=="test")
	{
		timer::tick("ORB_control","set_orb_tables",'B');
		return;
	}

    //=============================================================================
    // (2) FUNCTION : Generate Gaunt_Coefficients and S-table using UOT.init
    // 	   Must have 'Numerical Orbital' infomation
    // (2) RESULT : we have tabulated S table for use.
    //=============================================================================
    const int job0 = 3;
    // job0 :
    // 1: generate overlap table
    // 2: generate kinetic table
    // 3: generate overlap & kinetic table
    UOT.gen_tables(job0);
    // init lat0, in order to interpolated value from this table.
    UOT.set_unit(ucell.lat0);


	timer::tick("ORB_control","set_orb_tables",'B');
    return;
}

void ORB_control::clear_after_ions(void)
{
    TITLE("ORB_control","clear_after_ions");
    UOT.MOT.Destroy_Table();
    UOT.tbeta.Destroy_Table_Beta();
    //caoyu add 2021-03-18
    if (INPUT.out_descriptor && BASIS_TYPE == "lcao") {
        UOT.talpha.Destroy_Table_Alpha();
    }
    return;
}
