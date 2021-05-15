#include "sltk_atom_arrange.h"
#include "sltk_atom_input.h"
#include "sltk_grid.h"
#include "sltk_grid_driver.h"

// update the followig two includes in near future 
#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h" // mohan add 2021-01-30


atom_arrange::atom_arrange()
{
}

atom_arrange::~atom_arrange()
{
}

void atom_arrange::set_sr_NL(void)
{
	TITLE("atom_arrange","set_sr_NL");
	double longest_orb_rcut = 0.0;
	for(int it=0; it<ucell.ntype; it++)
	{
		longest_orb_rcut = std::max(longest_orb_rcut, ORB.Phi[it].getRcut() );	
	}
//	cout << " LONGEST ORB RCUT     : " << longest_orb_rcut << endl;

	if(OUT_LEVEL != "m") //xiaohui add 'OUT_LEVEL', 2015-09-16
	{
		ofs_running << "\n\n\n\n";
		ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
		ofs_running << " |                                                                    |" << endl;
		ofs_running << " | Search adjacent atoms:                                             |" << endl;
		ofs_running << " | Set the adjacent atoms for each atom and set the periodic boundary |" << endl;
		ofs_running << " | condition for the atoms on real space FFT grid. For k-dependent    |" << endl;  
		ofs_running << " | algorithm, we also need to set the sparse H and S matrix element   |" << endl;
		ofs_running << " | for each atom.                                                     |" << endl; 
		ofs_running << " |                                                                    |" << endl;
		ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
		ofs_running << "\n\n\n\n";
	}

	
	//xiaohui add 'OUT_LEVEL' line, 2015-09-16
	if(OUT_LEVEL != "m") ofs_running << "\n SETUP SEARCHING RADIUS FOR PROGRAM TO SEARCH ADJACENT ATOMS" << endl;
	if(OUT_LEVEL != "m") ofs_running << setprecision(3);
	if(OUT_LEVEL != "m") OUT(ofs_running,"longest orb rcut (Bohr)",longest_orb_rcut);

	double longest_nl_proj_rcut = 0.0;
	for(int it=0; it<ucell.ntype; it++)
	{
		longest_nl_proj_rcut = std::max(longest_nl_proj_rcut, ORB.Beta[it].get_rcut_max());
	}
//	cout << " LONGEST NL PROJ RCUT : " << longest_nl_proj_rcut << endl;
	if(OUT_LEVEL != "m") OUT(ofs_running,"longest nonlocal projector rcut (Bohr)", longest_nl_proj_rcut);

	// check in use_overlap_matrix, 
	double sr = 0.0;
	if(GAMMA_ONLY_LOCAL)
	{
		sr = 2 * longest_orb_rcut + 0.01;
	}
	else
	{
		sr = 2 * (longest_orb_rcut + longest_nl_proj_rcut) + 0.01; // 0.01 is added to make safe.
		//sr = 2 * longest_orb_rcut + 0.01;
	}

	// if use build_Nonlocal_mu (not GAMMA_ONLY_LOCAL) use 2*longest_orb_rcut
	// if use build_Nonlocal_beta ( K-point used ) use 2 * (longest_orb_rcut + longest_nl_proj_rcut) 
	SEARCH_RADIUS = sr;			
//	cout << " SEARCH RADIUS (BOHR) : " << SEARCH_RADIUS << endl;
//	OUT(ofs_running,"search radius (Bohr)", SEARCH_RADIUS);
}

// mohan update 2011-03-10
void atom_arrange::set_sr_OV(void)
{
	TITLE("atom_arrange","set_sr_OV");
	double longest_orb_rcut = 0.0;
	for(int it=0; it<ucell.ntype; it++)
	{
		longest_orb_rcut = std::max(longest_orb_rcut, ORB.Phi[it].getRcut() );	
	}
//	cout << " LONGEST ORB RCUT     : " << longest_orb_rcut << endl;
	OUT(ofs_running,"longest orb rcut (Bohr)",longest_orb_rcut);
	double sr = 2 * longest_orb_rcut + 0.01;
	// if use build_Nonlocal_mu (not GAMMA_ONLY_LOCAL) use 2*longest_orb_rcut
	// if use build_Nonlocal_beta ( K-point used ) use 2 * (longest_orb_rcut + longest_nl_proj_rcut) 
	SEARCH_RADIUS = sr;			
//	cout << " SEARCH RADIUS (BOHR) : " << SEARCH_RADIUS << endl;
	return;
}

void atom_arrange::search( const double &search_radius_bohr)
{
	TITLE("atom_arrange", "search");
	timer::tick("atom_arrange","search");

	assert( search_radius_bohr > 0.0 );

//	OUT(ofs_running,"Atom coordinates reading from",global_atom_card);
//	OUT(ofs_running,"The coordinate type",ucell.Coordinate);
//	OUT(ofs_running,"Use cartesian(unit:lat0) coordinate","TRUE");
	if(OUT_LEVEL != "m") OUT(ofs_running,"searching radius is (Bohr))", search_radius_bohr);
	if(OUT_LEVEL != "m") OUT(ofs_running,"searching radius unit is (Bohr))",ucell.lat0);

	assert(ucell.nat > 0);
	//=============================
	// Initial Atom information
	//=============================

	const double radius_lat0unit = search_radius_bohr / ucell.lat0;

	Atom_input at(ucell.nat, ucell.ntype, SEARCH_PBC, radius_lat0unit);
	//===========================================
	// Print important information in Atom_input
	//===========================================
//	at.print(cout);
//	at.print_xyz_format("1.xyz");
	//=========================================
	// Construct Grid , Cells , Adjacent atoms
	//=========================================
	GridD.init(at);

	// test the adjacent atoms and the box.
	//ofs_running << " " << setw(5) << "Type" << setw(5) << "Atom" << setw(8) << "AdjNum" << endl;
	for (int it = 0;it < ucell.ntype;it++)
	{
		for (int ia = 0;ia < ucell.atoms[it].na;ia++)
		{
	//		GridD.Find_atom(ucell.atoms[it].tau[ia]);
			
	//		ofs_running << " " << setw(5) << it << setw(5) << ia << setw(8) << GridD.getAdjacentNum()+1 << endl;
			/*
			for(int ad=0; ad < GridD.getAdjacentNum()+1; ad++)
			{
				Vector3<double> tau = GridD.getAdjacentTau(ad);
				Vector3<int> box = GridD.getBox(ad);
				cout << setw(8) << tau.x << setw(8) << tau.y << setw(8) << tau.z 
				<< setw(8) << box.x << setw(8) << box.y << setw(8) << box.z << endl;
			}
			*/
		}
	}
	
	timer::tick("atom_arrange","search");
	return;
}


//2015-05-07
void atom_arrange::delete_vector(const double &search_radius_bohr)
{
	const double radius_lat0unit2 = search_radius_bohr / ucell.lat0;

	Atom_input at2(ucell.nat, ucell.ntype, SEARCH_PBC, radius_lat0unit2);

	GridD.delete_vector(at2);

	if (GridD.init_cell_flag)
	{
		for (int i = 0;i < GridD.dx;i++)
		{
			for (int j = 0;j < GridD.dy;j++)
			{
				delete[] GridD.Cell[i][j];
			}
		}

		for (int i = 0;i < GridD.dx;i++)
		{
			delete[] GridD.Cell[i];
		}

		delete[] GridD.Cell;
	}
}
