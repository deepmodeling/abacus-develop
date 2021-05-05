#include "pseudo_atom.h"
#include "global.h"		// only rcut

pseudo_atom::pseudo_atom()
{
	r = new double[1];
	rab = new double[1];
	rho_at = new double[1];
	rho_atc = new double[1];
}

pseudo_atom::~pseudo_atom()
{
	delete[] r;
	delete[] rab;
	delete[] rho_at;
	delete[] rho_atc;
}

//---------------------------------------------------------------------
void pseudo_atom::set_pseudo_at(const Pseudopot_upf &upf)
{
	//-----------------------------------------------------------------

	set_pseudo_h(upf);
	int i, ir, j;

	// mohan 2009-12-15
	// mohan update again 2011-05-23, 
	// in order to calculate more accurate Vna.
	rcut = 15.0;//(a.u.);
	
	//if(!LOCAL_BASIS) xiaohui modify 2013-09-02 // mohan modified 2009-1-21
	if(BASIS_TYPE=="pw") //xiaohui add 2013-09-02
	{
		//rcut = winput::rcut;
	}

// remember to update here if you need it.
//	rcut = 25.0; 

	OUT(ofs_running,"PAO radial cut off (Bohr)",rcut);
	if(rcut <= 0.0)
	{
		WARNING_QUIT("pseudo_atom::set_pseudo_at","PAO rcut<=0.0");
	}

	chi.create(nchi, mesh);

	delete[] r;
	r = new double[mesh];
	assert(r != 0);
	ZEROS(r, mesh);

	delete[] rab;
	rab = new double[mesh];
	assert(rab != 0);
	ZEROS(rab, mesh);

	delete[] rho_at;
	rho_at  = new double[mesh];
	assert(rho_at != 0);
	ZEROS(rho_at,mesh);

	delete[] rho_atc;
	rho_atc = new double[mesh];
	assert(rho_atc != 0);
	ZEROS(rho_atc, mesh);

	for (i = 0;i < nchi;i++)
	{
		for (j = 0; j < mesh; j++)
		{
			chi(i, j) = upf.chi(i, j);
		}
	}

	for (i = 0;i < mesh;i++)
	{
		r  [i]     = upf.r  [i];
		rab[i]     = upf.rab[i];
		rho_at [i] = upf.rho_at [i];
	}

	if (nlcc)
	{
		for (i = 0;i < mesh;i++)
		{
			rho_atc[i] = upf.rho_atc[i];
		}
	}
	else
	{
		for (i = 0;i < upf.mesh;i++)
		{
			rho_atc[i] = 0.0;
		}
	} // end if

	bool br = false;

	msh = 0;

	for (ir = 0;ir < mesh;ir++) // do ir = 1, mesh [is]
	{
		if (r [ir] > rcut)
		{
			msh = ir + 1;
			br = true;
			break;
			//	goto 5;
		} // endif
	} // enddo

	if (br)
	{
		// force msh to be odd for simpson integration
		msh = 2 * (int)((msh + 1) / 2) - 1;	// 5
	}
	else
	{
		msh = mesh ;
	}
//	cout << " msh=" << msh << endl;
} // end subroutine set_pseudo

void pseudo_atom::print_pseudo_at(ofstream &ofs)
{
	print_pseudo_h(ofs);
	ofs << "\n pseudo_atom : ";
	ofs << "\n msh	" << msh;
//	ofs	<< "\n nchi	" << nchi;
	out.printr1_d(ofs, " r : ", r, mesh);
	out.printr1_d(ofs, " rab : ", rab, mesh);
	out.printr1_d(ofs, " rho_atc : ", rho_atc, mesh);
	out.printr1_d(ofs, " rho_at : ", rho_at, mesh);
	out.printr1_d(ofs," jchi : ", jchi, nchi);
	out.printrm(ofs, " chi : ", chi);
	ofs << "\n ----------------------";
}
