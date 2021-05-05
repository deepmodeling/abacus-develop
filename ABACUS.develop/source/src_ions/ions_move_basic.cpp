#include "ions_move_basic.h"
#include "../src_pw/global.h"

int Ions_Move_Basic::dim=0;
bool Ions_Move_Basic::converged=true;
double Ions_Move_Basic::largest_grad=0.0;
int Ions_Move_Basic::update_iter=0;
int Ions_Move_Basic::istep=0;

double Ions_Move_Basic::ediff=0.0;
double Ions_Move_Basic::etot=0.0;
double Ions_Move_Basic::etot_p=0.0;

double Ions_Move_Basic::trust_radius=0.0;
double Ions_Move_Basic::trust_radius_old=0.0;
double Ions_Move_Basic::trust_radius_max = -1.0; // default is 0.8
double Ions_Move_Basic::trust_radius_min = -1.0; // default is 1e-5
double Ions_Move_Basic::trust_radius_ini = -1.0; // default is 0.5
double Ions_Move_Basic::best_xxx = 1.0;

int Ions_Move_Basic::out_stru=0;

void Ions_Move_Basic::setup_gradient(double* pos, double *grad, const matrix &force)
{
	TITLE("Ions_Move_Basic","setup_gradient");
	
	assert(ucell.ntype>0);
	assert(pos!=NULL);
	assert(grad!=NULL);
	assert(dim == 3*ucell.nat);

	ZEROS(pos, dim);
	ZEROS(grad, dim);

	// (1) init gradient
	// the unit of pos: Bohr.
	// the unit of force: Ry/Bohr.
	// the unit of gradient: 
	int iat=0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		Atom* atom = &ucell.atoms[it];
		for(int ia =0;ia< ucell.atoms[it].na;ia++)
		{	
			pos[3*iat  ] = atom->tau[ia].x*ucell.lat0;
			pos[3*iat+1] = atom->tau[ia].y*ucell.lat0;
			pos[3*iat+2] = atom->tau[ia].z*ucell.lat0;

			// mohan remove mbl constrain 2010-04-26
			// mohan add mbl constrain 2010-07-11
			// mohan add ucell.lat0 2010-07-27
			if(atom->mbl[ia].x == 1)
			{
				grad[3*iat  ] = -force(iat, 0)*ucell.lat0;
				//this->grad[3*iat  ] = -force(iat, 0);
			}
			if(atom->mbl[ia].y == 1)
			{
				grad[3*iat+1] = -force(iat, 1)*ucell.lat0;
			}
			if(atom->mbl[ia].z == 1)
			{
				grad[3*iat+2] = -force(iat, 2)*ucell.lat0;
				//cout << " grad=" << grad[3*iat+2] << endl;
			}
			++iat;
		}
	}

	return;
}

void Ions_Move_Basic::move_atoms(double *move, double *pos)
{
	TITLE("Ions_Move_Basic","move_atoms");

	assert(move!=NULL);
	assert(pos!=NULL);

	// unit: Bohr
	int iat=0;

	//------------------------
	// for test only
	//------------------------
	if(test_ion_dynamics)
	{
		iat=0;
		ofs_running << "\n movement of ions (unit is Bohr) : " << endl;
		ofs_running << " " << setw(12) << "Atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
		for(int it = 0;it < ucell.ntype;it++)
		{
			for(int ia =0;ia< ucell.atoms[it].na;ia++)
			{
				stringstream ss;
				ss << "move_" << ucell.atoms[it].label << ia+1;
				ofs_running << " " 
					<< setw(12) << ss.str().c_str()
					<< setw(15) << move[3*iat+0] 
					<< setw(15) << move[3*iat+1] 
					<< setw(15) << move[3*iat+2] << endl;
				iat++;
			}
		}
		assert( iat == ucell.nat );
	}

	iat = 0;
	double move_threshold = 1.0e-10;
	for(int it = 0;it < ucell.ntype;it++)
	{
		Atom* atom = &ucell.atoms[it];
		for(int ia =0;ia< atom->na;ia++)
		{
			// mohan add 2010-08-06
			// otherwise, there might be bug for
			// sltk_grid, on system CO when C
			// atom is put on (0,0,0)
			for(int i=0; i<3; i++)
			{
				if( abs(move[3*iat+i]) < move_threshold )
				{
					move[3*iat+i] = 0.0;
				}
			}
		
			// mohan modify 2010-04-26
			if(atom->mbl[ia].x!=0)
			{
				atom->tau[ia].x = (move[3*iat]+pos[3*iat])/ucell.lat0;
			}
			if(atom->mbl[ia].y!=0)
			{
				atom->tau[ia].y = (move[3*iat+1]+pos[3*iat+1])/ucell.lat0;
			}
			if(atom->mbl[ia].z!=0)
			{
				atom->tau[ia].z = (move[3*iat+2]+pos[3*iat+2])/ucell.lat0;
			}

			// the direct coordinates also need to be updated.
			atom->taud[ia] = atom->tau[ia] * ucell.GT;
//			cout << " tau=" << atom->tau[ia].x << " " << atom->tau[ia].y << " " << atom->tau[ia].z << endl;
			iat++;
		}
	}
	assert(iat == ucell.nat);

	//----------------------------------------------
	// because of the periodic boundary condition
	// we need to adjust the atom positions,
	// first adjust direct coordinates,
	// then update them into cartesian coordinates,
	//----------------------------------------------
	for(int it=0; it<ucell.ntype; it++)
	{
		Atom* atom = &ucell.atoms[it];
		for(int ia=0; ia<atom->na; ia++)
		{
			// mohan update 2011-03-21
			if(atom->taud[ia].x<0) atom->taud[ia].x += 1.0;
			if(atom->taud[ia].y<0) atom->taud[ia].y += 1.0;
			if(atom->taud[ia].z<0) atom->taud[ia].z += 1.0;
			if(atom->taud[ia].x>=1.0) atom->taud[ia].x -= 1.0;
			if(atom->taud[ia].y>=1.0) atom->taud[ia].y -= 1.0;
			if(atom->taud[ia].z>=1.0) atom->taud[ia].z -= 1.0;

			if(atom->taud[ia].x<0 || atom->taud[ia].y<0
				|| atom->taud[ia].z<0 ||
				atom->taud[ia].x>=1.0 ||
				atom->taud[ia].y>=1.0 ||
				atom->taud[ia].z>=1.0)
			{
				ofs_warning << " it=" << it+1 << " ia=" << ia+1 << endl;
				ofs_warning << "d=" << atom->taud[ia].x << " " << 
				atom->taud[ia].y << " " << atom->taud[ia].z << endl;
				WARNING_QUIT("Ions_Move_Basic::move_ions","the movement of atom is larger than the length of cell.");
			}

			atom->tau[ia] = atom->taud[ia] * ucell.latvec;
		}
	}
//2015-09-16
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i=0;i<ucell.ntype;i++)
    {
        ucell.atoms[i].bcast_atom(); // bcast tau array
    }
#endif

	//--------------------------------------------
	// Print out the structure file.
	//--------------------------------------------
	ucell.print_tau();
	//xiaohui modify 2015-03-15, cancel outputfile "STRU_NOW.xyz"
	//ucell.print_cell_xyz("STRU_NOW.xyz");
	//xiaohui add out_stru, 2015-09-30
	if(out_stru==1) ucell.print_cell_cif("STRU_NOW.cif");
	return;
}

void Ions_Move_Basic::check_converged(const double *grad)
{
	TITLE("Ions_Move_Basic","check_converged");
	assert(dim>0);

	//------------------------------------------------
	// check the gradient value
	//------------------------------------------------
	Ions_Move_Basic::largest_grad = 0.0;
	for(int i=0;i<dim;i++)
	{
		if(Ions_Move_Basic::largest_grad < abs(grad[i]))
		{
			Ions_Move_Basic::largest_grad = abs(grad[i]);
		}
	}
	// mohan add 2010-08-06
	Ions_Move_Basic::largest_grad /= ucell.lat0;

	if(test_ion_dynamics)
	{	
		OUT(ofs_running,"old total energy (ry)", etot_p);
		OUT(ofs_running,"new total energy (ry)", etot);
		OUT(ofs_running,"energy difference (ry)", Ions_Move_Basic::ediff);
		OUT(ofs_running,"largest gradient (ry/bohr)",Ions_Move_Basic::largest_grad);
	}

	if(OUT_LEVEL=="ie")
	{
		cout << " ETOT DIFF (eV)       : " << Ions_Move_Basic::ediff*Ry_to_eV << endl;
		cout << " LARGEST GRAD (eV/A)  : " << Ions_Move_Basic::largest_grad * Ry_to_eV / 0.529177 << endl;
	}
	
	const double etot_diff = std::abs(Ions_Move_Basic::ediff);

	// need to update, mohan 2010-07-10
	const double etot_thr = 1.0e-3; // Rydeberg.

	if(Ions_Move_Basic::largest_grad == 0.0)
	{
		ofs_running << " largest force is 0, no movement is possible." << endl;
		ofs_running << " it may converged, otherwise no movement of atom is allowed." << endl;
		Ions_Move_Basic::converged = true;
	}
	// mohan update 2011-04-21
	else if(etot_diff < etot_thr && Ions_Move_Basic::largest_grad < FORCE_THR)
	{
		ofs_running << "\n Ion relaxation is converged!" << endl;
		ofs_running << "\n Energy difference (Ry) = " << etot_diff << endl;
		ofs_running << "\n Largest gradient is (eV/A) = " << largest_grad * Ry_to_eV / 0.529177 << endl;

		Ions_Move_Basic::converged = true;
		++ Ions_Move_Basic::update_iter;
	}
	else
	{
		ofs_running << "\n Ion relaxation is not converged yet (threshold is " 
		<< FORCE_THR * Ry_to_eV / 0.529177 << ")" << endl;
		//cout << "\n etot_diff=" << etot_diff << " etot_thr=" << etot_thr
		//<< " largest_grad=" << largest_grad << " force_thr=" << FORCE_THR << endl;
		Ions_Move_Basic::converged = false;
	}

	return;
}


void Ions_Move_Basic::terminate(void)
{
	TITLE("Ions_Move_Basic","terminate");
	if(Ions_Move_Basic::converged)
	{
		ofs_running << " end of geometry optimization"<<endl;
		OUT(ofs_running,"istep", Ions_Move_Basic::istep);
		OUT(ofs_running,"update iteration", Ions_Move_Basic::update_iter);
		/*
		ofs_running<<"Saving the approximate inverse hessian"<<endl;
		ofstream hess("hess.out");
		for(int i=0;i<dim;i++)
		{
			for(int j=0;j<dim;j++)
			{
				hess << inv_hess(i,j);
			}
		}
		hess.close();
		*/
	}
	else
	{
		ofs_running<<" the maximum number of steps has been reached." << endl;
		ofs_running<<" end of geometry optimization." << endl;
	}

	//-----------------------------------------------------------
	// Print the structure.
	//-----------------------------------------------------------
	ucell.print_tau();
	//xiaohui modify 2015-03-15, cancel outputfile "STRU_NOW.xyz"
	//ucell.print_cell_xyz("STRU_NOW.xyz");
	return;
}

void Ions_Move_Basic::setup_etot(const double &energy_in, const bool judgement)
{
	if(Ions_Move_Basic::istep==1)
	{
		// p == previous
		Ions_Move_Basic::etot_p = energy_in;
		Ions_Move_Basic::etot   = energy_in;
		ediff = etot - etot_p;
	}
	else
	{
		// mohan modify 2010-07-10
		// mohan modify again 2010-07-25
		if(judgement) // for sd
		{
			Ions_Move_Basic::etot = energy_in;
			if(Ions_Move_Basic::etot_p > etot)
			{
				ediff = etot - etot_p;
				Ions_Move_Basic::etot_p = etot;
			}
			else
			{
				// this step will not be accepted
				ediff = 0.0;
			}
		}
		// note: the equlibrium point is not
		// need to be the smallest energy point.	
		else // for bfgs
		{
			Ions_Move_Basic::etot_p = etot;
			Ions_Move_Basic::etot = energy_in;
			ediff = etot - etot_p;
		}
	}

	return;	
}

double Ions_Move_Basic::dot_func(const double* a, const double* b, const int &dim_in)
{
	double result = 0.0;
	for(int i=0; i<dim_in; i++)
	{
		result += a[i] * b[i];
	}
	return result;
}


//----------------------------------------------------------------------------
// second order interpolation scheme, 
//----------------------------------------------------------------------------
void Ions_Move_Basic::second_order(
	const double &e0, // energy of previous step
	const double &e1, // energy at this step
	const double *g0, // gradient at first step
	const double *x, // movement 
	const int &dim,
	double &best_x,
	double &best_e
	)
{
	TITLE("Ions_Move_Basic","second_order");

	// (1) E = ax^2 + bx + c
	// |x> is the movement,or the trust radius 
	// so c=e0,
	// and dE/dx = 2ax + b, so b=|g0>
	// ax^2+<g0|x>+e0=e1
	double bx = dot_func(g0, x, dim); 
	double xx = dot_func(x, x, dim);

	assert(xx!=0);
	double b =  bx / sqrt(xx);

	double a = ( (e1-e0) - bx ) / xx; 
	assert(a!=0);

	// (2) 2ax + b = 0; so best_x=-b/2a,
	best_x = -0.5 * b / a; 
	best_e = a*best_x*best_x+b*best_x+e0;
	ofs_running << " The next E should be ( 2nd order interpolation)" 
	<< best_e * Ry_to_eV << " eV" << endl;

	cout << " The next E should be ( 2nd order interpolation)" 
	<< best_e * Ry_to_eV << " eV" << endl;
	
	return;
}






