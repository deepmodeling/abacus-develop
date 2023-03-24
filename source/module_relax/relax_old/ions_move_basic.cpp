#include "ions_move_basic.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

int Ions_Move_Basic::dim=0;
bool Ions_Move_Basic::converged=false;
double Ions_Move_Basic::largest_grad=0.0;
int Ions_Move_Basic::update_iter=0;
int Ions_Move_Basic::istep=0;

double Ions_Move_Basic::ediff=0.0;
double Ions_Move_Basic::etot=0.0;
double Ions_Move_Basic::etot_p=0.0;

double Ions_Move_Basic::trust_radius=0.0;
double Ions_Move_Basic::trust_radius_old=0.0;
double Ions_Move_Basic::relax_bfgs_rmax = -1.0; // default is 0.8
double Ions_Move_Basic::relax_bfgs_rmin = -1.0; // default is 1e-5
double Ions_Move_Basic::relax_bfgs_init = -1.0; // default is 0.5
double Ions_Move_Basic::best_xxx = 1.0;

int Ions_Move_Basic::out_stru=0;

void Ions_Move_Basic::setup_gradient(double *grad, const ModuleBase::matrix &force)
{
	ModuleBase::TITLE("Ions_Move_Basic","setup_gradient");
	
	assert(GlobalC::ucell.ntype>0);
	assert(grad!=NULL);
	assert(dim == 3*GlobalC::ucell.nat);

	ModuleBase::GlobalFunc::ZEROS(grad, dim);

	// (1) init gradient
	// the unit of force: Ry/Bohr.
	// the unit of gradient: 
	int iat=0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		Atom* atom = &GlobalC::ucell.atoms[it];
		for(int ia =0;ia< GlobalC::ucell.atoms[it].na;ia++)
		{	
			if(atom->mbl[ia].x == 1)
			{
				grad[3*iat  ] = -force(iat, 0)*GlobalC::ucell.lat0;
				//this->grad[3*iat  ] = -force(iat, 0);
			}
			if(atom->mbl[ia].y == 1)
			{
				grad[3*iat+1] = -force(iat, 1)*GlobalC::ucell.lat0;
			}
			if(atom->mbl[ia].z == 1)
			{
				grad[3*iat+2] = -force(iat, 2)*GlobalC::ucell.lat0;
				//std::cout << " grad=" << grad[3*iat+2] << std::endl;
			}
			++iat;
		}
	}

	return;
}

void Ions_Move_Basic::move_atoms(double *move)
{
	ModuleBase::TITLE("Ions_Move_Basic","move_atoms");

	assert(move!=NULL);

	//------------------------
	// for test only
	//------------------------
	if(GlobalV::test_relax_method)
	{
		int iat=0;
		GlobalV::ofs_running << "\n movement of ions (unit is Bohr) : " << std::endl;
		GlobalV::ofs_running << " " << std::setw(12) << "Atom" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::endl;
		for(int it = 0;it < GlobalC::ucell.ntype;it++)
		{
			for(int ia =0;ia< GlobalC::ucell.atoms[it].na;ia++)
			{
				std::stringstream ss;
				ss << "move_" << GlobalC::ucell.atoms[it].label << ia+1;
				GlobalV::ofs_running << " " 
					<< std::setw(12) << ss.str().c_str()
					<< std::setw(15) << move[3*iat+0] 
					<< std::setw(15) << move[3*iat+1] 
					<< std::setw(15) << move[3*iat+2] << std::endl;
				iat++;
			}
		}
		assert( iat == GlobalC::ucell.nat );
	}

    int iat = 0;
    const double move_threshold = 1.0e-10;
    ModuleBase::Vector3<double> *move_ion = new ModuleBase::Vector3<double> [GlobalC::ucell.nat];
    for(int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        Atom* atom = &GlobalC::ucell.atoms[it];
        for(int ia = 0; ia < atom->na; ia++)
        {
            for (int k = 0; k < 3; ++k)
            {
                if( abs(move[3*iat + k]) > move_threshold && atom->mbl[ia][k])
                {
                    move_ion[iat][k] = move[3*iat + k] / GlobalC::ucell.lat0;
                }
                else
                {
                    move_ion[iat][k] = 0;
                }
            }
            move_ion[iat] = move_ion[iat] * GlobalC::ucell.GT;
            iat++;
        }
    }
    assert( iat == GlobalC::ucell.nat );
    GlobalC::ucell.update_pos_taud(move_ion);

	//--------------------------------------------
	// Print out the structure file.
	//--------------------------------------------
	GlobalC::ucell.print_tau();
	//xiaohui modify 2015-03-15, cancel outputfile "STRU_NOW.xyz"
	//GlobalC::ucell.print_cell_xyz("STRU_NOW.xyz");
	//xiaohui add out_stru, 2015-09-30
	if(out_stru==1)
	{
		GlobalC::ucell.print_cell_cif("STRU_NOW.cif");
	}
    delete[] move_ion;
	return;
}

void Ions_Move_Basic::check_converged(const double *grad)
{
	ModuleBase::TITLE("Ions_Move_Basic","check_converged");
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
	Ions_Move_Basic::largest_grad /= GlobalC::ucell.lat0;

	if(GlobalV::test_relax_method)
	{	
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"old total energy (ry)", etot_p);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"new total energy (ry)", etot);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"energy difference (ry)", Ions_Move_Basic::ediff);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"largest gradient (ry/bohr)",Ions_Move_Basic::largest_grad);
	}

	if(GlobalV::OUT_LEVEL=="ie")
	{
		std::cout << " ETOT DIFF (eV)       : " << Ions_Move_Basic::ediff*ModuleBase::Ry_to_eV << std::endl;
		std::cout << " LARGEST GRAD (eV/A)  : " << Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / 0.529177 << std::endl;
	}
	
	const double etot_diff = std::abs(Ions_Move_Basic::ediff);

	// need to update, mohan 2010-07-10
	const double etot_thr = 1.0e-3; // Rydeberg.

	if(Ions_Move_Basic::largest_grad == 0.0)
	{
		GlobalV::ofs_running << " largest force is 0, no movement is possible." << std::endl;
		GlobalV::ofs_running << " it may converged, otherwise no movement of atom is allowed." << std::endl;
		Ions_Move_Basic::converged = true;
	}
	// mohan update 2011-04-21
	else if(etot_diff < etot_thr && Ions_Move_Basic::largest_grad < GlobalV::FORCE_THR)
	{
		GlobalV::ofs_running << "\n Ion relaxation is converged!" << std::endl;
		GlobalV::ofs_running << "\n Energy difference (Ry) = " << etot_diff << std::endl;
		GlobalV::ofs_running << "\n Largest gradient is (eV/A) = " << largest_grad * ModuleBase::Ry_to_eV / 0.529177 << std::endl;

		Ions_Move_Basic::converged = true;
		++ Ions_Move_Basic::update_iter;
	}
	else
	{
		GlobalV::ofs_running << "\n Ion relaxation is not converged yet (threshold is " 
		<< GlobalV::FORCE_THR * ModuleBase::Ry_to_eV / 0.529177 << ")" << std::endl;
		//std::cout << "\n etot_diff=" << etot_diff << " etot_thr=" << etot_thr
		//<< " largest_grad=" << largest_grad << " force_thr=" << GlobalV::FORCE_THR << std::endl;
		Ions_Move_Basic::converged = false;
	}

	return;
}


void Ions_Move_Basic::terminate(void)
{
	ModuleBase::TITLE("Ions_Move_Basic","terminate");
	if(Ions_Move_Basic::converged)
	{
		GlobalV::ofs_running << " end of geometry optimization"<<std::endl;
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"istep", Ions_Move_Basic::istep);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"update iteration", Ions_Move_Basic::update_iter);
		/*
		GlobalV::ofs_running<<"Saving the approximate inverse hessian"<<std::endl;
		std::ofstream hess("hess.out");
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
		GlobalV::ofs_running<<" the maximum number of steps has been reached." << std::endl;
		GlobalV::ofs_running<<" end of geometry optimization." << std::endl;
	}

	//-----------------------------------------------------------
	// Print the structure.
	//-----------------------------------------------------------
	GlobalC::ucell.print_tau();
	//xiaohui modify 2015-03-15, cancel outputfile "STRU_NOW.xyz"
	//GlobalC::ucell.print_cell_xyz("STRU_NOW.xyz");
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





