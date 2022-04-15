#include "charge_extra.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/memory.h"
#include "global.h"
#ifdef __LCAO
#include "../src_lcao/global_fp.h"
#endif

Charge_Extra::Charge_Extra()
{
	init_rho = false;

	// for first-order extrapolation
	this->delta_rho1 = new double*[GlobalV::NSPIN];
	this->delta_rho2 = new double*[GlobalV::NSPIN];
	this->delta_rho = new double*[GlobalV::NSPIN];
	// for second-order extrapolation
	this->delta_rho3 = new double*[GlobalV::NSPIN];

	// PLEASE update the following lines, because
	// the GlobalC::pw.nrxx may not be initialized yet
	// since Charge_Extra is a member of LOOP_ions
	// you can move the initialization of the following 
	// arrays to somewhere else
	// mohan add 2021-03-30
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		delta_rho1[is] = new double[GlobalC::pw.nrxx];
		delta_rho2[is] = new double[GlobalC::pw.nrxx];
		delta_rho[is] = new double[GlobalC::pw.nrxx];

		// for second-order extrapolation
		delta_rho3[is] = new double[GlobalC::pw.nrxx];

		ModuleBase::GlobalFunc::ZEROS(delta_rho1[is], GlobalC::pw.nrxx);
		ModuleBase::GlobalFunc::ZEROS(delta_rho2[is], GlobalC::pw.nrxx);
		ModuleBase::GlobalFunc::ZEROS(delta_rho[is], GlobalC::pw.nrxx);
		ModuleBase::GlobalFunc::ZEROS(delta_rho3[is], GlobalC::pw.nrxx);
	}

	pos_old1 = new double[1];
	pos_old2 = new double[1];
	pos_now = new double[1];
	pos_next = new double[1];
	alpha = 0.0;
	beta = 0.0;
}


Charge_Extra::~Charge_Extra()
{
	if(init_rho)
	{
		for(int i=0; i<dim; i++)
		{
			for(int is=0; is<GlobalV::NSPIN; is++)
			{
				delete[] rho_ion[i][is];
			}
			delete[] rho_ion[i];
		}	
		delete[] rho_ion;
	}

	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		delete[] delta_rho1[is];
		delete[] delta_rho2[is];
		delete[] delta_rho[is];
		delete[] delta_rho3[is];
	}	
	delete[] delta_rho1;
	delete[] delta_rho2;
	delete[] delta_rho;
	delete[] delta_rho3;

	delete[] pos_old1;
	delete[] pos_old2;
	delete[] pos_now;
	delete[] pos_next;
}


void Charge_Extra::allocate_ions(void)
{
	ModuleBase::TITLE("Charge_Extra","allocate_ions");

	// 1: first order extrapolation.
	// 2: second order extrapolation.

	this->dim = 0;
	
	// for the second-order extrapolation	
	pos_dim = GlobalC::ucell.nat * 3;

	delete[] this->pos_old1;
	delete[] this->pos_old2;
	delete[] this->pos_now;
	delete[] this->pos_next;

	this->pos_old1 = new double[pos_dim];
	this->pos_old2 = new double[pos_dim];
	this->pos_now = new double[pos_dim];
	this->pos_next = new double[pos_dim];

	ModuleBase::GlobalFunc::ZEROS(pos_old1, pos_dim);
	ModuleBase::GlobalFunc::ZEROS(pos_old2, pos_dim);
	ModuleBase::GlobalFunc::ZEROS(pos_now, pos_dim);
	ModuleBase::GlobalFunc::ZEROS(pos_next, pos_dim);

	if(init_rho)
	{
		ModuleBase::WARNING_QUIT("Charge_Extra::allocate","rho_ion has been allocated, pls check.");
	}

	this->rho_ion = new double**[dim];

	for(int i=0; i<dim; i++)
	{
		rho_ion[i] = new double*[GlobalV::NSPIN];
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			rho_ion[i][is] = new double[GlobalC::pw.nrxx];
			// first value from charge density.
			for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				rho_ion[i][is][ir] = GlobalC::CHR.rho[is][ir];	
			}
		}
	}	

	init_rho = true;

	ModuleBase::Memory::record("charge_extra","rho_ion",dim*GlobalV::NSPIN*GlobalC::pw.nrxx,"double");

	return;
}


void Charge_Extra::extrapolate_charge()
{
    ModuleBase::TITLE("Charge_Extra","extrapolate_charge");
	//-------------------------------------------------------
    // charge density expolation:
    // pot_order = 0 copy the old potential(nothing is done);
    // pot_order = 3 substrate old atomic charge density and
    // sum the new.
    // If the dynamics is done, this routine extrapolates also
    // the difference between the scf and the atomic one.
    // pot_order = 1 first order extrapolation:
    // rho(t+dt) = 2*rho(t) - rho(t-dt)
    // pot_order = 2 second order extrapolation:
    // rho(t+dt) = rho+
    // alpha*(rho(t) - rho(t-dt))
    // + beta*( rho(t-dt) - rho(t-2*dt) )
	// just use atomic charge.
	//-------------------------------------------------------

	if(GlobalC::pot.chg_extrap == "dm")//xiaohui modify 2015-02-01
	{
		if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			ModuleBase::WARNING_QUIT("Charge_Extra","charge extrapolation method is not available");
		}
		else
		{
			GlobalC::pw.setup_structure_factor();
		}
	}
	// "atomic" extrapolation
	else if(GlobalC::pot.chg_extrap == "atomic")
	{
		double** rho_atom_old = new double*[GlobalV::NSPIN];
		double** rho_atom_new = new double*[GlobalV::NSPIN];

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			rho_atom_old[is] = new double[GlobalC::pw.nrxx];
			rho_atom_new[is] = new double[GlobalC::pw.nrxx];

			ModuleBase::GlobalFunc::ZEROS(rho_atom_old[is], GlobalC::pw.nrxx);
			ModuleBase::GlobalFunc::ZEROS(rho_atom_new[is], GlobalC::pw.nrxx);
		}
		GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_old);
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				delta_rho[is][ir] = GlobalC::CHR.rho[is][ir] - rho_atom_old[is][ir];
			}
		}

		if(GlobalV::OUT_LEVEL != "m") 
		{
			GlobalV::ofs_running << " Setup the structure factor in plane wave basis." << std::endl;
		}
		GlobalC::pw.setup_structure_factor();

		GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_new);

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				GlobalC::CHR.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
			}
		}
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			delete[] rho_atom_old[is];
			delete[] rho_atom_new[is];
		}	
		delete[] rho_atom_old;
		delete[] rho_atom_new;

	}
	// "first-order" extrapolation
	else if(GlobalC::pot.chg_extrap == "first-order")
	{
		double** rho_atom_old = new double*[GlobalV::NSPIN];
		double** rho_atom_new = new double*[GlobalV::NSPIN];

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			rho_atom_old[is] = new double[GlobalC::pw.nrxx];
			rho_atom_new[is] = new double[GlobalC::pw.nrxx];

			ModuleBase::GlobalFunc::ZEROS(rho_atom_old[is], GlobalC::pw.nrxx);
			ModuleBase::GlobalFunc::ZEROS(rho_atom_new[is], GlobalC::pw.nrxx);
		}

		// generate atomic rho
		GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_old);

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				delta_rho2[is][ir] = delta_rho1[is][ir];
				delta_rho1[is][ir] = GlobalC::CHR.rho[is][ir] - rho_atom_old[is][ir];
				delta_rho[is][ir] = 2*delta_rho1[is][ir] - delta_rho2[is][ir];
			}
		}

		if(GlobalV::OUT_LEVEL != "m") 
		{
			GlobalV::ofs_running << " Setup the structure factor in plane wave basis." << std::endl;
		}
		GlobalC::pw.setup_structure_factor();

		GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_new);
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				if(istep == 1)
				{
					GlobalC::CHR.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
				}
				else
				{
					GlobalC::CHR.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
				}
			}
		}
		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			delete[] rho_atom_old[is];
			delete[] rho_atom_new[is];
		}	
		delete[] rho_atom_old;
		delete[] rho_atom_new;
	}

	// "second-order" extrapolation of charge density
	else if(GlobalC::pot.chg_extrap == "second-order")
	{
		double** rho_atom_old = new double*[GlobalV::NSPIN];
		double** rho_atom_new = new double*[GlobalV::NSPIN];

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			rho_atom_old[is] = new double[GlobalC::pw.nrxx];
			rho_atom_new[is] = new double[GlobalC::pw.nrxx];

			ModuleBase::GlobalFunc::ZEROS(rho_atom_old[is], GlobalC::pw.nrxx);
			ModuleBase::GlobalFunc::ZEROS(rho_atom_new[is], GlobalC::pw.nrxx);
		}

		// generate atomic_rho
		GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_old);

		// compute alpha and beta
		find_alpha_and_beta();

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				delta_rho3[is][ir] = delta_rho2[is][ir];
				delta_rho2[is][ir] = delta_rho1[is][ir];
				delta_rho1[is][ir] = GlobalC::CHR.rho[is][ir] - rho_atom_old[is][ir];
				delta_rho[is][ir] = delta_rho1[is][ir] + 
					alpha * (delta_rho1[is][ir] - delta_rho2[is][ir]) +
					beta * (delta_rho2[is][ir] - delta_rho3[is][ir]);
			}
		}

		//xiaohui add 'GlobalV::OUT_LEVEL', 2015-09-16
		if(GlobalV::OUT_LEVEL != "m") 
		{
			GlobalV::ofs_running << " Setup the structure factor in plane wave basis." << std::endl;
		}

		// setup the structure factor
		GlobalC::pw.setup_structure_factor();

		// generate atomic rho
		GlobalC::CHR.atomic_rho(GlobalV::NSPIN,rho_atom_new);

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
			{
				if(istep == 1)
				{
					GlobalC::CHR.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
				}
				else if(istep == 2)
				{
					delta_rho[is][ir] = 2*delta_rho1[is][ir] - delta_rho2[is][ir];
					GlobalC::CHR.rho[is][ir] = delta_rho1[is][ir] + rho_atom_new[is][ir];
				}
				else
				{
					GlobalC::CHR.rho[is][ir] = delta_rho[is][ir] + rho_atom_new[is][ir];
				}
			}
		}

		for(int is=0; is<GlobalV::NSPIN; is++)
		{
			delete[] rho_atom_old[is];
			delete[] rho_atom_new[is];
		}	

		delete[] rho_atom_old;
		delete[] rho_atom_new;
	}
	else
	{
		ModuleBase::WARNING_QUIT("potential::init_pot","chg_extrap parameter is wrong!");
	}

    return;
}


void Charge_Extra::find_alpha_and_beta(void)
{
	double a11,a12,a21,a22;
	double b1,b2;
	double detA;

	a11 = 0.0;
	a12 = 0.0;
	a21 = 0.0;
	a22 = 0.0;
	b1 = 0.0;
	b2 = 0.0;
	detA = 0.0;

	if(istep >= 3)
	{
		int iat=0;
		for(int it = 0;it < GlobalC::ucell.ntype;it++)
		{
			//Atom* atom = &GlobalC::ucell.atoms[it];
			for(int ia =0;ia< GlobalC::ucell.atoms[it].na;ia++)
			{
				a11 += (pos_now[3*iat  ] - pos_old1[3*iat  ]) * (pos_now[3*iat  ] - pos_old1[3*iat  ]) + 
					(pos_now[3*iat+1] - pos_old1[3*iat+1]) * (pos_now[3*iat+1] - pos_old1[3*iat+1]) + 
					(pos_now[3*iat+2] - pos_old1[3*iat+2]) * (pos_now[3*iat+2] - pos_old1[3*iat+2]);

				a12 += (pos_now[3*iat  ] - pos_old1[3*iat  ]) * (pos_old1[3*iat  ] - pos_old2[3*iat  ]) + 
					(pos_now[3*iat+1] - pos_old1[3*iat+1]) * (pos_old1[3*iat+1] - pos_old2[3*iat+1]) + 
					(pos_now[3*iat+2] - pos_old1[3*iat+2]) * (pos_old1[3*iat+2] - pos_old2[3*iat+2]);

				a22 += (pos_old1[3*iat  ] - pos_old2[3*iat  ]) * (pos_old1[3*iat  ] - pos_old2[3*iat  ]) + 
					(pos_old1[3*iat+1] - pos_old2[3*iat+1]) * (pos_old1[3*iat+1] - pos_old2[3*iat+1]) + 
					(pos_old1[3*iat+2] - pos_old2[3*iat+2]) * (pos_old1[3*iat+2] - pos_old2[3*iat+2]);

				b1 += -((pos_now[3*iat  ] - pos_next[3*iat  ]) * (pos_now[3*iat  ] - pos_old1[3*iat  ]) +
					(pos_now[3*iat+1] - pos_next[3*iat+1]) * (pos_now[3*iat+1] - pos_old1[3*iat+1]) +
					(pos_now[3*iat+2] - pos_next[3*iat+2]) * (pos_now[3*iat+2] - pos_old1[3*iat+2]));

				b2 += -((pos_now[3*iat  ] - pos_next[3*iat  ]) * (pos_old1[3*iat  ] - pos_old2[3*iat  ]) +
					(pos_now[3*iat+1] - pos_next[3*iat+1]) * (pos_old1[3*iat+1] - pos_old2[3*iat+1]) +
					(pos_now[3*iat+2] - pos_next[3*iat+2]) * (pos_old1[3*iat+2] - pos_old2[3*iat+2]));

				iat++;
			}
		}
	}
	a21 = a12;
	detA = a11 * a22 - a21 * a12;

	if(detA > 1.0E-12)
	{
		alpha = (b1 * a22 - b2 * a12) / detA;
		beta = (b2 * a11 - b1 * a21) / detA;
		if(abs(alpha) >10) alpha=1.0;
		if(abs(beta)>10) beta=0.0;
	}
	else
	{
		alpha = 0.0;
		beta = 0.0;
		if(a11 != 0)
		{
			alpha = b1 / a11 ;
		}
		if(abs(alpha) >10) alpha=1.0;
	}
	return;
}

void Charge_Extra::save_pos_next(const UnitCell_pseudo& ucell)
{
	ucell.save_cartesian_position(this->pos_next);
	return;
}

void Charge_Extra::update_istep(const int &step)
{
	this->istep = step;
	return;
}

void Charge_Extra::update_all_pos(const UnitCell_pseudo& ucell)
{
	const int total_freedom = ucell.nat * 3;
	for(int i=0;i<total_freedom;i++)
	{
		this->pos_old2[i] = this->pos_old1[i];
		this->pos_old1[i] = this->pos_now[i];
	}
	ucell.save_cartesian_position(this->pos_now);
	return;
}
