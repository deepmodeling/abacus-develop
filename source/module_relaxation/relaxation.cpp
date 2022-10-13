#include "ions.h"
#include "../src_pw/global.h" // use chr.
#include "../src_io/print_info.h"
#include "variable_cell.h" // mohan add 2021-02-01

// The interface for relaxation
bool Ions::relaxation(ModuleBase::matrix force, ModuleBase::matrix stress, const int &istep, int &force_step, int &stress_step)
{
	ModuleBase::TITLE("Ions","after_scf");

	// should not do it this way, will change after the refactor of ucell class
	GlobalC::ucell.ionic_position_updated = false;
	GlobalC::ucell.cell_parameter_updated = false;

	//stop in last step
	if(istep==GlobalV::RELAX_NMAX)
	{
		return 1;
	}
	//choose what to do next
	if(GlobalV::CALCULATION!="cell-relax") force_step = istep;
	if(this->if_do_relax()) 
	{
		//do relax calculation and generate next structure 
		bool converged = 0;
		converged = this->do_relax(istep, force_step, force, GlobalC::en.etot);
		if(!converged) 
		{
			this->reset_after_relax(istep);
			GlobalC::ucell.ionic_position_updated = true;
			return converged;
		}
		else if(GlobalV::CALCULATION!="cell-relax")
		{
			return converged;
		}
	}
	if(this->if_do_cellrelax())
	{
		//do cell relax calculation and generate next structure
		bool converged = 0;
		converged = this->do_cellrelax(istep,stress_step, stress, GlobalC::en.etot);
		if(!converged)
		{
			GlobalC::ucell.cell_parameter_updated = true;
			this->reset_after_cellrelax(force_step, stress_step);
		}
		return converged;
	}

    return 1;
}

bool Ions::if_do_relax()
{
	ModuleBase::TITLE("Ions","if_do_relax");
	if(GlobalV::CALCULATION=="relax"||GlobalV::CALCULATION=="cell-relax")
	{
		if(!GlobalC::ucell.if_atoms_can_move()) 
		{
			ModuleBase::WARNING("Ions","No atom is allowed to move!");
			return 0;
		}
//		if(!IMM.get_converged()) return 1;
		else 
		{
			assert(GlobalV::CAL_FORCE==1);
			return 1;
		}
	}
	else return 0;
}
bool Ions::if_do_cellrelax()
{
	ModuleBase::TITLE("Ions","if_do_cellrelax");
	if(GlobalV::CALCULATION=="cell-relax")
	{
		if(!GlobalC::ucell.if_cell_can_change()) 
		{
			ModuleBase::WARNING("Ions", "Lattice vectors are not allowed to change!");
			return 0;
		}
		else if(GlobalC::ucell.if_atoms_can_move()&&!IMM.get_converged())
		{
			GlobalV::ofs_running<<"Note: Need to wait for atomic relaxation first!";
			return 0;
		}
		else 
		{
			assert(GlobalV::CAL_STRESS==1);
			return 1;
		}
	}
	else return 0;
}
bool Ions::do_relax(const int& istep, int& jstep, const ModuleBase::matrix& ionic_force, const double& total_energy)
{
	ModuleBase::TITLE("Ions","do_relax");
	IMM.cal_movement(istep, jstep, ionic_force, total_energy);
	++jstep;
	return IMM.get_converged();
}
bool Ions::do_cellrelax(const int&istep, const int& stress_step, const ModuleBase::matrix& stress, const double& total_energy)
{
	ModuleBase::TITLE("Ions","do_cellrelax");
	LCM.cal_lattice_change(istep, stress_step, stress, total_energy);
    return LCM.get_converged();
}
void Ions::reset_after_relax(const int& istep)
{
	ModuleBase::TITLE("Ions","reset_after_relax");
	GlobalV::ofs_running << " Setup the structure factor in plane wave basis." << std::endl;
	GlobalC::sf.setup_structure_factor(&GlobalC::ucell,GlobalC::rhopw);
}

void Ions::reset_after_cellrelax(int& f_step, int& s_step)
{
	ModuleBase::TITLE("Ions","reset_after_cellrelax");
	Variable_Cell::init_after_vc();
	GlobalC::pot.init_pot(s_step, GlobalC::sf.strucFac); //LiuXh add 20180619

	f_step = 1;
	++s_step;
}
