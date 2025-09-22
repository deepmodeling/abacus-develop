#include "source_io/ctrl_output_fp.h" // use ctrl_output_fp() 

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_output_fp(UnitCell& ucell, 
		elecstate::ElecStateLCAO<TK>* pelec, 
		const int istep)
{
/*
    ModuleBase::TITLE("ModuleIO", "ctrl_output_fp");
    ModuleBase::timer::tick("ModuleIO", "ctrl_output_fp");

    const bool out_app_flag = PARAM.inp.out_app_flag;
    const bool gamma_only = PARAM.globalv.gamma_only_local;
    const int nspin = PARAM.inp.nspin;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;

	// 1) write charge density
	if (PARAM.inp.out_chg[0] > 0)
	{
		for (int is = 0; is < PARAM.inp.nspin; is++)
		{
			this->pw_rhod->real2recip(this->chr.rho_save[is], this->chr.rhog_save[is]);
			std::string fn =PARAM.globalv.global_out_dir + "/chgs" + std::to_string(is + 1) + ".cube";
			ModuleIO::write_vdata_palgrid(Pgrid,
					this->chr.rho_save[is],
					is,
					PARAM.inp.nspin,
					istep,
					fn,
					this->pelec->eferm.get_efval(is),
					&(ucell),
					PARAM.inp.out_chg[1],
					1);

			if (XC_Functional::get_ked_flag())
			{
				fn =PARAM.globalv.global_out_dir + "/taus" + std::to_string(is + 1) + ".cube";
				ModuleIO::write_vdata_palgrid(Pgrid,
						this->chr.kin_r_save[is],
						is,
						PARAM.inp.nspin,
						istep,
						fn,
						this->pelec->eferm.get_efval(is),
						&(ucell));
			}
		}
	}


	// 2) write potential
	if (PARAM.inp.out_pot == 1 || PARAM.inp.out_pot == 3)
	{
		for (int is = 0; is < PARAM.inp.nspin; is++)
		{
			std::string fn =PARAM.globalv.global_out_dir + "/pots" + std::to_string(is + 1) + ".cube";

			ModuleIO::write_vdata_palgrid(Pgrid,
					this->pelec->pot->get_effective_v(is),
					is,
					PARAM.inp.nspin,
					istep,
					fn,
					0.0, // efermi
					&(ucell),
					3,  // precision
					0); // out_fermi
		}
	}
	else if (PARAM.inp.out_pot == 2)
	{
		std::string fn =PARAM.globalv.global_out_dir + "/pot_es.cube";
		ModuleIO::write_elecstat_pot(
#ifdef __MPI
				this->pw_big->bz,
				this->pw_big->nbz,
#endif
				fn,
				istep,
				this->pw_rhod,
				&this->chr,
				&(ucell),
				this->pelec->pot->get_fixed_v(),
				this->solvent);
	}


	// 3) write ELF
	if (PARAM.inp.out_elf[0] > 0)
	{
		this->chr.cal_elf = true;
		Symmetry_rho srho;
		for (int is = 0; is < PARAM.inp.nspin; is++)
		{
			srho.begin(is, this->chr, this->pw_rhod, ucell.symm);
		}

		std::string out_dir =PARAM.globalv.global_out_dir;
		ModuleIO::write_elf(
#ifdef __MPI
				this->pw_big->bz,
				this->pw_big->nbz,
#endif
				out_dir,
				istep,
				PARAM.inp.nspin,
				this->chr.rho,
				this->chr.kin_r,
				this->pw_rhod,
				this->Pgrid,
				&(ucell),
				PARAM.inp.out_elf[1]);
	}

    ModuleBase::timer::tick("ModuleIO", "ctrl_output_fp");

*/

}

} // End ModuleIO


// For gamma only
template void ModuleIO::ctrl_output_lcao<double, double>(UnitCell& ucell, 
		const int istep);

// For multiple k-points
template void ModuleIO::ctrl_output_lcao<std::complex<double>, double>(UnitCell& ucell, 
		const int istep);

template void ModuleIO::ctrl_output_lcao<std::complex<double>, std::complex<double>>(UnitCell& ucell, 
		const int istep);

