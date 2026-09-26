#include "write_dos_pw.h"
#include "cal_dos.h"
#include "source_base/parallel_reduce.h"

void ModuleIO::write_dos_pw(
		const UnitCell& ucell,
		const ModuleBase::matrix& ekb,
		const ModuleBase::matrix& wg,
		const K_Vectors& kv,
		const int nbands,
		const int istep_in,
		const elecstate::Efermi &energy_fermi,
		const double& dos_edelta_ev,
		const double& dos_scale,
		const double& bcoeff,
		const int nspin,
		const int out_dos,
		const bool dos_setemax,
		const double dos_emax_ev,
		const bool dos_setemin,
		const double dos_emin_ev,
		const bool two_fermi,
		const bool out_app_flag,
		const int bndpar,
		const std::string& global_out_dir,
		std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_dos_pw");

    const int nspin0 = (nspin == 2) ? 2 : 1;

    double emax = 0.0;
    double emin = 0.0;

	prepare_dos(ofs_running, 
			energy_fermi, 
			ekb,
			kv.get_nks(),
			nbands,
			dos_edelta_ev,
            dos_scale,
			emax, 
			emin,
			dos_setemax,
			dos_emax_ev,
			dos_setemin,
			dos_emin_ev,
			two_fermi);

    for (int is = 0; is < nspin0; ++is)
    {
        // DOS_ispin contains not smoothed dos
		std::stringstream ss;
		ss << global_out_dir << "dos";

		if(nspin0==2)
		{
			ss << "s" << is + 1;
		}
		else
		{
			// do nothing;
		}

		ss << ".txt";

        ModuleBase::GlobalFunc::OUT(ofs_running, "DOS file", ss.str());

		ModuleIO::cal_dos(is,
				ss.str(),
				dos_edelta_ev,
				emax,
				emin,
				bcoeff,
				kv.get_nks(),
				kv.get_nkstot(),
				kv.wk,
				kv.isk,
				nbands,
				ekb,
				wg,
				istep_in,
				out_app_flag,
				bndpar);
	}


    if (out_dos == 2)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dos_pw","PW basis do not support PDOS calculations yet.");
    }

    ofs_running << " #DOS CALCULATION ENDS# " << std::endl;
}
