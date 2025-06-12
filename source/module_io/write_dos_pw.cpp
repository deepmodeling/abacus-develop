#include "write_dos_pw.h"
#include "cal_dos.h"
#include "nscf_fermi_surf.h"
#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"

void ModuleIO::write_dos_pw(
		const UnitCell& ucell,
		const ModuleBase::matrix& ekb,
		const ModuleBase::matrix& wg,
		const K_Vectors& kv,
		const int nbands,
		const elecstate::efermi &energy_fermi,
		const double& dos_edelta_ev,
		const double& dos_scale,
		const double& bcoeff,
		std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_dos_pw");

    const int nspin0 = (PARAM.inp.nspin == 2) ? 2 : 1;

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
			emin);

    for (int is = 0; is < nspin0; ++is)
    {
        // DOS_ispin contains not smoothed dos
        std::stringstream ss;
        ss << PARAM.globalv.global_out_dir << "doss" << is + 1 << "_pw.txt";

        std::stringstream ss1;
        ss1 << PARAM.globalv.global_out_dir << "doss" << is + 1 << "s_pw.txt";

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
				wg);
	}


    if (PARAM.inp.out_dos == 2)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dos_pw","PW basis do not support PDOS calculations yet.");
    }

    if(PARAM.inp.out_dos == 3)
    {
        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss3;
            ss3 << PARAM.globalv.global_out_dir << "fermi" << is << ".bxsf";
            nscf_fermi_surface(ss3.str(), nbands, energy_fermi.ef, kv, ucell, ekb);
        }
    }

    ofs_running << " DOS CALCULATIONS ENDS." << std::endl; 
}
