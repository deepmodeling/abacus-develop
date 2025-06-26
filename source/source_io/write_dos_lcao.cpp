#include "write_dos_lcao.h"
#include "cal_dos.h"
#include "cal_pdos_gamma.h"
#include "cal_pdos_multik.h"
#include "nscf_fermi_surf.h"
#include "module_parameter/parameter.h"

namespace ModuleIO
{

template <typename T>
void write_dos_lcao(
        const psi::Psi<T>* psi,
		hamilt::Hamilt<T>* p_ham,
        const Parallel_Orbitals &pv, 
        const UnitCell& ucell,
		const K_Vectors& kv,
		const int nbands,
		const elecstate::efermi &energy_fermi,
        const ModuleBase::matrix& ekb,
        const ModuleBase::matrix& wg,
        const double& dos_edelta_ev,
        const double& dos_scale,
        const double& bcoeff,
        const bool out_app_flag,
        const int istep,
        std::ofstream &ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_dos_lcao");
    
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

    // output the DOS file.
    for (int is = 0; is < nspin0; ++is)
    {
        std::stringstream ss;

        ss << PARAM.globalv.global_out_dir << "doss" << is + 1;

		if(istep>=0)
		{
            ss << "g" << istep+1;
		}

        ss << "_nao.txt";

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
		cal_pdos(psi,
				p_ham,
				pv,
				ucell,
				kv,
				nspin0,
				nbands,
				ekb,
				emax,
				emin,
				dos_edelta_ev,
				bcoeff);
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

    ofs_running << " #DOS CALCULATIONS ENDS# " << std::endl;

    return;
}


template void write_dos_lcao(
        const psi::Psi<double>* psi,
		hamilt::Hamilt<double>* p_ham,
        const Parallel_Orbitals &pv, 
        const UnitCell& ucell,
		const K_Vectors& kv,
		const int nbands,
		const elecstate::efermi &energy_fermi,
        const ModuleBase::matrix& ekb,
        const ModuleBase::matrix& wg,
        const double& dos_edelta_ev,
        const double& dos_scale,
        const double& bcoeff,
        const bool out_app_flag,
        const int istep,
        std::ofstream &ofs_running);


template void write_dos_lcao(
        const psi::Psi<std::complex<double>>* psi,
		hamilt::Hamilt<std::complex<double>>* p_ham,
        const Parallel_Orbitals &pv, 
        const UnitCell& ucell,
		const K_Vectors& kv,
		const int nbands,
		const elecstate::efermi &energy_fermi,
        const ModuleBase::matrix& ekb,
        const ModuleBase::matrix& wg,
        const double& dos_edelta_ev,
        const double& dos_scale,
        const double& bcoeff,
        const bool out_app_flag,
        const int istep,
        std::ofstream &ofs_running);

}
