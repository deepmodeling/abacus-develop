#include "write_dos_lcao.h"
#include "cal_dos.h"
#include "cal_pdos_gamma.h"
#include "cal_pdos_multik.h"

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
		const elecstate::Efermi &energy_fermi,
        const ModuleBase::matrix& ekb,
        const ModuleBase::matrix& wg,
        const double& dos_edelta_ev,
        const double& dos_scale,
        const double& bcoeff,
        const bool out_app_flag,
        const int istep,
        const int nspin,
        const int nlocal,
        const bool dos_setemax,
        const double dos_emax_ev,
        const bool dos_setemin,
        const double dos_emin_ev,
        const bool two_fermi,
        const int bndpar,
        const std::string& global_out_dir,
        std::ofstream &ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_dos_lcao");
    
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

    // output the DOS file.
    for (int is = 0; is < nspin0; ++is)
    {
        std::stringstream ss;

        ss << global_out_dir << "doss" << is + 1;

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
				wg,
				istep,
				out_app_flag,
				bndpar);
	}


    // out_dos >= 1: always compute PDOS alongside DOS (LCAO only)
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
            bcoeff,
            istep,
            nlocal,
            nspin,
            global_out_dir);

    ofs_running << " #DOS CALCULATION ENDS# " << std::endl;

    return;
}


template void write_dos_lcao(
        const psi::Psi<double>* psi,
		hamilt::Hamilt<double>* p_ham,
        const Parallel_Orbitals &pv, 
        const UnitCell& ucell,
		const K_Vectors& kv,
		const int nbands,
		const elecstate::Efermi &energy_fermi,
        const ModuleBase::matrix& ekb,
        const ModuleBase::matrix& wg,
        const double& dos_edelta_ev,
        const double& dos_scale,
        const double& bcoeff,
        const bool out_app_flag,
        const int istep,
        const int nspin,
        const int nlocal,
        const bool dos_setemax,
        const double dos_emax_ev,
        const bool dos_setemin,
        const double dos_emin_ev,
        const bool two_fermi,
        const int bndpar,
        const std::string& global_out_dir,
        std::ofstream &ofs_running);


template void write_dos_lcao(
        const psi::Psi<std::complex<double>>* psi,
		hamilt::Hamilt<std::complex<double>>* p_ham,
        const Parallel_Orbitals &pv, 
        const UnitCell& ucell,
		const K_Vectors& kv,
		const int nbands,
		const elecstate::Efermi &energy_fermi,
        const ModuleBase::matrix& ekb,
        const ModuleBase::matrix& wg,
        const double& dos_edelta_ev,
        const double& dos_scale,
        const double& bcoeff,
        const bool out_app_flag,
        const int istep,
        const int nspin,
        const int nlocal,
        const bool dos_setemax,
        const double dos_emax_ev,
        const bool dos_setemin,
        const double dos_emin_ev,
        const bool two_fermi,
        const int bndpar,
        const std::string& global_out_dir,
        std::ofstream &ofs_running);

}
