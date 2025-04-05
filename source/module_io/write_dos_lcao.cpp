#include "write_dos_lcao.h"
#include "cal_dos.h"

#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "write_orb_info.h"

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"

#include "module_base/parallel_reduce.h"
#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"


// for gamma only
template <>
void ModuleIO::write_dos_lcao(const UnitCell& ucell,
                              const psi::Psi<double>* psi,
                              const Parallel_Orbitals& pv,
                              const ModuleBase::matrix& ekb,
                              const ModuleBase::matrix& wg,
                              const double& dos_edelta_ev,
                              const double& dos_scale,
                              const double& bcoeff,
                              const K_Vectors& kv,
                              const int nbands,
                              const elecstate::efermi &energy_fermi,
                              hamilt::Hamilt<double>* p_ham,
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
        ss << PARAM.globalv.global_out_dir << "DOS" << is + 1 << ".dat";
        std::stringstream ss1;
        ss1 << PARAM.globalv.global_out_dir << "DOS" << is + 1 << "_smear.dat";

		ModuleIO::cal_dos(is,
				ss.str(),
				ss1.str(),
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
        cal_pdos_gamma(nspin0,
                emax,
                emin,
                dos_edelta_ev,
                bcoeff,
                p_ham,
                psi,
                pv);
    }

    ofs_running << " DOS CALCULATIONS ENDS." << std::endl;

    return;
}

// for multi-k case
template <>
void ModuleIO::write_dos_lcao(const UnitCell& ucell,
                              const psi::Psi<std::complex<double>>* psi,
                              const Parallel_Orbitals& pv,
                              const ModuleBase::matrix& ekb,
                              const ModuleBase::matrix& wg,
                              const double& dos_edelta_ev,
                              const double& dos_scale,
                              const double& bcoeff,
                              const K_Vectors& kv,
                              const int nbands,
                              const elecstate::efermi &energy_fermi,
                              hamilt::Hamilt<std::complex<double>>* p_ham,
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

    if (PARAM.inp.out_dos == 2)
    {
        cal_pdos_gamma(nspin0,
                emax,
                emin,
                dos_edelta_ev,
                bcoeff,
                p_ham,
                psi,
                pv);
    }

    // output the DOS file.
    for (int is = 0; is < nspin0; ++is)
    {
        std::stringstream ss;
        ss << PARAM.globalv.global_out_dir << "DOS" << is + 1 << ".dat";
        std::stringstream ss1;
        ss1 << PARAM.globalv.global_out_dir << "DOS" << is + 1 << "_smear.dat";

		ModuleIO::cal_dos(is,
				ss.str(),
				ss1.str(),
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

    ofs_running << " DOS CALCULATIONS ENDS." << std::endl; 

    return;
}
