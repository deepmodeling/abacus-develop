#include "source_pw/module_ofdft/of_print_info.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"
#include "source_base/formatter.h"

/**
 * @brief Print nessecary information to the screen,
 * and write the components of the total energy into running_log.
 */
void OFDFT::print_info(const int iter,
	ModuleBase::TimePoint &iter_time,
	const double &energy_current,
	const double &energy_last,
	const double &normdLdphi,
	const elecstate::ElecState *pelec,
	KEDF_Manager *kedf_manager,
	const bool conv_esolver)
{
    if (iter == 0)
    {
        std::cout << " ============================= Running OFDFT "
                     "=============================="
                  << std::endl;
        std::cout << " ITER       ETOT/eV           EDIFF/eV        EFERMI/eV    POTNORM   TIME/s"
                  << std::endl;
    }

    std::map<std::string, std::string> prefix_map = {
        {"cg1", "CG"},
        {"cg2", "CG"},
        {"tn", "TN"}
    };
    std::string iteration = prefix_map[PARAM.inp.of_method] + std::to_string(iter);
    double duration = ModuleBase::get_duration(iter_time, ModuleBase::get_time());
    std::cout << " " << std::setw(8) << iteration
              << std::setw(18) << std::scientific << std::setprecision(8) << energy_current * ModuleBase::Ry_to_eV
              << std::setw(18) << (energy_current - energy_last) * ModuleBase::Ry_to_eV
              << std::setw(13) << std::setprecision(4) << pelec->eferm.get_efval(0) * ModuleBase::Ry_to_eV
              << std::setw(13) << std::setprecision(4) << normdLdphi
              << std::setw(6) << std::fixed << std::setprecision(2) << duration << std::endl;

    GlobalV::ofs_running << std::setprecision(12);
    GlobalV::ofs_running << std::setiosflags(std::ios::right);

    GlobalV::ofs_running << "\nIter" << iter << ": the norm of potential is " << normdLdphi << std::endl;

    std::vector<std::string> titles;
    std::vector<double> energies_Ry;
    std::vector<double> energies_eV;
	if ((PARAM.inp.out_band[0] > 0 && 
				((iter + 1) % PARAM.inp.out_band[0] == 0 || 
				 conv_esolver || 
				 iter == PARAM.inp.scf_nmax)) || 
			PARAM.inp.init_chg == "file")
    {
        titles.push_back("E_Total");
        energies_Ry.push_back(pelec->f_en.etot);
        titles.push_back("E_Kinetic");
        energies_Ry.push_back(pelec->f_en.ekinetic);
        titles.push_back("E_Hartree");
        energies_Ry.push_back(pelec->f_en.hartree_energy);
        titles.push_back("E_xc");
        energies_Ry.push_back(pelec->f_en.etxc - pelec->f_en.etxcc);
        titles.push_back("E_LocalPP");
        energies_Ry.push_back(pelec->f_en.e_local_pp);
        titles.push_back("E_Ewald");
        energies_Ry.push_back(pelec->f_en.ewald_energy);

        kedf_manager->record_energy(titles, energies_Ry);
        
        std::string vdw_method = PARAM.inp.vdw_method;
        if (vdw_method == "d2") // Peize Lin add 2014-04, update 2021-03-09
        {
            titles.push_back("E_vdwD2");
            energies_Ry.push_back(pelec->f_en.evdw);
        }
        else if (vdw_method == "d3_0" || vdw_method == "d3_bj") // jiyy add 2019-05, update 2021-05-02
        {
            titles.push_back("E_vdwD3");
            energies_Ry.push_back(pelec->f_en.evdw);
        }
        if (PARAM.inp.imp_sol)
        {
            titles.push_back("E_sol_el");
            energies_Ry.push_back(pelec->f_en.esol_el);
            titles.push_back("E_sol_cav");
            energies_Ry.push_back(pelec->f_en.esol_cav);
        }
        if (PARAM.inp.efield_flag)
        {
            titles.push_back("E_efield");
            energies_Ry.push_back(elecstate::Efield::etotefield);
        }
        if (PARAM.inp.gate_flag)
        {
            titles.push_back("E_gatefield");
            energies_Ry.push_back(elecstate::Gatefield::etotgatefield);
        }
        if (PARAM.inp.ml_exx)
        {
            titles.push_back("E_ML-EXX");
            energies_Ry.push_back(pelec->f_en.ml_exx);
        }
    }
    else
    {
        titles.push_back("E_Total");
        energies_Ry.push_back(pelec->f_en.etot);
    }

    if (PARAM.globalv.two_fermi)
    {
        titles.push_back("E_Fermi_up");
        energies_Ry.push_back(pelec->eferm.get_efval(0));
        titles.push_back("E_Fermi_dw");
        energies_Ry.push_back(pelec->eferm.get_efval(1));
    }
    else
    {
        titles.push_back("E_Fermi");
        energies_Ry.push_back(pelec->eferm.get_efval(0));
    }
    energies_eV.resize(energies_Ry.size());
    std::transform(energies_Ry.begin(), energies_Ry.end(), energies_eV.begin(), [](double energy) {
        return energy * ModuleBase::Ry_to_eV;
    });
    FmtTable table(/*titles=*/{"Energy", "Rydberg", "eV"}, 
                   /*nrows=*/titles.size(), 
                   /*formats=*/{"%20s", "%20.12f", "%20.12f"}, 0);
    table << titles << energies_Ry << energies_eV;
    GlobalV::ofs_running << table.str() << std::endl;

    // reset the iter_time for the next iteration
    iter_time = ModuleBase::get_time();
}
