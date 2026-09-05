#include "relax_driver.h"
#include "source_base/formatter.h"
#include "source_base/global_file.h"
#include "source_base/version.h"
#include "source_cell/cif_io.h"
#include "source_io/module_json/output_info.h"
#include "source_cell/output_log.h"
#include "source_io/module_output/print_info.h"
#include "source_base/module_out/read_exit_file.h"
#include "source_io/module_parameter/parameter.h"
#include "source_cell/print_cell.h"

#include <ctime>

void Relax_Driver::relax_driver(
        ModuleESolver::ESolver* p_esolver,
        UnitCell& ucell,
        const Input_para& inp,
        std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Relax_Driver", "relax_driver");
    ModuleBase::timer::start("Relax_Driver", "relax_driver");

    this->init_relax(ucell.nat, inp);

    // steps[0]: istep (main iteration step)
    // steps[1]: force_step
    // steps[2]: stress_step
    std::vector<int> steps = {0, 1, 1};

    // Main iteration loop for relaxation calculations
    // For scf/nscf calculations, relax_step returns true immediately,
    // so the loop exits after one iteration
    double etot = 0.0;
    ModuleBase::matrix stress(3, 3);

    while (steps[0] < inp.relax_nmax)
    {
        ModuleBase::matrix force(ucell.nat, 3);

        this->iter_info(steps, inp);
        this->esolve(steps[0], p_esolver, ucell, inp, force, stress, etot);
        this->stru_out(steps[0], ucell, inp, etot, stress);
        bool converged = this->relax_step(steps, p_esolver, ucell, inp, force, stress, etot, ofs_running);
        this->json_out(p_esolver, ucell, inp, force, stress);

        // Check stop conditions
        if (converged)
        {
            // Relaxation converged, exit loop immediately
            break;
        }
        else if (ModuleIO::read_exit_file(GlobalV::MY_RANK, "EXIT", ofs_running))
        {
            // EXIT file detected, exit loop
            break;
        }

        ++steps[0];
    }

    this->final_out(steps[0], ucell, inp, etot, stress);

    ModuleBase::timer::end("Relax_Driver", "relax_driver");
    return;
}

void Relax_Driver::init_relax(const int nat, const Input_para& inp)
{
    if (inp.calculation == "relax" || inp.calculation == "cell-relax")
    {
        if (!inp.uses_simultaneous_relaxation())
        {
            this->rl_old.init_relax(nat, inp);
        }
        else
        {
            this->rl.init_relax(nat, inp);
        }
    }
}

void Relax_Driver::iter_info(const std::vector<int>& steps, const Input_para& inp)
{
    if (inp.out_level == "ie"
            && (inp.calculation == "relax"
                || inp.calculation == "cell-relax"
                || inp.calculation == "scf"
                || inp.calculation == "nscf")
            && (inp.esolver_type != "lr"))
    {
        ModuleIO::print_screen(steps[2], steps[1], steps[0]+1);
    }

#ifdef __RAPIDJSON
    Json::init_output_array_obj();
#endif
}

void Relax_Driver::esolve(const int istep,
		ModuleESolver::ESolver* p_esolver,
		UnitCell& ucell,
		const Input_para& inp,
		ModuleBase::matrix& force,
		ModuleBase::matrix& stress,
		double& etot)
{
    p_esolver->runner(ucell, istep);

    etot = p_esolver->cal_energy();

    if (inp.cal_force)
    {
        p_esolver->cal_force(ucell, force);
    }

    if (inp.cal_stress)
    {
        p_esolver->cal_stress(ucell, stress);
    }
}

bool Relax_Driver::relax_step(std::vector<int>& steps,
		ModuleESolver::ESolver* p_esolver,
		UnitCell& ucell,
		const Input_para& inp,
		const ModuleBase::matrix& force,
		const ModuleBase::matrix& stress,
		const double etot,
		std::ofstream& ofs_running)
{
    // Guard: For non-relaxation calculations (scf, nscf, etc.), return true immediately
    // to ensure the main loop exits after one iteration. This provides robustness
    // even if relax_nmax is set to a large value.
    if (inp.calculation != "relax" && inp.calculation != "cell-relax")
    {
        return true;
    }

    bool converged = false;

    if (inp.uses_simultaneous_relaxation())
    {
        converged = this->rl.relax_step(ucell, force, stress, etot, ofs_running);
	// stress step +1
        steps[2]++;
	// fix force step to 1
        steps[1] = 1;
    }
    else
    {
        converged = this->rl_old.relax_step(steps[0]+1, etot, ucell, force,
			stress, steps[1], steps[2], ofs_running);
    }

    ModuleIO::output_after_relax(converged, p_esolver->conv_esolver, ofs_running);

    return converged;
}

void Relax_Driver::stru_out(const int istep, UnitCell& ucell, const Input_para& inp, const double etot, const ModuleBase::matrix& stress)
{
    // Guard: only output structure files for relaxation calculations
    if (inp.calculation != "relax" && inp.calculation != "cell-relax")
    {
        return;
    }

    // out_stru: -1 no output, 0 final only, 1 STRU format, 2 CIF format
    // For -1 and 0, no per-step structure output
    if (inp.out_stru <= 0)
    {
        return;
    }

    // cache global parameters to reduce repeated PARAM access
    const std::string& out_dir = PARAM.globalv.global_out_dir;
    const bool deepks_setorb = PARAM.globalv.deepks_setorb;

    // Build header comment with version, timestamp, energy and stress
    std::time_t now = std::time(nullptr);
    char time_buf[64];
    std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    std::string header = FmtCore::format("# ABACUS version: %s\n# Written at %s\n# RELAX STEP %d, Energy: %.8f eV\n",
                                          VERSION,
                                          time_buf,
                                          istep + 1,
                                          etot * ModuleBase::Ry_to_eV);
    // stress in kbar: Ry/Bohr^3 -> kbar, 3 rows
    const double stress_transform = ModuleBase::RYDBERG_SI
                                    / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI)
                                    * 1.0e-8;
    for (int i = 0; i < 3; i++)
    {
        header += FmtCore::format("# Stress (kbar): %.6f %.6f %.6f\n",
                                  stress(i, 0) * stress_transform,
                                  stress(i, 1) * stress_transform,
                                  stress(i, 2) * stress_transform);
    }

    bool need_orb = inp.basis_type == "pw";
    need_orb = need_orb && inp.init_wfc.substr(0, 3) == "nao";
    need_orb = need_orb || inp.basis_type == "lcao";
    need_orb = need_orb || inp.basis_type == "lcao_in_pw";

    const bool freq_ok = (inp.out_freq_ion > 0 && istep % inp.out_freq_ion == 0);

    // STRU_NOW: overwrite each step (for out_stru 1 and 2)
    if (inp.out_stru == 1)
    {
        unitcell::print_stru_file(ucell,
                              ucell.atoms,
                              ucell.latvec,
                              out_dir + "STRU_NOW",
                              header,
                              inp.nspin,
                              true,
                              inp.calculation == "md",
                              inp.out_mul,
                              need_orb,
                              deepks_setorb,
                              GlobalV::MY_RANK);
    }
    else if (inp.out_stru == 2)
    {
        ModuleIO::CifParser::write(out_dir + "STRU_NOW.cif",
                                   ucell,
                                   header,
                                   "data_?",
                                   GlobalV::MY_RANK);
    }

    // Numbered files per out_freq_ion (for out_stru 1 and 2 only)
    if (freq_ok)
    {
        if (inp.out_stru == 1)
        {
            unitcell::print_stru_file(ucell,
                                  ucell.atoms,
                                  ucell.latvec,
                                  out_dir + "STRU" + std::to_string(istep + 1),
                                  header,
                                  inp.nspin,
                                  true,
                                  inp.calculation == "md",
                                  inp.out_mul,
                                  need_orb,
                                  deepks_setorb,
                                  GlobalV::MY_RANK);
        }
        else if (inp.out_stru == 2)
        {
            ModuleIO::CifParser::write(out_dir + "STRU" + std::to_string(istep + 1) + ".cif",
                                       ucell,
                                       header,
                                       "data_?",
                                       GlobalV::MY_RANK);
        }
    }
}

void Relax_Driver::json_out(ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp, const ModuleBase::matrix& force, const ModuleBase::matrix& stress)
{
#ifdef __RAPIDJSON
    Json::add_output_energy(p_esolver->cal_energy() * ModuleBase::Ry_to_eV);

    double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double fac = ModuleBase::Ry_to_eV / 0.529177;
    Json::add_output_cell_coo_stress_force(&ucell, force, fac, stress, unit_transform);
#endif
}

void Relax_Driver::final_out(const int istep, UnitCell& ucell, const Input_para& inp, const double etot, const ModuleBase::matrix& stress)
{
    if (inp.calculation != "relax" && inp.calculation != "cell-relax")
    {
        return;
    }

    // out_stru: 0 no output, 1 STRU format, 2 CIF format
    // 1: write STRU_FINAL; 2: write STRU_FINAL.cif
    if (inp.out_stru == 1 || inp.out_stru == 2)
    {
        // cache global parameters to reduce repeated PARAM access
        const std::string& out_dir = PARAM.globalv.global_out_dir;
        const bool deepks_setorb = PARAM.globalv.deepks_setorb;

        // Build header comment for STRU_FINAL
        std::time_t now = std::time(nullptr);
        char time_buf[64];
        std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
        std::string header = FmtCore::format("# ABACUS version: %s\n# Written at %s\n# RELAX STEP %d (FINAL), Energy: %.8f eV\n",
                                              VERSION,
                                              time_buf,
                                              istep + 1,
                                              etot * ModuleBase::Ry_to_eV);
        const double stress_transform = ModuleBase::RYDBERG_SI
                                        / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI * ModuleBase::BOHR_RADIUS_SI)
                                        * 1.0e-8;
        for (int i = 0; i < 3; i++)
        {
            header += FmtCore::format("# Stress (kbar): %.6f %.6f %.6f\n",
                                      stress(i, 0) * stress_transform,
                                      stress(i, 1) * stress_transform,
                                      stress(i, 2) * stress_transform);
        }

        if (inp.out_stru == 1)
        {
            bool need_orb = inp.basis_type == "pw";
            need_orb = need_orb && inp.init_wfc.substr(0, 3) == "nao";
            need_orb = need_orb || inp.basis_type == "lcao";
            need_orb = need_orb || inp.basis_type == "lcao_in_pw";

            unitcell::print_stru_file(ucell,
                                      ucell.atoms,
                                      ucell.latvec,
                                      out_dir + "STRU_FINAL",
                                      header,
                                      inp.nspin,
                                      true,
                                      inp.calculation == "md",
                                      inp.out_mul,
                                      need_orb,
                                      deepks_setorb,
                                      GlobalV::MY_RANK);
        }
        else if (inp.out_stru == 2)
        {
            ModuleIO::CifParser::write(out_dir + "STRU_FINAL.cif",
                                       ucell,
                                       header,
                                       "data_?",
                                       GlobalV::MY_RANK);
        }
    }

    if (istep == inp.relax_nmax)
    {
        std::cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
        std::cout << " Geometry relaxation stops here due to reaching the maximum      " << std::endl;
        std::cout << " relaxation steps. More steps are needed to converge the results " << std::endl;
        std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
    }
    else
    {
        std::cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
        std::cout << " Geometry relaxation thresholds are reached within " << istep << " steps." << std::endl; 
        std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
    }

    if (inp.relax_nmax == 0)
    {
        std::cout << "-----------------------------------------------" << std::endl;
        std::cout << " relax_nmax = 0, DRY RUN TEST SUCCEEDS :)" << std::endl;
        std::cout << "-----------------------------------------------" << std::endl;
    }
}
