#include "relax_driver.h"
#include "source_base/global_file.h"
#include "source_io/module_output/cif_io.h"
#include "source_io/module_json/output_info.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_output/print_info.h"
#include "source_io/module_output/read_exit_file.h"
#include "source_io/module_parameter/parameter.h"
#include "source_cell/print_cell.h"

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
    while (steps[0] < inp.relax_nmax)
    {
        ModuleBase::matrix force(ucell.nat, 3);
        ModuleBase::matrix stress(3, 3);
        double etot = 0.0;

        this->iter_info(steps, inp);
        this->esolve(steps[0], p_esolver, ucell, inp, force, stress, etot);
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

    this->final_out(steps[0], ucell, inp);

    ModuleBase::timer::end("Relax_Driver", "relax_driver");
    return;
}

void Relax_Driver::init_relax(const int nat, const Input_para& inp)
{
    if (inp.calculation == "relax" || inp.calculation == "cell-relax")
    {
        if (!inp.relax_new)
        {
            this->rl_old.init_relax(nat);
        }
        else
        {
            this->rl.init_relax(nat);
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

    if (inp.relax_new)
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

    this->stru_out(steps[0], ucell, inp);

    ModuleIO::output_after_relax(converged, p_esolver->conv_esolver, ofs_running);

    return converged;
}

void Relax_Driver::stru_out(const int istep, UnitCell& ucell, const Input_para& inp)
{
    bool need_orb = inp.basis_type == "pw";
    need_orb = need_orb && inp.init_wfc.substr(0, 3) == "nao";
    need_orb = need_orb || inp.basis_type == "lcao";
    need_orb = need_orb || inp.basis_type == "lcao_in_pw";

    std::stringstream ss, ss1;
    ss << PARAM.globalv.global_out_dir << "STRU_ION_D";

    unitcell::print_stru_file(ucell,
                          ucell.atoms,
                          ucell.latvec,
                          ss.str(),
                          inp.nspin,
                          true,
                          inp.calculation == "md",
                          inp.out_mul,
                          need_orb,
                          PARAM.globalv.deepks_setorb,
                          GlobalV::MY_RANK);

    if (inp.out_stru)
    {
        if (inp.out_freq_ion == 0 || istep % inp.out_freq_ion == 0)
        {
            ss1 << PARAM.globalv.global_out_dir << "STRU_ION";
            ss1 << istep+1 << "_D";
            unitcell::print_stru_file(ucell,
                                  ucell.atoms,
                                  ucell.latvec,
                                  ss1.str(),
                                  inp.nspin,
                                  true,
                                  inp.calculation == "md",
                                  inp.out_mul,
                                  need_orb,
                                  PARAM.globalv.deepks_setorb,
                                  GlobalV::MY_RANK);

            ModuleIO::CifParser::write(PARAM.globalv.global_out_dir + "STRU_NOW.cif",
                                       ucell,
                                       "# Generated by ABACUS ModuleIO::CifParser",
                                       "data_?");
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

void Relax_Driver::final_out(const int istep, UnitCell& ucell, const Input_para& inp)
{
    if (inp.calculation != "relax" && inp.calculation != "cell-relax")
    {
        return;
    }

    ModuleIO::CifParser::write(PARAM.globalv.global_out_dir + "STRU_FINAL.cif",
                               ucell,
                               "# Generated by ABACUS ModuleIO::CifParser",
                               "data_?");

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
