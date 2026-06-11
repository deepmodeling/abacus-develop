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
        const Input_para& inp)
{ 
    ModuleBase::TITLE("Relax_Driver", "relax_driver");
    ModuleBase::timer::start("Relax_Driver", "relax_driver");

    this->init_relax(ucell.nat, inp);

    this->istep = 1;

    // Main iteration loop for relaxation calculations
    // For scf/nscf calculations, relax_step returns true immediately,
    // so the loop exits after one iteration
    while (this->istep <= inp.relax_nmax)
    {
        this->iter_info(inp);
        this->esolve(p_esolver, ucell);
        bool converged = this->relax_step(p_esolver, ucell, inp);
        this->json_out(p_esolver, ucell, inp);

        // Check stop conditions
        if (converged)
        {
            // Relaxation converged, exit loop immediately
            break;
        }
        else if (ModuleIO::read_exit_file(GlobalV::MY_RANK, "EXIT", GlobalV::ofs_running))
        {
            // EXIT file detected, exit loop
            break;
        }

        ++this->istep;
    }

    this->final_out(ucell, inp);

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

void Relax_Driver::iter_info(const Input_para& inp)
{
    if (inp.out_level == "ie"
            && (inp.calculation == "relax" 
                || inp.calculation == "cell-relax" 
                || inp.calculation == "scf"
                || inp.calculation == "nscf")
            && (inp.esolver_type != "lr"))
    {
        ModuleIO::print_screen(this->stress_step, this->force_step, this->istep);
    }

#ifdef __RAPIDJSON
    Json::init_output_array_obj();
#endif
}

void Relax_Driver::esolve(ModuleESolver::ESolver* p_esolver, UnitCell& ucell)
{
    p_esolver->runner(ucell, this->istep - 1);

    this->etot = p_esolver->cal_energy();

    if (PARAM.inp.cal_force)
    {
        p_esolver->cal_force(ucell, this->force_);
    }

    if (PARAM.inp.cal_stress)
    {
        p_esolver->cal_stress(ucell, this->stress_);
    }
}

bool Relax_Driver::relax_step(ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp)
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
        converged = this->rl.relax_step(ucell, this->force_, this->stress_, this->etot);
        this->stress_step = this->istep + 1;
        this->force_step = 1;
    }
    else
    {
        converged = this->rl_old.relax_step(this->istep, this->etot, ucell, this->force_, 
			this->stress_, this->force_step, this->stress_step);
    }

    this->stru_out(ucell, inp);

    ModuleIO::output_after_relax(converged, p_esolver->conv_esolver, GlobalV::ofs_running);

    return converged;
}

void Relax_Driver::stru_out(UnitCell& ucell, const Input_para& inp)
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

    if (Ions_Move_Basic::out_stru)
    {
        ss1 << PARAM.globalv.global_out_dir << "STRU_ION";
        ss1 << this->istep << "_D";
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

void Relax_Driver::json_out(ModuleESolver::ESolver* p_esolver, UnitCell& ucell, const Input_para& inp)
{
#ifdef __RAPIDJSON
    Json::add_output_energy(p_esolver->cal_energy() * ModuleBase::Ry_to_eV);

    double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double fac = ModuleBase::Ry_to_eV / 0.529177;
    Json::add_output_cell_coo_stress_force(&ucell, this->force_, fac, this->stress_, unit_transform);
#endif
}

void Relax_Driver::final_out(UnitCell& ucell, const Input_para& inp)
{
    if (inp.calculation != "relax" && inp.calculation != "cell-relax")
    {
        return;
    }

    ModuleIO::CifParser::write(PARAM.globalv.global_out_dir + "STRU_FINAL.cif",
                               ucell,
                               "# Generated by ABACUS ModuleIO::CifParser",
                               "data_?");

    if (this->istep - 1 == inp.relax_nmax)
    {
        std::cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
        std::cout << " Geometry relaxation stops here due to reaching the maximum      " << std::endl;
        std::cout << " relaxation steps. More steps are needed to converge the results " << std::endl;
        std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
    }
    else
    {
        std::cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
        std::cout << " Geometry relaxation thresholds are reached within " << this->istep - 1 << " steps." << std::endl; 
        std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl; 
    }

    if (inp.relax_nmax == 0)
    {
        std::cout << "-----------------------------------------------" << std::endl;
        std::cout << " relax_nmax = 0, DRY RUN TEST SUCCEEDS :)" << std::endl;
        std::cout << "-----------------------------------------------" << std::endl;
    }
}
