#include "ctrl_output_fp.h" // use ctrl_output_fp()

#include "source_basis/module_pw/pw_basis_big.h"   // use PW_Basis_Big
#include "source_estate/elecstate.h"               // use elecstate::ElecState
#include "source_estate/module_charge/charge.h"    // use Charge
#include "source_estate/module_charge/chg_symm.h"  // use module_charge::cal_rhog_symm
#include "source_estate/write_elecstat_pot.h"      // use write_elecstat_pot
#include "source_hamilt/module_xc/xc_functional.h" // use XC_Functional
#include "source_io/module_dipole/dipole_io.h"     // use write_dipole
#include "source_io/module_elf/write_elf.h"
#include "source_io/module_output/cube_io.h"            // use write_vdata_palgrid
#include "source_io/module_parameter/input_parameter.h" // use Input_para

#ifdef __LIBXC
#include "source_io/module_chgpot/write_libxc_r.h"
#endif

namespace ModuleIO
{

void ctrl_output_fp(UnitCell& ucell,
                    const Input_para& inp,
                    elecstate::ElecState* pelec,
                    ModulePW::PW_Basis_Big* pw_big,
                    ModulePW::PW_Basis* pw_rhod,
                    Charge& chr,
                    surchem& solvent,
                    Parallel_Grid& para_grid,
                    const int istep,
                    std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_output_fp");
    ModuleBase::timer::start("ModuleIO", "ctrl_output_fp");

    const int out_app_flag = inp.out_app_flag;
    const bool gamma_only = PARAM.globalv.gamma_only_local;
    const int nspin = inp.nspin;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const std::string global_chg_dir = PARAM.globalv.global_chg_dir;
    const std::string global_pot_dir = PARAM.globalv.global_pot_dir;
    const std::string global_elf_dir = PARAM.globalv.global_elf_dir;
    const std::string global_dip_dir = PARAM.globalv.global_dip_dir;

    const int output_frequency = inp.esolver_type == "tddft" ? inp.out_freq_td : inp.out_freq_ion;
    // Print the global electronic step for TDDFT and the ionic step otherwise.
    int istep_in = -1;
    if (output_frequency > 0)
    {
        if (istep % output_frequency == 0)
        {
            istep_in = istep;
        }
    }

    std::string geom_block;
    bool should_output = (output_frequency == 0);
    if (istep_in >= 0)
    {
        // In append mode (out_app_flag == 1) each cube output reuses a fixed
        // filename, so every written step appends one complete frame to that
        // file. Modes 0 and 2 keep a per-step g{index} suffix, which makes
        // every file hold exactly one frame.
        if (out_app_flag != 1)
        {
            geom_block = "g" + std::to_string(istep + 1);
        }
        should_output = true;
    }

    // Mirror the out_wfc_lcao convention: in append mode the first ionic step
    // (istep == 0, always written because every out_freq_ion divides it)
    // truncates the fixed-name cube files and later steps append one complete
    // frame each. Single-point SCF/NSCF runs stop at istep == 0, so they keep
    // overwriting instead of stacking on a stale file from an earlier run.
    const bool append_cube = (out_app_flag == 1) && (istep > 0);

    // 4) write charge density
    if (inp.out_chg[0] > 0 && should_output)
    {
        for (int is = 0; is < nspin; ++is)
        {
            std::string fn = global_chg_dir + "chg";

            std::string spin_block;
            if (nspin == 2 || nspin == 4)
            {
                spin_block = "s" + std::to_string(is + 1);
            }
            else if (nspin == 1)
            {
                // do nothing
            }

            fn += spin_block + geom_block + ".cube";

            const std::string desc_chg = ModuleIO::make_data_desc("charge density", "magnetization density m", is, nspin);

            ModuleIO::write_vdata_palgrid(para_grid,
                                          chr.rho[is],
                                          is,
                                          nspin,
                                          istep, // change istep_in to istep, mohan 20260222
                                          fn,
                                          pelec->eferm.get_efval(is),
                                          &(ucell),
                                          inp.out_chg[1],
                                          1,
                                          PARAM.globalv.two_fermi,
                                          false,
                                          append_cube,
                                          ofs_running,
                                          desc_chg);

            if (XC_Functional::get_ked_flag())
            {
                fn = global_chg_dir + "tau";

                fn += spin_block + geom_block + ".cube";

                const std::string desc_tau = ModuleIO::make_data_desc("kinetic energy density", "kinetic energy density magnetization m", is, nspin);

                ModuleIO::write_vdata_palgrid(para_grid,
                                              chr.kin_r[is],
                                              is,
                                              nspin,
                                              istep,
                                              fn,
                                              pelec->eferm.get_efval(is),
                                              &(ucell),
                                              11, // default precision
                                              1,  // default out_fermi
                                              PARAM.globalv.two_fermi,
                                              false,
                                              append_cube,
                                              ofs_running,
                                              desc_tau);
            }
        }
    }

    // 5) write potential
    if ((inp.out_pot[0] == 1 || inp.out_pot[0] == 3) && should_output)
    {
        for (int is = 0; is < nspin; is++)
        {
            std::string fn = global_pot_dir + "pot";

            std::string spin_block;
            if (nspin == 2 || nspin == 4)
            {
                spin_block = "s" + std::to_string(is + 1);
            }
            else if (nspin == 1)
            {
                // do nothing
            }

            fn += spin_block + geom_block + ".cube";

            const std::string desc_pot = ModuleIO::make_data_desc("effective potential", "effective magnetic field b", is, nspin);

            ModuleIO::write_vdata_palgrid(para_grid,
                                          pelec->pot->get_eff_v(is),
                                          is,
                                          nspin,
                                          istep,
                                          fn,
                                          0.0, // efermi
                                          &(ucell),
                                          inp.out_pot[1], // precision
                                          0,              // out_fermi
                                          PARAM.globalv.two_fermi,
                                          false,
                                          append_cube,
                                          ofs_running,
                                          desc_pot);
        }
    }
    else if (inp.out_pot[0] == 2 && should_output)
    {
        std::string fn = global_pot_dir + "potes";
        fn += geom_block + ".cube";

        const std::vector<double>* pcc_potential = pelec->pot->pcc_electron_potential();
        const std::vector<double>* solvent_potential = pelec->pot->find_sol_pot();
        ModuleIO::write_elecstat_pot(
#ifdef __MPI
            pw_big->bz,
            pw_big->nbz,
#endif
            fn,
            istep,
            pw_rhod,
            &chr,
            &(ucell),
            pelec->pot->get_fixed_v(),
            solvent,
            pcc_potential,
            solvent_potential,
            inp.out_pot[1],
            append_cube,
            ofs_running);
    }

    // 6) write ELF
    if (inp.out_elf[0] > 0 && should_output)
    {
        // write_elf() consumes chr.kin_r; the final tau written by sum_band
        // must be symmetrized before output (the tau buffer is allocated
        // because ELF output was requested).
        for (int is = 0; is < nspin; is++)
        {
            module_charge::cal_rhog_symm(is, chr, pw_rhod, ucell.symm);
        }

        std::string out_dir = global_elf_dir;
        ModuleIO::write_elf(out_dir,
                            istep,
                            nspin,
                            chr.rho,
                            chr.kin_r,
                            pw_rhod,
                            para_grid,
                            &(ucell),
                            inp.out_elf[1],
                            geom_block,
                            PARAM.globalv.two_fermi,
                            append_cube,
                            ofs_running);
    }

#ifdef __LIBXC
    // 7) write xc(r)
    if (inp.out_xc_r[0] >= 0 && should_output)
    {
        ModuleIO::write_libxc_r(inp.out_xc_r[0],
                                XC_Functional::get_func_id(),
                                pw_rhod->nrxx, // number of real-space grid
                                ucell.omega,   // volume of cell
                                ucell.tpiba,
                                chr,
                                *pw_big,
                                *pw_rhod);
    }
#endif

    // RT-TDDFT response analysis requires every dipole sample, independently of snapshot output.
    const bool output_dipole = inp.esolver_type == "tddft" || should_output;
    if (inp.out_dipole == 1 && output_dipole)
    {
        for (int is = 0; is < nspin; ++is)
        {
            std::stringstream ss_dipole;
            ss_dipole << global_dip_dir << "dipole_s" << is + 1 << ".txt";
            ModuleIO::write_dipole(ucell, chr.rho[is], pw_rhod, istep, ss_dipole.str(), ofs_running);
        }
    }

    ModuleBase::timer::end("ModuleIO", "ctrl_output_fp");
}

} // namespace ModuleIO
