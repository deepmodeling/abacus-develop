#include "ctrl_output_fp.h" // use ctrl_output_fp()

#include "../module_output/cube_io.h"                                  // use write_vdata_palgrid
#include "source_estate/module_charge/symmetry_rho.h" // use Symmetry_rho
#include "source_hamilt/module_xc/xc_functional.h"    // use XC_Functional
#include "source_io/module_chgpot/write_elecstat_pot.h"             // use write_elecstat_pot
#include "source_io/module_elf/write_elf.h"

#ifdef USE_LIBXC
#include "source_io/module_chgpot/write_libxc_r.h"
#endif

namespace ModuleIO
{

void ctrl_output_fp(UnitCell& ucell,
                    elecstate::ElecState* pelec,
                    ModulePW::PW_Basis_Big* pw_big,
                    ModulePW::PW_Basis* pw_rhod,
                    Charge& chr,
                    surchem& solvent,
                    Parallel_Grid& para_grid,
                    const int istep)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_output_fp");
    ModuleBase::timer::start("ModuleIO", "ctrl_output_fp");

    const bool out_app_flag = PARAM.inp.out_app_flag;
    const bool gamma_only = PARAM.globalv.gamma_only_local;
    const int nspin = PARAM.inp.nspin;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;

    // print out the 'g' index when istep_in != -1
    int istep_in = -1;
    if (PARAM.inp.esolver_type != "tddft" && PARAM.inp.out_freq_ion > 0) // default value of out_freq_ion is 0
    {
        if (istep % PARAM.inp.out_freq_ion == 0)
        {
            istep_in = istep;
        }
    }
    else if (PARAM.inp.esolver_type == "tddft" && PARAM.inp.out_freq_td > 0) // default value of out_freq_td is 0
    {
        if (istep % PARAM.inp.out_freq_td == 0)
        {
            istep_in = istep;
        }
    }

    std::string geom_block;
    if (istep_in == -1)
    {
        // do nothing
    }
    else if (istep_in >= 0)
    {
        geom_block = "g" + std::to_string(istep + 1);
    }

    // 4) write charge density
    if (PARAM.inp.out_chg[0] > 0)
    {
        for (int is = 0; is < nspin; ++is)
        {
            pw_rhod->real2recip(chr.rho_save[is], chr.rhog_save[is]);

            std::string fn = PARAM.globalv.global_out_dir + "chg";

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

            ModuleIO::write_vdata_palgrid(para_grid,
                                          chr.rho_save[is],
                                          is,
                                          nspin,
                                          istep, // change istep_in to istep, mohan 20260222
                                          fn,
                                          pelec->eferm.get_efval(is),
                                          &(ucell),
                                          PARAM.inp.out_chg[1],
                                          1);

            if (XC_Functional::get_ked_flag())
            {
                fn = PARAM.globalv.global_out_dir + "tau";

                fn += spin_block + geom_block + ".cube";

                ModuleIO::write_vdata_palgrid(para_grid,
                                              chr.kin_r_save[is],
                                              is,
                                              nspin,
                                              istep,
                                              fn,
                                              pelec->eferm.get_efval(is),
                                              &(ucell));
            }
        }
    }

    // 5) write potential
    if (PARAM.inp.out_pot[0] == 1 || PARAM.inp.out_pot[0] == 3)
    {
        for (int is = 0; is < nspin; is++)
        {
            std::string fn = PARAM.globalv.global_out_dir + "pot";

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

            ModuleIO::write_vdata_palgrid(para_grid,
                                          pelec->pot->get_eff_v(is),
                                          is,
                                          nspin,
                                          istep_in,
                                          fn,
                                          0.0, // efermi
                                          &(ucell),
                                          PARAM.inp.out_pot[1],  // precision
                                          0); // out_fermi
        }
    }
    else if (PARAM.inp.out_pot[0] == 2)
    {
        std::string fn = PARAM.globalv.global_out_dir + "potes";
        fn += geom_block + ".cube";

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
	    PARAM.inp.out_pot[1]);
    }

    // 6) write ELF
    if (PARAM.inp.out_elf[0] > 0)
    {
        chr.cal_elf = true;
        Symmetry_rho srho;
        for (int is = 0; is < nspin; is++)
        {
            srho.begin(is, chr, pw_rhod, ucell.symm);
        }

        std::string out_dir = PARAM.globalv.global_out_dir;
        ModuleIO::write_elf(
#ifdef __MPI
            pw_big->bz,
            pw_big->nbz,
#endif
            out_dir,
            istep,
            nspin,
            chr.rho,
            chr.kin_r,
            pw_rhod,
            para_grid,
            &(ucell),
            PARAM.inp.out_elf[1]);
    }

#ifdef USE_LIBXC
    // 7) write xc(r)
    if (PARAM.inp.out_xc_r[0] >= 0)
    {
        ModuleIO::write_libxc_r(PARAM.inp.out_xc_r[0],
                                XC_Functional::get_func_id(),
                                pw_rhod->nrxx, // number of real-space grid
                                ucell.omega,   // volume of cell
                                ucell.tpiba,
                                chr,
                                *pw_big,
                                *pw_rhod);
    }
#endif

    ModuleBase::timer::end("ModuleIO", "ctrl_output_fp");
}

} // namespace ModuleIO
