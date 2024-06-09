#include "potential_io.h"

#include "module_base/element_name.h"
#include "module_base/timer.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cube_io.h"

namespace ModuleIO
{

void write_pot(
    const int &out_pot,
    const int &nspin, 
    const std::string &global_out_dir,
#ifdef __MPI
    const int& bz,
    const int& nbz,
    const int& nplane,
    const int& startz_current,
#endif
    const int& nx,
    const int& ny,
    const int& nz,
    const ModuleBase::matrix& v)
{
    ModuleBase::TITLE("ModuleIO", "write_pot");
    if(out_pot == 3)
    {
        for(int is = 0; is < nspin; is++)
        {
            std::stringstream ss;
            ss << global_out_dir << "SPIN" << is+1 << "_POT_INI.cube";
            ModuleIO::write_pot_spin(
                    out_pot,
#ifdef __MPI
					bz,
					nbz,
					nplane,
					startz_current,
#endif
                    is,
                    0, // iter
                    ss.str(),
                    nx,
                    ny,
                    nz,
                    v,
                    11); // precsion
        }
    }

    ModuleBase::TITLE("ModuleIO", "write_pot");
    return;
}



void write_pot_spin(
    const int& out_pot,
#ifdef __MPI
    const int& bz,
    const int& nbz,
    const int& nplane,
    const int& startz_current,
#endif
    const int& is,
    const int& iter,
    const std::string& fn,
    const int& nx,
    const int& ny,
    const int& nz,
    const ModuleBase::matrix& v,
    const int& precision,
    const int& hartree)
{
    ModuleBase::TITLE("ModuleIO", "write_pot_spin");
    if (out_pot != 1 && out_pot != 3)
    {
        return;
    }
    ModuleBase::timer::tick("ModuleIO", "write_pot_spin");

    double* temp_v = nullptr;
    if (is == 0)
    {
        temp_v = v.c;
    }
    else if (is == 1)
    {
        temp_v = &(v.c[nx * ny * nz]);
    }

    double ef_tmp = 0.;
    int out_fermi = 0;

    ModuleIO::write_cube(
#ifdef __MPI
        bz,
        nbz,
        nplane,
        startz_current,
#endif
        temp_v,
        is,
        GlobalV::NSPIN,
        iter,
        fn,
        nx,
        ny,
        nz,
        ef_tmp,
        &(GlobalC::ucell),
        precision,
        out_fermi);

    ModuleBase::timer::tick("ModuleIO", "write_pot_spin");
    return;
}

void write_elecstat_pot(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& fn,
    ModulePW::PW_Basis* rho_basis,
    const Charge* const chr,
    const UnitCell* ucell_,
    const double* v_effective_fixed)
{
    ModuleBase::TITLE("Potential", "write_elecstat_pot");
    ModuleBase::timer::tick("Potential", "write_elecstat_pot");

    double* v_elecstat = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(v_elecstat, rho_basis->nrxx);

    //==========================================
    // Hartree potential
    //==========================================
    ModuleBase::matrix vh(GlobalV::NSPIN, rho_basis->nrxx);
    vh = elecstate::H_Hartree_pw::v_hartree(*ucell_, rho_basis, GlobalV::NSPIN, chr->rho);

    //==========================================
    // Dipole correction
    //==========================================
    ModuleBase::matrix v_efield;
    if (GlobalV::EFIELD_FLAG && GlobalV::DIP_COR_FLAG)
    {
        v_efield.create(GlobalV::NSPIN, rho_basis->nrxx);
        v_efield = elecstate::Efield::add_efield(*(ucell_),
                                                 const_cast<ModulePW::PW_Basis*>(rho_basis),
                                                 GlobalV::NSPIN,
                                                 chr->rho,
                                                 GlobalC::solvent_model);
    }

    //==========================================
    // Add hartree potential and local pseudopot
    //==========================================
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        v_elecstat[ir] = vh(0, ir) + v_effective_fixed[ir];

        if (GlobalV::EFIELD_FLAG && GlobalV::DIP_COR_FLAG)
        {
            v_elecstat[ir] += v_efield(0, ir);
        }
        if (GlobalV::imp_sol)
        {
            v_elecstat[ir] += GlobalC::solvent_model.delta_phi[ir];
        }
    }

    //-------------------------------------------
    // output the electrostatic potential into a file.
    //-------------------------------------------
    int precision = 9;
    int is = -1;
    double ef_tmp = 0.;
    int out_fermi = 0;
    ModuleIO::write_cube(
#ifdef __MPI
        bz,
        nbz,
        rho_basis->nplane,
        rho_basis->startz_current,
#endif
        v_elecstat,
        is,
        GlobalV::NSPIN,
        0,
        fn,
        rho_basis->nx,
        rho_basis->ny,
        rho_basis->nz,
        ef_tmp,
        &(GlobalC::ucell),
        precision,
        out_fermi);

    delete[] v_elecstat;

    ModuleBase::timer::tick("Potential", "write_elecstat_pot");
    return;
}

} // namespace ModuleIO
