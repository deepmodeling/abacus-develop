#include "module_base/element_name.h"
#include "module_base/timer.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cube_io.h"
#include "potential_io.h"

namespace ModuleIO
{

void write_potential(
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
    ModuleBase::TITLE("potential", "write_potential");
    if (GlobalV::out_pot != 1)
    {
        return;
    }
    ModuleBase::timer::tick("Potential", "write_potential");

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

    ModuleBase::timer::tick("Potential", "write_potential");
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