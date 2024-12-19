#include "module_base/element_name.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cube_io.h"
#include "module_io/output_log.h"
#include "write_elecstat_pot.h"

namespace ModuleIO
{

void write_elecstat_pot(
#ifdef __MPI
    const int& bz,
    const int& nbz,
#endif
    const std::string& fn,
    const int& istep,
    ModulePW::PW_Basis* rho_basis,
    const Charge* const chr,
    const UnitCell* ucell,
    const double* v_eff,
    const surchem& solvent)
{
    ModuleBase::TITLE("ModuleIO", "write_elecstat_pot");
    ModuleBase::timer::tick("ModuleIO", "write_elecstat_pot");

    std::vector<double> v_elecstat(rho_basis->nrxx, 0.0);

    const int nspin = PARAM.inp.nspin;
    const int efield = PARAM.inp.efield_flag;
    const int dip_corr = PARAM.inp.dip_cor_flag;
    const bool imp_sol = PARAM.inp.imp_sol;

    //==========================================
    // Hartree potential
    //==========================================
    ModuleBase::matrix vh(nspin, rho_basis->nrxx);
    vh = elecstate::H_Hartree_pw::v_hartree(*ucell, rho_basis, nspin, chr->rho);

    //==========================================
    //! Dipole correction
    //==========================================
    ModuleBase::matrix v_efield;
    if (efield>0 && dip_corr>0)
    {
        v_efield.create(nspin, rho_basis->nrxx);
        v_efield = elecstate::Efield::add_efield(*ucell,
                                                 const_cast<ModulePW::PW_Basis*>(rho_basis),
                                                 nspin,
                                                 chr->rho,
                                                 solvent);
    }

    //==========================================
    //! Add hartree potential and local pseudopot
    //==========================================
    for (int ir = 0; ir < rho_basis->nrxx; ir++)
    {
        // the spin index is 0
        v_elecstat[ir] = vh(0, ir) + v_eff[ir];

        if (efield>0 && dip_corr>0)
        {
            v_elecstat[ir] += v_efield(0, ir);
        }
        if(imp_sol == true)
        {
            v_elecstat[ir] += solvent.delta_phi[ir];
        }
    }

    //-------------------------------------------
    //! Get the vacuum level of the system
    //-------------------------------------------
    ModuleIO::output_vacuum_level(ucell,
                                  chr->rho,
                                  v_elecstat.data(),
                                  rho_basis->nx,
                                  rho_basis->ny,
                                  rho_basis->nz,
                                  rho_basis->nxyz,
                                  rho_basis->nrxx,
                                  rho_basis->nplane,
                                  rho_basis->startz_current);

    //-------------------------------------------
    //! Write down the electrostatic potential
    //-------------------------------------------
    int precision = 9;
    int is = -1;
    double ef_tmp = 0.0;
    int out_fermi = 0;

    ModuleIO::write_vdata_palgrid(*chr->pgrid,
        v_elecstat.data(),
        is,
        nspin,
        istep,
        fn,
        ef_tmp,
        ucell,
        precision,
        out_fermi);

    ModuleBase::timer::tick("ModuleIO", "write_elecstat_pot");
    return;
}

} // namespace ModuleIO
