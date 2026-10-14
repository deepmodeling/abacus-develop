#include "source_base/element_name.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#include "source_estate/module_pot/electrostatic_potential.h"
#include "source_io/module_output/cube_io.h"
#include "source_io/module_output/output_log.h"
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
    const surchem& solvent,
    const std::vector<double>* pcc_potential,
    const int precision,
    std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_elecstat_pot");
    ModuleBase::timer::start("ModuleIO", "write_elecstat_pot");

    auto v_elecstat = elecstate::electrostatic_potential(*ucell, *rho_basis, *chr,
                                                       v_eff, solvent, PARAM.inp, pcc_potential);

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
                                  rho_basis->startz_current,
                                  ofs_running);

    //-------------------------------------------
    //! Write down the electrostatic potential
    //-------------------------------------------
    int is = -1;
    double ef_tmp = 0.0;
    int out_fermi = 0;

    ModuleIO::write_vdata_palgrid(*chr->pgrid,
        v_elecstat.data(),
        is,
        PARAM.inp.nspin,
        istep,
        fn,
        ef_tmp,
        ucell,
        precision,
        out_fermi,
        PARAM.globalv.two_fermi,
        false,
        ofs_running,
        "electrostatic potential");

    ModuleBase::timer::end("ModuleIO", "write_elecstat_pot");
    return;
}

} // namespace ModuleIO
