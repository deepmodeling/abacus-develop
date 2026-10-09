#include "chg_atomic.h"
#include "chg_atomic_detail.h"

#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_cell/unitcell.h"
#include "source_cell/magnetism.h"

#include <iomanip>

namespace module_charge
{

void atomic_rho(const int spin_number_need,
                const double& omega,
                double** rho_in,
                const ModuleBase::ComplexMatrix& strucFac,
                const UnitCell& ucell,
                const ModulePW::PW_Basis* rhopw,
                const AtomicRhoCfg& cfg)
{
    ModuleBase::TITLE("module_charge", "atomic_rho");
    ModuleBase::timer::start("module_charge", "atomic_rho");

    std::ostream& ofs_warning = cfg.ofs_warning;
    const int test_charge = cfg.test_charge;
    const bool domag = cfg.domag;
    const bool domag_z = cfg.domag_z;

    ModuleBase::ComplexMatrix rho_g3d(spin_number_need, rhopw->npw);

    for (int it = 0; it < ucell.ntype; it++)
    {
        // check the start magnetization
        const int startmag_type = (ucell.magnet.start_mag[it] != 0.0) ? 1 : 2;
        ofs_warning << " " << std::setw(40) << "startmag_type"
                    << " = " << startmag_type << std::endl;

        const Atom* const atom = &ucell.atoms[it];

        if (!atom->flag_empty_element) // Peize Lin add for bsse 2021.04.07
        {
            const int mesh = atom->ncpp.msh;
            const std::vector<double> rhoatm
                = detail::compute_rhoatm(*atom, mesh, ofs_warning);
            const std::vector<double> rho_lgl
                = detail::compute_rho_lgl(*atom, rhopw, ucell, rhoatm,
                                          test_charge, omega, ofs_warning);

            detail::RhoG3dCtx ctx{rho_g3d, strucFac, rho_lgl, rhopw, it};

            if (spin_number_need == 1)
            {
                detail::fill_rho_g3d_nspin1(ctx);
            }
            else if (spin_number_need == 2)
            {
                detail::fill_rho_g3d_nspin2(ctx, startmag_type,
                                            ucell.magnet.start_mag[it], *atom);
            }
            else if (spin_number_need == 4)
            {
                if (startmag_type == 1)
                {
                    detail::fill_rho_g3d_nspin4_type1(ctx,
                                                       ucell.magnet.start_mag[it],
                                                       *atom, domag, domag_z);
                }
                else
                {
                    detail::fill_rho_g3d_nspin4_type2(ctx, *atom, domag, domag_z);
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("module_charge::atomic_rho",
                                          " Either 1 or 2 or 4, check SPIN number !");
            }
        }
    }

    detail::normalize_and_check(rho_in, rho_g3d, rhopw, spin_number_need,
                               omega, ofs_warning, cfg.nelec);

    ModuleBase::timer::end("module_charge", "atomic_rho");
}

} // namespace module_charge
