#include "source_pw/module_pwdft/td_pw.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/cell_tools.h"

#include <algorithm>

namespace pw
{

void check_td_input(const Input_para& input, const UnitCell& cell, const bool needs_ked, const bool cal_exx)
{
    if (input.init_vecpot_file)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT does not support init_vecpot_file=true; set init_vecpot_file=false.");
    }
    if (input.out_current != 0 && input.out_current != 1)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT requires out_current=0 or 1.");
    }
    if (input.scf_nmax < 2)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT",
                                 "PW RT-TDDFT requires scf_nmax >= 2: one predictor and at least one corrector iteration.");
    }
    if (unitcell::if_atoms_can_move(cell.atoms, cell.ntype))
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT",
                                 "PW RT-TDDFT requires fixed ions; set every atom's movement flags to 0 0 0 in STRU.");
    }
    if (input.calculation == "cell-relax" || (input.calculation == "md" && (input.mdp.md_type == "npt" || input.mdp.md_type == "msst")))
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT requires a fixed cell; cell-relax, NPT and MSST are not supported.");
    }
    if (input.socket_driver)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT",
                                 "PW RT-TDDFT requires fixed ions and a fixed cell; socket_driver is not supported.");
    }
    if (input.dft_plus_u != 0)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT does not support DFT+U; set dft_plus_u=0.");
    }
    if (input.sc_mag_switch)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT does not support spin constraints; set sc_mag_switch=false.");
    }
    if (cal_exx)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT does not support exact exchange (EXX); use a non-hybrid functional.");
    }
    if ((input.nspin != 1 && input.nspin != 2) || input.bndpar != 1 || (input.td_stype != 0 && input.td_stype != 1))
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT requires nspin=1 or 2, bndpar=1 and td_stype=0 or 1.");
    }
    if (input.estep_per_md != 1)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT currently requires estep_per_md to be 1.");
    }
    if (input.mdp.md_restart || input.restart_load)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT", "PW RT-TDDFT restart is not supported yet.");
    }
    if (input.td_stype == 1 && needs_ked)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW_TDDFT",
                                 "PW RT-TDDFT with a kinetic-energy-density functional requires td_stype=0; "
                                 "the velocity-gauge coupling of the functional is not implemented.");
    }
}

double td_momentum_bound(const ModulePW::PW_Basis_K& basis, const double tpiba)
{
    ModuleBase::timer::start("pw", "td_momentum_bound");
    double bound = 0.0;
    for (int ik = 0; ik < basis.nks; ++ik)
    {
        for (int ig = 0; ig < basis.npwk[ik]; ++ig)
        {
            bound = std::max(bound, basis.getgpluskcar(ik, ig).norm() * tpiba);
        }
    }
    Parallel_Reduce::reduce_max(bound);
    ModuleBase::timer::end("pw", "td_momentum_bound");
    return bound;
}

void ensure_td_vnl(const UnitCell& cell,
                   const double unshifted_bound,
                   const ModuleBase::Vector3<double>& A_right_ha,
                   const ModuleBase::Vector3<double>& A_prop_ha,
                   pseudopot_cell_vnl* projectors)
{
    ModuleBase::timer::start("pw", "ensure_td_vnl");
    const double field_bound = std::max(A_right_ha.norm(), A_prop_ha.norm());
    projectors->ensure_vnl_range(cell, unshifted_bound + field_bound);
    ModuleBase::timer::end("pw", "ensure_td_vnl");
}

} // namespace pw
