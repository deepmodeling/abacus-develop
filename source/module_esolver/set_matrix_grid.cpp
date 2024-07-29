#include "esolver_ks_lcao.h"

#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_esolver/esolver_ks_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_parameter/parameter.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_io/rho_io.h"
#include "module_io/write_pot.h"
#include "module_io/write_wfc_nao.h"
#include "module_io/read_wfc_nao.h"


namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::set_matrix_grid(Record_adj& ra)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "set_matrix_grid");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "set_matrix_grid");

    // (1) Find adjacent atoms for each atom.
    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                     GlobalV::OUT_LEVEL,
                                                     GlobalC::ORB.get_rcutmax_Phi(),
                                                     GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                                                     GlobalV::GAMMA_ONLY_LOCAL);

    atom_arrange::search(PARAM.inp.search_pbc,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         GlobalC::ucell,
                         GlobalV::SEARCH_RADIUS,
                         GlobalV::test_atom_input);

    // ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"SEARCH ADJACENT
    // ATOMS");

    // (3) Periodic condition search for each grid.
    double dr_uniform = 0.001;
    std::vector<double> rcuts;
    std::vector<std::vector<double>> psi_u;
    std::vector<std::vector<double>> dpsi_u;
    std::vector<std::vector<double>> d2psi_u;

    Gint_Tools::init_orb(dr_uniform, rcuts, GlobalC::ucell, psi_u, dpsi_u, d2psi_u);

    this->GridT.set_pbc_grid(this->pw_rho->nx,
                             this->pw_rho->ny,
                             this->pw_rho->nz,
                             this->pw_big->bx,
                             this->pw_big->by,
                             this->pw_big->bz,
                             this->pw_big->nbx,
                             this->pw_big->nby,
                             this->pw_big->nbz,
                             this->pw_big->nbxx,
                             this->pw_big->nbzp_start,
                             this->pw_big->nbzp,
                             this->pw_rho->ny,
                             this->pw_rho->nplane,
                             this->pw_rho->startz_current,
                             GlobalC::ucell,
                             dr_uniform,
                             rcuts,
                             psi_u,
                             dpsi_u,
                             d2psi_u,
                             PARAM.inp.nstream);
    psi_u.clear();
    psi_u.shrink_to_fit();
    dpsi_u.clear();
    dpsi_u.shrink_to_fit();
    d2psi_u.clear();
    d2psi_u.shrink_to_fit();

    // (2)For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    ra.for_2d(this->pv, GlobalV::GAMMA_ONLY_LOCAL);

    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        // need to first calculae lgd.
        // using GridT.init.
        this->GridT.cal_nnrg(&this->pv);
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "set_matrix_grid");
    return;
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;

}
