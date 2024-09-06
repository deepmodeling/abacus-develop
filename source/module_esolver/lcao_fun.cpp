#include "esolver_ks_lcao.h"

#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_io/dos_nao.h"
#include "module_io/nscf_band.h"
#include "module_io/output_dmk.h"
#include "module_io/output_log.h"
#include "module_io/output_mulliken.h"
#include "module_io/output_sk.h"
#include "module_io/to_qo.h"
#include "module_io/write_HS.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_proj_band_lcao.h"
#include "module_parameter/parameter.h"

//--------------temporary----------------------------
#include <memory>

#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h" // need divide_HS_in_frag
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"

//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
// function used by deepks
#include "module_elecstate/cal_dm.h"
//---------------------------------------------------

#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_io/io_dmk.h"
#include "module_io/write_dmr.h"
#include "module_io/write_wfc_nao.h"

namespace ModuleESolver
{

//------------------------------------------------------------------------------
//! the 18th function of ESolver_KS_LCAO: create_Output_Mat_Sparse
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
ModuleIO::Output_Mat_Sparse<TK> ESolver_KS_LCAO<TK, TR>::create_Output_Mat_Sparse(int istep)
{
    return ModuleIO::Output_Mat_Sparse<TK>(hsolver::HSolverLCAO<TK>::out_mat_hsR,
        hsolver::HSolverLCAO<TK>::out_mat_dh,
        hsolver::HSolverLCAO<TK>::out_mat_t,
        PARAM.inp.out_mat_r,
        istep,
        this->pelec->pot->get_effective_v(),
        this->pv,
        this->GK, // mohan add 2024-04-01
        two_center_bundle_,
        orb_,
        GlobalC::GridD, // mohan add 2024-04-06
        this->kv,
        this->p_hamilt);
}

//------------------------------------------------------------------------------
//! the 19th function of ESolver_KS_LCAO: md_skip_out
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template <typename TK, typename TR>
bool ESolver_KS_LCAO<TK, TR>::md_skip_out(std::string calculation, int istep, int interval)
{
    if (calculation == "md")
    {
        if (istep % interval != 0)
        {
            return true;
        }
    }
    return false;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_mag(const int istep, const bool print)
{
    auto cell_index = CellIndex(GlobalC::ucell.get_atomLabels(),
                                GlobalC::ucell.get_atomCounts(),
                                GlobalC::ucell.get_lnchiCounts(),
                                GlobalV::NSPIN);
    auto out_sk = ModuleIO::Output_Sk<TK>(this->p_hamilt,
                                          &(this->pv),
                                          GlobalV::NSPIN,
                                          this->kv.get_nks());
    auto out_dmk = ModuleIO::Output_DMK<TK>(dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                            &(this->pv),
                                            GlobalV::NSPIN,
                                            this->kv.get_nks());
    auto mulp = ModuleIO::Output_Mulliken<TK>(&(out_sk),
                                              &(out_dmk),
                                              &(this->pv),
                                              &cell_index,
                                              this->kv.isk,
                                              GlobalV::NSPIN);
    auto atom_chg = mulp.get_atom_chg();
    /// used in updating mag info in STRU file
    GlobalC::ucell.atom_mulliken = mulp.get_atom_mulliken(atom_chg);
    if (print && GlobalV::MY_RANK == 0)
    {
        /// write the Orbital file
        cell_index.write_orb_info(GlobalV::global_out_dir);
        /// write mulliken.txt
        mulp.write(istep, GlobalV::global_out_dir);
        /// write atomic mag info in running log file
        mulp.print_atom_mag(atom_chg, GlobalV::ofs_running);
    }
}

//------------------------------------------------------------------------------
//! the 20th,21th,22th functions of ESolver_KS_LCAO
//! mohan add 2024-05-11
//------------------------------------------------------------------------------
template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
