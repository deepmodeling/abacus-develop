#include "veff_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

namespace hamilt
{

template class Veff<OperatorLCAO<double, double>>;

template class Veff<OperatorLCAO<std::complex<double>, double>>;

template class Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>;

// initialize_HR()
template <typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::initialize_HR(const UnitCell* ucell_in,
                                        Grid_Driver* GridD,
                                        const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("Veff", "initialize_HR");
    ModuleBase::timer::tick("Veff", "initialize_HR");

    for (int iat1 = 0; iat1 < ucell_in->nat; iat1++)
    {
        auto tau1 = ucell_in->get_tau(iat1);
        int T1, I1;
        ucell_in->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell_in, tau1, T1, I1, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell_in->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad1];
            // choose the real adjacent atoms
            const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell_in->cal_dtau(iat1, iat2, R_index2).norm() * ucell_in->lat0
                < orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut())
            {
                hamilt::AtomPair<TR> tmp(iat1, iat2, R_index2.x, R_index2.y, R_index2.z, paraV);
                this->hR->insert_pair(tmp);
            }
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    this->hR->allocate(nullptr, true);

    ModuleBase::timer::tick("Veff", "initialize_HR");
}



template<typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::tick("Veff", "contributeHR");
    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    double* vr_eff1 = this->pot->get_effective_v(GlobalV::CURRENT_SPIN);
    double* vofk_eff1 = this->pot->get_effective_vofk(GlobalV::CURRENT_SPIN);

    //--------------------------------------------
    //(2) check if we need to calculate
    // pvpR = < phi0 | v(spin) | phiR> for a new spin.
    //--------------------------------------------
    // GlobalV::ofs_running << " (spin change)" << std::endl;
    this->GK->reset_spin(GlobalV::CURRENT_SPIN);

    // if you change the place of the following code,
    // rememeber to delete the #include
    if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
    {
        Gint_inout inout(vr_eff1, vofk_eff1, 0, Gint_Tools::job_type::vlocal_meta);
        this->GK->cal_gint(&inout);
    }
    else
    {
        // vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
        Gint_inout inout(vr_eff1, 0, Gint_Tools::job_type::vlocal);
        this->GK->cal_gint(&inout);
    }

    // added by zhengdy-soc, for non-collinear case
    // integral 4 times, is there any method to simplify?
    if (GlobalV::NSPIN == 4)
    {
        for (int is = 1; is < 4; is++)
        {
            vr_eff1 = this->pot->get_effective_v(is);
            if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
            {
                vofk_eff1 = this->pot->get_effective_vofk(is);
            }
            
            if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
            {
                Gint_inout inout(vr_eff1, vofk_eff1, is, Gint_Tools::job_type::vlocal_meta);
                this->GK->cal_gint(&inout);
            }
            else
            {
                Gint_inout inout(vr_eff1, is, Gint_Tools::job_type::vlocal);
                this->GK->cal_gint(&inout);
            }
        }
    }

    this->GK->transfer_pvpR(this->hR);

    ModuleBase::timer::tick("Veff", "contributeHR");
    return;
}

// special case of gamma-only
template<>
void Veff<OperatorLCAO<double, double>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::tick("Veff", "contributeHR");

    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    const double* vr_eff1 = this->pot->get_effective_v(GlobalV::CURRENT_SPIN);
    const double* vofk_eff1 = this->pot->get_effective_vofk(GlobalV::CURRENT_SPIN);

    //--------------------------------------------
    // (3) folding matrix,
    // and diagonalize the H matrix (T+Vl+Vnl).
    //--------------------------------------------

    if(XC_Functional::get_func_type()==3 || XC_Functional::get_func_type()==5)
    {
        Gint_inout inout(vr_eff1, vofk_eff1, Gint_Tools::job_type::vlocal_meta);
        this->GG->cal_vlocal(&inout, this->LM, this->new_e_iteration);
    }
    else
    {
        Gint_inout inout(vr_eff1, Gint_Tools::job_type::vlocal);
        this->GG->cal_vlocal(&inout, this->LM, this->new_e_iteration);
    }
    this->GG->transfer_pvpR(this->hR);

    this->new_e_iteration = false;
}

}