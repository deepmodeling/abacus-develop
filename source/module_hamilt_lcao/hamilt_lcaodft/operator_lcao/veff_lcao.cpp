#include "veff_lcao.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "module_base/tool_title.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_gint/temp_gint/gint_interface.h"
namespace hamilt
{


// initialize_HR()
template <typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::initialize_HR(const UnitCell* ucell_in, const Grid_Driver* GridD)
{
    ModuleBase::TITLE("Veff", "initialize_HR");
    ModuleBase::timer::tick("Veff", "initialize_HR");

    this->nspin = PARAM.inp.nspin;
    auto* paraV = this->hR->get_paraV();// get parallel orbitals from HR
    // TODO: if paraV is nullptr, AtomPair can not use paraV for constructor, I will repair it in the future.

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
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell_in->cal_dtau(iat1, iat2, R_index2).norm() * ucell_in->lat0
                < orb_cutoff_[T1] + orb_cutoff_[T2])
            {
                hamilt::AtomPair<TR> tmp(iat1, iat2, R_index2, paraV);
                this->hR->insert_pair(tmp);
            }
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    this->hR->allocate(nullptr, true);

    ModuleBase::timer::tick("Veff", "initialize_HR");
}

template<>
void Veff<OperatorLCAO<double, double>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::tick("Veff", "contributeHR");
    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    double* vr_eff1 = this->pot->get_effective_v(this->current_spin);
    double* vofk_eff1 = this->pot->get_effective_vofk(this->current_spin);

#ifndef __NEW_GINT
    if(XC_Functional::get_ked_flag())
    {
        Gint_inout inout(vr_eff1, vofk_eff1, Gint_Tools::job_type::vlocal_meta);
        this->GG->cal_vlocal(&inout,  this->new_e_iteration);
    }
    else
    {
        Gint_inout inout(vr_eff1, Gint_Tools::job_type::vlocal);
        this->GG->cal_vlocal(&inout,  this->new_e_iteration);
    }
    this->GG->transfer_pvpR(this->hR,this->ucell);
    this->new_e_iteration = false;
#else
    if(XC_Functional::get_ked_flag())
    {
        ModuleGint::cal_gint_vl_metagga(vr_eff1, vofk_eff1, this->hR);
    }
    else
    {
        ModuleGint::cal_gint_vl(vr_eff1, this->hR);
    }
#endif

    if(this->nspin == 2) 
    { 
        this->current_spin = 1 - this->current_spin;
    }

    ModuleBase::timer::tick("Veff", "contributeHR");
    return;
}

template<>
void Veff<OperatorLCAO<std::complex<double>, double>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::tick("Veff", "contributeHR");
    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    double* vr_eff1 = this->pot->get_effective_v(this->current_spin);
    double* vofk_eff1 = this->pot->get_effective_vofk(this->current_spin);

#ifndef __NEW_GINT
    // if you change the place of the following code,
    // rememeber to delete the #include
    if(XC_Functional::get_ked_flag())
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

    this->GK->transfer_pvpR(this->hR,this->ucell,this->gd);
#else
    if(XC_Functional::get_ked_flag())
    {
        ModuleGint::cal_gint_vl_metagga(vr_eff1, vofk_eff1, this->hR);
    }
    else
    {
        ModuleGint::cal_gint_vl(vr_eff1, this->hR);
    }
#endif

    if(this->nspin == 2) 
    { 
        this->current_spin = 1 - this->current_spin;
    }

    ModuleBase::timer::tick("Veff", "contributeHR");
    return;
}

template<>
void Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::tick("Veff", "contributeHR");

#ifndef __NEW_GINT
    double* vr_eff1 = nullptr;
    double* vofk_eff1 = nullptr;
    for (int is = 0; is < 4; is++)
    {
        vr_eff1 = this->pot->get_effective_v(is);
        if(XC_Functional::get_ked_flag())
        {
            vofk_eff1 = this->pot->get_effective_vofk(is);
        }
        
        if(XC_Functional::get_ked_flag())
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
    this->GK->transfer_pvpR(this->hR,this->ucell,this->gd);
#else
    std::vector<const double*> vr_eff(4, nullptr);
    std::vector<const double*> vofk_eff(4, nullptr);
    for (int is = 0; is < 4; is++)
    {
        vr_eff[is] = this->pot->get_effective_v(is);
        if(XC_Functional::get_ked_flag())
        {
            vofk_eff[is] = this->pot->get_effective_vofk(is);
            if(is == 3)
            {
                ModuleGint::cal_gint_vl_metagga(vr_eff, vofk_eff, this->hR);
            }
        }
        else
        {
            if(is == 3)
            {
                ModuleGint::cal_gint_vl(vr_eff, this->hR);
            }
        }
    }
#endif

    ModuleBase::timer::tick("Veff", "contributeHR");
    return;
}

// definition of class template should in the end of file to avoid compiling warning 
template class Veff<OperatorLCAO<double, double>>;

template class Veff<OperatorLCAO<std::complex<double>, double>>;

template class Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>;
}