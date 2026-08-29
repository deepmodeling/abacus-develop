#include "veff_lcao.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/tool_title.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include <cassert>
namespace hamilt
{

namespace
{
HContainer<double> make_real_hcontainer(const HContainer<std::complex<double>>& complex_hr)
{
    HContainer<double> real_hr(complex_hr.get_paraV());
    for (int iap = 0; iap < complex_hr.size_atom_pairs(); ++iap)
    {
        const AtomPair<std::complex<double>>& complex_ap = complex_hr.get_atom_pair(iap);
        for (int ir = 0; ir < complex_ap.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> r_index = complex_ap.get_R_index(ir);
            real_hr.insert_pair(AtomPair<double>(complex_ap.get_atom_i(), complex_ap.get_atom_j(), r_index,
                                                 complex_hr.get_paraV()));
        }
    }
    real_hr.allocate(nullptr, true);
    return real_hr;
}

void add_real_hcontainer(const HContainer<double>& real_hr, HContainer<std::complex<double>>& complex_hr)
{
    for (int iap = 0; iap < real_hr.size_atom_pairs(); ++iap)
    {
        const AtomPair<double>& real_ap = real_hr.get_atom_pair(iap);
        AtomPair<std::complex<double>>* complex_ap
            = complex_hr.find_pair(real_ap.get_atom_i(), real_ap.get_atom_j());
        assert(complex_ap != nullptr);
        for (int ir = 0; ir < real_ap.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> r_index = real_ap.get_R_index(ir);
            const BaseMatrix<double>& real_matrix = real_ap.get_HR_values(ir);
            BaseMatrix<std::complex<double>>* complex_matrix = complex_ap->find_matrix(r_index);
            assert(complex_matrix != nullptr);
            for (int index = 0; index < real_ap.get_size(); ++index)
            {
                complex_matrix->get_pointer()[index] += real_matrix.get_pointer()[index];
            }
        }
    }
}
} // namespace


// initialize_HR()
template <typename TK, typename TR>
void Veff<OperatorLCAO<TK, TR>>::initialize_HR(const UnitCell* ucell_in, const Grid_Driver* GridD)
{
    ModuleBase::TITLE("Veff", "initialize_HR");
    ModuleBase::timer::start("Veff", "initialize_HR");

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
            if (paraV->is_invalid_atom_pair(iat1, iat2))
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

    ModuleBase::timer::end("Veff", "initialize_HR");
}

template<>
void Veff<OperatorLCAO<double, double>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::start("Veff", "contributeHR");
    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    double* vr_eff1 = this->pot->get_eff_v(this->current_spin);
    double* vofk_eff1 = this->pot->get_eff_vofk(this->current_spin);

    if(XC_Functional::get_ked_flag())
    {
        ModuleGint::cal_gint_vl_metagga(vr_eff1, vofk_eff1, this->hR);
    }
    else
    {
        ModuleGint::cal_gint_vl(vr_eff1, this->hR);
    }

    if(this->nspin == 2) 
    { 
        this->current_spin = 1 - this->current_spin;
    }

    ModuleBase::timer::end("Veff", "contributeHR");
    return;
}

template<>
void Veff<OperatorLCAO<std::complex<double>, double>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::start("Veff", "contributeHR");
    //-----------------------------------------
    //(1) prepare data for this k point.
    // copy the local potential from array.
    //-----------------------------------------
    double* vr_eff1 = this->pot->get_eff_v(this->current_spin);
    double* vofk_eff1 = this->pot->get_eff_vofk(this->current_spin);

    if(XC_Functional::get_ked_flag())
    {
        ModuleGint::cal_gint_vl_metagga(vr_eff1, vofk_eff1, this->hR);
    }
    else
    {
        ModuleGint::cal_gint_vl(vr_eff1, this->hR);
    }

    if(this->nspin == 2) 
    { 
        this->current_spin = 1 - this->current_spin;
    }

    ModuleBase::timer::end("Veff", "contributeHR");
    return;
}

template<>
void Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::start("Veff", "contributeHR");

    if (this->nspin < 4)
    {
        HContainer<double> real_hr = make_real_hcontainer(*this->hR);
        double* vr_eff = this->pot->get_eff_v(this->current_spin);
        double* vofk_eff = this->pot->get_eff_vofk(this->current_spin);
        if (XC_Functional::get_ked_flag())
        {
            ModuleGint::cal_gint_vl_metagga(vr_eff, vofk_eff, &real_hr);
        }
        else
        {
            ModuleGint::cal_gint_vl(vr_eff, &real_hr);
        }
        add_real_hcontainer(real_hr, *this->hR);
        if (this->nspin == 2)
        {
            this->current_spin = 1 - this->current_spin;
        }
    }
    else
    {
        std::vector<const double*> vr_eff(4, nullptr);
        std::vector<const double*> vofk_eff(4, nullptr);
        for (int is = 0; is < 4; is++)
        {
            vr_eff[is] = this->pot->get_eff_v(is);
            if(XC_Functional::get_ked_flag())
            {
                vofk_eff[is] = this->pot->get_eff_vofk(is);
            }
        }
        if(XC_Functional::get_ked_flag())
        {
            ModuleGint::cal_gint_vl_metagga(vr_eff, vofk_eff, this->hR);
        }
        else
        {
            ModuleGint::cal_gint_vl(vr_eff, this->hR);
        }
    }

    ModuleBase::timer::end("Veff", "contributeHR");
    return;
}

// definition of class template should in the end of file to avoid compiling warning 
template class Veff<OperatorLCAO<double, double>>;

template class Veff<OperatorLCAO<std::complex<double>, double>>;

template class Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>;
}
