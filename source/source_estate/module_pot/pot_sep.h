#ifndef POTSEP_H
#define POTSEP_H

#include "pot_base.h"
#include "source_base/matrix.h"
#include "source_pw/module_pwdft/VSep_in_pw.h"

namespace elecstate
{

class PotSep : public PotBase
{
  public:
    // PotSep(const ModuleBase::matrix* vsep_in,
    //        const ModuleBase::ComplexMatrix* sf_in,
    //        const ModulePW::PW_Basis* rho_basis_in,
    //        const bool* sep_enable_in)
    //     : vsep_(vsep_in), sf_(sf_in), sep_enable_(sep_enable_in)
    // {
    //     assert(this->vsep_->nr == this->sf_->nr);
    //     this->rho_basis_ = rho_basis_in;
    //     this->ntype_ = this->vsep_->nr;
    //     this->fixed_mode = true;
    //     this->dynamic_mode = false;
    // }
    PotSep(const ModuleBase::ComplexMatrix* sf_in, const ModulePW::PW_Basis* rho_basis_in, const VSep* vsep_cell_in)
        : sf_(sf_in), vsep_cell(vsep_cell_in)
    {
        assert(vsep_cell->vsep_form.nr == this->sf_->nr);
        // assert(this->vsep_->vsep_form.nr == this->sf_->nr);
        this->rho_basis_ = rho_basis_in;
        // this->ntype_ = this->vsep_->vsep_form.nr;
        this->fixed_mode = true;
        this->dynamic_mode = false;
    }

    void cal_fixed_v(double* vl_pseudo) override;

    const VSep* vsep_cell = nullptr;
    const ModuleBase::ComplexMatrix* sf_ = nullptr;
    // int ntype_ = 0;
    // const bool* sep_enable_;
};

} // namespace elecstate

#endif /* ifndef POTSEP_H */
