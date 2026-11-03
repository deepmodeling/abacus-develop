#include "pot_surchem.h"
#include "sol_grid_field.h"

namespace elecstate
{

PotSurChem::PotSurChem(const ModulePW::PW_Basis* rho_basis_in,
                       Structure_Factor* structure_factors_in,
                       const double* vlocal_in,
                       surchem* surchem_in)
    : vlocal(vlocal_in), surchem_(surchem_in)
{
    this->rho_basis_ = rho_basis_in;
    this->structure_factors_ = structure_factors_in;
    this->dynamic_mode = true;
    this->fixed_mode = false;
}

PotSurChem::~PotSurChem()
{
    if (this->allocated)
    {
        this->surchem_->clear();
    }
}

void PotSurChem::cal_v_eff(const Charge* const chg, const UnitCell* const ucell, ModuleBase::matrix& v_eff)
{
    // nspin from the charge: v_eff has no rows on a rank without real-space grid points.
    const int nspin = chg->nspin;
    if (!this->allocated)
    {
        this->surchem_->allocate(this->rho_basis_->nrxx, nspin);
        this->allocated = true;
    }
    ModuleBase::matrix v_sol_correction(nspin, this->rho_basis_->nrxx);
    this->surchem_->v_correction(*ucell,
                                 *chg->pgrid,
                                 const_cast<ModulePW::PW_Basis*>(this->rho_basis_),
                                 nspin,
                                 chg->rho,
                                 this->vlocal,
                                 this->structure_factors_,
                                 v_sol_correction);
    v_eff += v_sol_correction;
}

void PotSurChem::add_sol_fields(std::vector<SolventGridField>& fields) const
{
    // Fields are written collectively; last_epsilon() is empty on a rank without
    // grid points, has_last_eps() is set on every rank.
    if (!this->surchem_->has_last_eps())
    {
        return;
    }
    const std::vector<double>& epsilon = this->surchem_->last_epsilon();
    SolventGridField dielectric;
    dielectric.name = "eps";
    dielectric.values = epsilon;
    fields.push_back(dielectric);
    // epsilon = 1 + (eb_k - 1) S with the solvent shape function S, so the
    // cavity s = 1 - S follows from epsilon; without a dielectric it is undefined.
    const double contrast = this->surchem_->bulk_epsilon() - 1.0;
    if (contrast <= 0.0)
    {
        return;
    }
    SolventGridField cavity;
    cavity.name = "cavity";
    cavity.values.resize(epsilon.size());
    for (std::size_t ir = 0; ir < epsilon.size(); ++ir)
    {
        const double shape = (epsilon[ir] - 1.0) / contrast;
        cavity.values[ir] = 1.0 - shape;
    }
    fields.push_back(cavity);
}

} // namespace elecstate
