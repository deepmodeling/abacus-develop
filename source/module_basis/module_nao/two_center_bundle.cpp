#include "module_basis/module_nao/two_center_bundle.h"
#include "module_base/ylm.h"
#include "module_basis/module_nao/real_gaunt_table.h"
#include <memory>

TwoCenterBundle::~TwoCenterBundle()
{
}

void TwoCenterBundle::build(const int nfile_orb,
                            const std::string* file_orb,
                            const int nfile_pp,
                            const std::string* file_pp,
                            const int nfile_desc,
                            const std::string* file_desc)
{
    // build RadialCollection objects
    orb_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    orb_->build(nfile_orb, file_orb, 'o');

    beta_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    beta_->build(nfile_pp, file_pp, 'p');

    double rmax = std::max(orb_->rcut_max(), beta_->rcut_max());

    if (nfile_desc > 0)
    {
        alpha_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        alpha_->build(nfile_desc, file_desc, 'o');
        rmax = std::max(rmax, alpha_->rcut_max());
    }

    // set up a universal radial grid
    double dr = 0.01;
    double cutoff = 2.0 * rmax;
    int nr = static_cast<int>(rmax / dr) + 1;
    orb_->set_uniform_grid(true, nr, cutoff, 'i', true);
    beta_->set_uniform_grid(true, nr, cutoff, 'i', true);
    if (nfile_desc > 0) alpha_->set_uniform_grid(true, nr, cutoff, 'i', true);

    // build TwoCenterIntegrator objects
    kinetic_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    kinetic_orb->tabulate(*orb_, *orb_, 'T', nr, cutoff, true);

    overlap_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb->tabulate(*orb_, *orb_, 'S', nr, cutoff, true);

    overlap_orb_beta = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb_beta->tabulate(*orb_, *beta_, 'S', nr, cutoff, true);

    if (nfile_desc > 0)
    {
        overlap_orb_alpha = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        overlap_orb_alpha->tabulate(*orb_, *alpha_, 'S', nr, cutoff, true);
    }

    // init Ylm (this shall be done by Ylm automatically! to be done later...)
    ModuleBase::Ylm::set_coefficients();
}
