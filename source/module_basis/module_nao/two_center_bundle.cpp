#include "module_basis/module_nao/two_center_bundle.h"

#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/parallel_common.h"
#include "module_base/ylm.h"
#include "module_basis/module_nao/real_gaunt_table.h"
#include "module_parameter/parameter.h"

#include <memory>

void TwoCenterBundle::build_orb(int ntype, const std::string* file_orb0)
{
    std::vector<std::string> file_orb(ntype);
    if (GlobalV::MY_RANK == 0)
    {
        std::transform(file_orb0, file_orb0 + ntype, file_orb.begin(), [](const std::string& file) {
            return PARAM.inp.orbital_dir + file;
        });
    }
#ifdef __MPI
    Parallel_Common::bcast_string(file_orb.data(), ntype);
#endif

    orb_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    orb_->build(ntype, file_orb.data()); // automatically detect file type
}

void TwoCenterBundle::build_beta(int ntype, Numerical_Nonlocal* nl)
{
    beta_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    beta_->build(ntype, nl);
}

void TwoCenterBundle::build_alpha(int ndesc, std::string* file_desc0)
{
    if (PARAM.globalv.deepks_setorb)
    {
        std::vector<std::string> file_desc(ndesc);
        if (GlobalV::MY_RANK == 0)
        {
            std::copy(file_desc0, file_desc0 + ndesc, file_desc.begin());
        }
#ifdef __MPI
        Parallel_Common::bcast_string(file_desc.data(), ndesc);
#endif

        alpha_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        alpha_->build(ndesc, file_desc.data());
    }
}

void TwoCenterBundle::build_orb_onsite(const double& radius)
{
    if (radius > 0)
    {
        orb_onsite_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        orb_onsite_->build(orb_.get(), radius);
    }
}

void TwoCenterBundle::tabulate()
{
    ModuleBase::SphericalBesselTransformer sbt(true);
    orb_->set_transformer(sbt);
    if (beta_) { beta_->set_transformer(sbt); }
    if (alpha_) {
        alpha_->set_transformer(sbt);
}
    if (orb_onsite_) {
        orb_onsite_->set_transformer(sbt);
}

    //================================================================
    //              build two-center integration tables
    //================================================================
    // set up a universal radial grid
    double rmax = orb_->rcut_max();
    if (beta_) { rmax = std::max(rmax, beta_->rcut_max()); }
    if (alpha_) { rmax = std::max(rmax, alpha_->rcut_max()); }
    double dr = 0.01;
    double cutoff = 2.0 * rmax;
    int nr = static_cast<int>(rmax / dr) + 1;

    orb_->set_uniform_grid(true, nr, cutoff, 'i', true);
    if (beta_) { beta_->set_uniform_grid(true, nr, cutoff, 'i', true); }
    if (alpha_) { alpha_->set_uniform_grid(true, nr, cutoff, 'i', true);}
    if (orb_onsite_) {  orb_onsite_->set_uniform_grid(true, nr, cutoff, 'i', true);}

    // build TwoCenterIntegrator objects
    kinetic_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    kinetic_orb->tabulate(*orb_, *orb_, 'T', nr, cutoff);
    ModuleBase::Memory::record("TwoCenterTable: Kinetic", kinetic_orb->table_memory());

    overlap_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb->tabulate(*orb_, *orb_, 'S', nr, cutoff);
    ModuleBase::Memory::record("TwoCenterTable: Overlap", overlap_orb->table_memory());

    if (beta_)
    {
        overlap_orb_beta = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        overlap_orb_beta->tabulate(*orb_, *beta_, 'S', nr, cutoff);
        ModuleBase::Memory::record("TwoCenterTable: Nonlocal", overlap_orb_beta->table_memory());
    }

    if (alpha_)
    {
        overlap_orb_alpha = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        overlap_orb_alpha->tabulate(*orb_, *alpha_, 'S', nr, cutoff);
        ModuleBase::Memory::record("TwoCenterTable: Descriptor", overlap_orb_beta->table_memory());
    }

    if (orb_onsite_)
    {
        overlap_orb_onsite = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        overlap_orb_onsite->tabulate(*orb_, *orb_onsite_, 'S', nr, cutoff);
    }

    ModuleBase::Memory::record("RealGauntTable", RealGauntTable::instance().memory());

    sbt.clear();
}

void TwoCenterBundle::tabulate(const double lcao_ecut,
                               const double lcao_dk,
                               const double lcao_dr,
                               const double lcao_rmax)
{
    ModuleBase::SphericalBesselTransformer sbt(true);
    orb_->set_transformer(sbt);
    beta_->set_transformer(sbt);
    if (alpha_) {
        alpha_->set_transformer(sbt);
}
    if (orb_onsite_) {
        orb_onsite_->set_transformer(sbt);
}

    //================================================================
    //              build two-center integration tables
    //================================================================

    // old formula for the number of k-space grid points
    int nk = static_cast<int>(sqrt(lcao_ecut) / lcao_dk) + 4;
    nk += 1 - nk % 2; // make nk odd

    std::vector<double> kgrid(nk);
    for (int ik = 0; ik < nk; ++ik)
    {
        kgrid[ik] = ik * lcao_dk;
    }

    orb_->set_grid(false, nk, kgrid.data(), 't');
    beta_->set_grid(false, nk, kgrid.data(), 't');
    if (alpha_)
    {
        alpha_->set_grid(false, nk, kgrid.data(), 't');
    }
    if (orb_onsite_)
    {
        orb_onsite_->set_grid(false, nk, kgrid.data(), 't');
    }

    // "st" stands for overlap (s) and kinetic (t)
    const double cutoff_st = std::min(lcao_rmax, 2.0 * orb_->rcut_max());
    const int nr_st = static_cast<int>(cutoff_st / lcao_dr) + 5;

    kinetic_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    kinetic_orb->tabulate(*orb_, *orb_, 'T', nr_st, cutoff_st);
    ModuleBase::Memory::record("TwoCenterTable: Kinetic", kinetic_orb->table_memory());

    overlap_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb->tabulate(*orb_, *orb_, 'S', nr_st, cutoff_st);
    ModuleBase::Memory::record("TwoCenterTable: Overlap", overlap_orb->table_memory());

    // overlap between orbital and beta (for nonlocal potential)
    const double cutoff_nl = std::min(lcao_rmax, orb_->rcut_max() + beta_->rcut_max());
    const int nr_nl = static_cast<int>(cutoff_nl / lcao_dr) + 5;
    overlap_orb_beta = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb_beta->tabulate(*orb_, *beta_, 'S', nr_nl, cutoff_nl);
    ModuleBase::Memory::record("TwoCenterTable: Nonlocal", overlap_orb_beta->table_memory());

    // overlap between orbital and deepks projector
    if (alpha_)
    {
        const double cutoff_alpha = std::min(lcao_rmax, orb_->rcut_max() + alpha_->rcut_max());
        const int nr_alpha = static_cast<int>(cutoff_alpha / lcao_dr) + 5;
        overlap_orb_alpha = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        overlap_orb_alpha->tabulate(*orb_, *alpha_, 'S', nr_alpha, cutoff_alpha);
        ModuleBase::Memory::record("TwoCenterTable: Descriptor", overlap_orb_beta->table_memory());
    }

    // overlap between orbital and "onsite orbital" (for DFT+U)
    if (orb_onsite_)
    {
        const double cutoff_onsite = std::min(lcao_rmax, orb_->rcut_max() + orb_onsite_->rcut_max());
        const int nr_onsite = static_cast<int>(cutoff_onsite / lcao_dr) + 5;
        overlap_orb_onsite = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        overlap_orb_onsite->tabulate(*orb_, *orb_onsite_, 'S', nr_onsite, cutoff_onsite);
    }

    ModuleBase::Memory::record("RealGauntTable", RealGauntTable::instance().memory());

    sbt.clear();
}

void TwoCenterBundle::to_LCAO_Orbitals(LCAO_Orbitals& ORB,
                                       const double lcao_ecut,
                                       const double lcao_dk,
                                       const double lcao_dr,
                                       const double lcao_rmax) const
{
    ORB.ntype = orb_->ntype();
    ORB.lmax = orb_->lmax();
    ORB.nchimax = orb_->nzeta_max();
    ORB.rcutmax_Phi = orb_->rcut_max();
    ORB.dR = lcao_dr;
    ORB.Rmax = lcao_rmax;
    ORB.dr_uniform = 0.001;

    // Due to algorithmic difference in the spherical Bessel transform
    // (FFT vs. Simpson's integration), k grid of FFT is not appropriate
    // for Simpson's integration. The k grid for Simpson's integration is
    // specifically set by the two variables below to mimick the original
    // behavior.
    ORB.ecutwfc = lcao_ecut;
    ORB.dk = lcao_dk;

    if (ORB.ecutwfc < 20)
    {
        ORB.kmesh = static_cast<int>(2 * sqrt(ORB.ecutwfc) / ORB.dk) + 4;
    }
    else
    {
        ORB.kmesh = static_cast<int>(sqrt(ORB.ecutwfc) / ORB.dk) + 4;
    }
    ORB.kmesh += 1 - ORB.kmesh % 2;

    delete[] ORB.Phi;
    ORB.Phi = new Numerical_Orbital[orb_->ntype()];
    for (int itype = 0; itype < orb_->ntype(); ++itype)
    {
        (*orb_)(itype).to_numerical_orbital(ORB.Phi[itype], ORB.kmesh, ORB.dk);
    }

    if (PARAM.globalv.deepks_setorb)
    {
        ORB.lmax_d = alpha_->lmax();
        ORB.nchimax_d = alpha_->nzeta_max();

        delete[] ORB.Alpha;
        ORB.Alpha = new Numerical_Orbital[alpha_->ntype()];
        for (int itype = 0; itype < alpha_->ntype(); ++itype)
        {
            (*alpha_)(itype).to_numerical_orbital(ORB.Alpha[itype], ORB.kmesh, ORB.dk);
        }
    }
}
