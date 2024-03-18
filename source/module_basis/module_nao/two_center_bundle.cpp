#include "module_basis/module_nao/two_center_bundle.h"
#include "module_base/memory.h"
#include "module_base/ylm.h"
#include "module_basis/module_nao/real_gaunt_table.h"
#include <memory>
#include "module_base/parallel_common.h"
#include "module_base/global_variable.h"

void TwoCenterBundle::build_orb(int ntype, const std::string* file_orb0)
{
    std::vector<std::string> file_orb(ntype);
    if (GlobalV::MY_RANK == 0)
    {
        std::transform(file_orb0, file_orb0 + ntype, file_orb.begin(),
                [](const std::string& file) { return GlobalV::global_orbital_dir + file; });
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
    if (GlobalV::deepks_setorb)
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

void TwoCenterBundle::build_orb_onsite(int ntype, double radius)
{
    if(GlobalV::onsite_radius > 0)
    {
        orb_onsite_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        orb_onsite_->build(orb_.get(), GlobalV::onsite_radius);
    }
}

void TwoCenterBundle::tabulate()
{
    ModuleBase::SphericalBesselTransformer sbt;
    orb_->set_transformer(sbt);
    beta_->set_transformer(sbt);
    if (alpha_) alpha_->set_transformer(sbt);
    if (orb_onsite_) orb_onsite_->set_transformer(sbt);

    //================================================================
    //              build two-center integration tables
    //================================================================
    // set up a universal radial grid
    double rmax = std::max(orb_->rcut_max(), beta_->rcut_max());
    if (alpha_) rmax = std::max(rmax, alpha_->rcut_max());
    double dr = 0.01;
    double cutoff = 2.0 * rmax;
    int nr = static_cast<int>(rmax / dr) + 1;

    orb_->set_uniform_grid(true, nr, cutoff, 'i', true);
    beta_->set_uniform_grid(true, nr, cutoff, 'i', true);
    if (alpha_) alpha_->set_uniform_grid(true, nr, cutoff, 'i', true);
    if (orb_onsite_) orb_onsite_->set_uniform_grid(true, nr, cutoff, 'i', true);

    // build TwoCenterIntegrator objects
    kinetic_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    kinetic_orb->tabulate(*orb_, *orb_, 'T', nr, cutoff);
    ModuleBase::Memory::record("TwoCenterTable: Kinetic", kinetic_orb->table_memory());

    overlap_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb->tabulate(*orb_, *orb_, 'S', nr, cutoff);
    ModuleBase::Memory::record("TwoCenterTable: Overlap", overlap_orb->table_memory());

    overlap_orb_beta = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb_beta->tabulate(*orb_, *beta_, 'S', nr, cutoff);
    ModuleBase::Memory::record("TwoCenterTable: Nonlocal", overlap_orb_beta->table_memory());

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

    // init Ylm (this shall be done by Ylm automatically! to be done later...)
    ModuleBase::Ylm::set_coefficients();

    sbt.fft_clear();
}

void TwoCenterBundle::to_LCAO_Orbitals(LCAO_Orbitals& ORB, const double lcao_ecut, const double lcao_dk, const double lcao_dr, const double lcao_rmax) const
{
    ORB.ntype = orb_->ntype();
    ORB.lmax = orb_->lmax();
    ORB.nchimax = orb_->nzeta_max();
    ORB.rcutmax_Phi = orb_->rcut_max();
    ORB.dR = lcao_dr;
    ORB.Rmax = lcao_rmax; // lcao_rmax, see ORB_control.cpp
    ORB.dr_uniform = 0.001;
    
    // Due to algorithmic difference in the spherical Bessel transform
    // (FFT vs. Simpson's integration), k grid of FFT is not appropriate
    // for Simpson's integration. The k grid for Simpson's integration is
    // specifically set by the two variables below to mimick the original
    // behavior.
    ORB.ecutwfc = lcao_ecut;
    ORB.dk = lcao_dk;

	if(ORB.ecutwfc < 20)
	{
		ORB.kmesh = static_cast<int>( 2 * sqrt(ORB.ecutwfc) / ORB.dk )  + 4;
	}
	else
	{
		ORB.kmesh = static_cast<int>( sqrt(ORB.ecutwfc) / ORB.dk )  + 4;
	}
    ORB.kmesh += 1 - ORB.kmesh % 2;


    delete[] ORB.Phi;
    ORB.Phi = new Numerical_Orbital[orb_->ntype()];
    for (int itype = 0; itype < orb_->ntype(); ++itype)
    {
        (*orb_)(itype).to_numerical_orbital(ORB.Phi[itype], ORB.kmesh, ORB.dk);
    }

    if (GlobalV::deepks_setorb)
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
