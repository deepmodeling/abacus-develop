#include "module_basis/module_nao/two_center_bundle.h"
#include "module_base/memory.h"
#include "module_base/ylm.h"
#include "module_basis/module_nao/real_gaunt_table.h"
#include <memory>
#include "module_base/parallel_common.h"
#include "module_base/global_variable.h"

void TwoCenterBundle::build(int ntype,
                            const std::string* file_orb0,
                            Numerical_Nonlocal* nl,
                            const int nfile_desc,
                            const std::string* file_desc0)
{

    ModuleBase::SphericalBesselTransformer sbt;

    //================================================================
    //                      read in the files
    //================================================================
    // NOTE: only RANK-0 has the complete file name information; a broadcast is necessary
    // NOTE: the passed-in file names do not contain the directory information

#ifdef __MPI
    Parallel_Common::bcast_int(ntype);
#endif

    std::string* file_orb = new std::string[ntype];
    if (GlobalV::MY_RANK == 0) {
        for (int it = 0; it < ntype; ++it)
        {
            file_orb[it] = GlobalV::global_orbital_dir + file_orb0[it];
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_string(file_orb, ntype);
#endif

    // build RadialCollection objects
    orb_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    orb_->build(ntype, file_orb);
    orb_->set_transformer(sbt);

    beta_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    beta_->build(ntype, nl);
    beta_->set_transformer(sbt);

    double rmax = std::max(orb_->rcut_max(), beta_->rcut_max());

    delete[] file_orb;

    //========== DeePKS =========
    bool deepks_on = GlobalV::deepks_setorb;
#ifdef __MPI
    Parallel_Common::bcast_bool(deepks_on);
#endif

    if (deepks_on)
    {
        int ndesc = nfile_desc;
#ifdef __MPI
        Parallel_Common::bcast_int(ndesc);
#endif

        std::string* file_desc = new std::string[ndesc];
        if (GlobalV::MY_RANK == 0) {
            for (int it = 0; it < ndesc; ++it)
            {
                file_desc[it] = file_desc0[it];;
            }
        }
#ifdef __MPI
        Parallel_Common::bcast_string(file_desc, ndesc);
#endif

        alpha_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        alpha_->build(nfile_desc, file_desc);
        alpha_->set_transformer(sbt);
        rmax = std::max(rmax, alpha_->rcut_max());

        delete[] file_desc;
    }

    //================================================================
    //              build two-center integration tables
    //================================================================
    // set up a universal radial grid
    double dr = 0.01;
    double cutoff = 2.0 * rmax;
    int nr = static_cast<int>(rmax / dr) + 1;
    orb_->set_uniform_grid(true, nr, cutoff, 'i', true);
    beta_->set_uniform_grid(true, nr, cutoff, 'i', true);
    if (deepks_on) alpha_->set_uniform_grid(true, nr, cutoff, 'i', true);

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

    if (deepks_on)
    {
        overlap_orb_alpha = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        overlap_orb_alpha->tabulate(*orb_, *alpha_, 'S', nr, cutoff);
        ModuleBase::Memory::record("TwoCenterTable: Descriptor", overlap_orb_beta->table_memory());
    }

    ModuleBase::Memory::record("RealGauntTable", RealGauntTable::instance().memory());

    // init Ylm (this shall be done by Ylm automatically! to be done later...)
    ModuleBase::Ylm::set_coefficients();

    sbt.fft_clear();
}

void TwoCenterBundle::to_LCAO_Orbitals(LCAO_Orbitals& ORB) const
{
    ORB.ntype = orb_->ntype();
    ORB.kmesh = (*orb_)(0, 0, 0).nk();
    ORB.lmax = orb_->lmax();
    ORB.nchimax = orb_->nzeta_max();
    ORB.rcutmax_Phi = orb_->rcut_max();
    
    delete[] ORB.Phi;
    ORB.Phi = new Numerical_Orbital[orb_->ntype()];
    for (int itype = 0; itype < orb_->ntype(); ++itype)
    {
        (*orb_)(itype).to_numerical_orbital(ORB.Phi[itype]);
    }

    if (GlobalV::deepks_setorb)
    {
        ORB.lmax_d = alpha_->lmax();
        ORB.nchimax_d = alpha_->nzeta_max();

        delete[] ORB.Alpha;
        ORB.Alpha = new Numerical_Orbital[alpha_->ntype()];
        for (int itype = 0; itype < alpha_->ntype(); ++itype)
        {
            (*alpha_)(itype).to_numerical_orbital(ORB.Alpha[itype]);
        }
    }
}
