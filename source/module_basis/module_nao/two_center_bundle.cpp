#include "module_basis/module_nao/two_center_bundle.h"
#include "module_base/ylm.h"
#include "module_basis/module_nao/real_gaunt_table.h"
#include <memory>
#include "module_base/parallel_common.h"
#include "module_base/global_variable.h"

TwoCenterBundle::~TwoCenterBundle()
{
}

void TwoCenterBundle::build(const int nfile_orb,
                            const std::string* file_orb0,
                            const int nfile_pp,
                            const std::string* file_pp0,
                            const int nfile_desc,
                            const std::string* file_desc0)
{
    //================================================================
    //                      read in the files
    //================================================================
    // NOTE: only RANK-0 has the complete file name information; a broadcast is necessary
    // NOTE: the passed-in file names do not contain the directory information

    // currently the number of orbital files and pseudo files are the same
    int ntype = nfile_orb; 
#ifdef __MPI
    Parallel_Common::bcast_int(ntype);
#endif

    std::string* file_orb = new std::string[ntype];
    std::string* file_pp = new std::string[ntype];
    if (GlobalV::MY_RANK == 0) {
        for (int it = 0; it < ntype; ++it)
        {
            file_orb[it] = GlobalV::global_orbital_dir + file_orb0[it];
            file_pp[it] = GlobalV::global_pseudo_dir + file_pp0[it];;
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_string(file_orb, ntype);
    Parallel_Common::bcast_string(file_pp, ntype);
#endif

    // build RadialCollection objects
    orb_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    orb_->build(nfile_orb, file_orb, 'o');

    beta_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    beta_->build(nfile_pp, file_pp, 'p');

    double rmax = std::max(orb_->rcut_max(), beta_->rcut_max());

    delete[] file_orb;
    delete[] file_pp;

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
        alpha_->build(nfile_desc, file_desc, 'o');
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
    kinetic_orb->tabulate(*orb_, *orb_, 'T', nr, cutoff, true);

    overlap_orb = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb->tabulate(*orb_, *orb_, 'S', nr, cutoff, true);

    overlap_orb_beta = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    overlap_orb_beta->tabulate(*orb_, *beta_, 'S', nr, cutoff, true);

    if (deepks_on)
    {
        overlap_orb_alpha = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        overlap_orb_alpha->tabulate(*orb_, *alpha_, 'S', nr, cutoff, true);
    }

    // init Ylm (this shall be done by Ylm automatically! to be done later...)
    ModuleBase::Ylm::set_coefficients();
}
