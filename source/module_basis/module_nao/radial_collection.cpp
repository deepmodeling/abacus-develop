#include "module_basis/module_nao/radial_collection.h"
#include <memory>

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/atomic_radials.h"
#include "module_basis/module_nao/beta_radials.h"
#include "module_basis/module_nao/hydrogen_radials.h"

RadialCollection::RadialCollection(const RadialCollection& other) :
    ntype_(other.ntype_),
    lmax_(other.lmax_),
    nchi_(other.nchi_),
    nzeta_max_(other.nzeta_max_),
    rcut_max_(other.rcut_max_),
    radset_(nullptr),
    iter_(nullptr),
    nl_(nullptr)
{
    if (ntype_ == 0)
    {
        return;
    }

    radset_ = new RadialSet*[ntype_];
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = other.radset_[itype]->clone();
    }

    iter_build();
}

RadialCollection& RadialCollection::operator=(const RadialCollection& rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    cleanup();

    ntype_ = rhs.ntype_;
    lmax_ = rhs.lmax_;
    nchi_ = rhs.nchi_;
    nzeta_max_ = rhs.nzeta_max_;
    rcut_max_ = rhs.rcut_max_;

    radset_ = new RadialSet*[ntype_];
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = rhs.radset_[itype]->clone();
    }

    iter_build();

    return *this;
}

RadialCollection::~RadialCollection()
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        delete radset_[itype];
    }
    delete[] radset_;
    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    delete[] nl_;
}

void RadialCollection::set_rcut_max()
{
    rcut_max_ = 0.0;
    for (int itype = 0; itype < ntype_; ++itype)
    {
        rcut_max_ = std::max(rcut_max_, radset_[itype]->rcut_max());
    }
}

void RadialCollection::cleanup()
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        delete radset_[itype];
    }
    delete[] radset_;
    radset_ = nullptr;

    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    iter_ = nullptr;

    delete[] nl_;
    nl_ = nullptr;

    ntype_ = 0;
    lmax_ = -1;
    nchi_ = 0;
    nzeta_max_ = 0;
}

void RadialCollection::iter_build()
{
    /*
     * collect the addresses of NumericalRadial objects from different RadialSet objects
     * so that all NumericalRadial objects can be iterated over in a single loop
     *
     * objects are sorted by l first, by itype next, by izeta last.
     *                                                                                      */
    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    delete[] nl_;

    nl_ = new int[lmax_ + 1];
    iter_ = new const NumericalRadial*[nchi_];

    int i = 0;
    std::fill(nl_, nl_ + lmax_ + 1, 0);
    for (int l = 0; l <= lmax_; ++l)
    {
        for (int itype = 0; itype != ntype_; ++itype)
        {
            for (int izeta = 0; izeta < radset_[itype]->nzeta(l); ++izeta)
            {
                iter_[i] = &radset_[itype]->chi(l, izeta);
                ++i;
                ++nl_[l];
            }
        }
    }
}

void RadialCollection::build(const int ntype, Numerical_Nonlocal* const nls)
{
    cleanup();
    ntype_ = ntype;
    radset_ = new RadialSet*[ntype_];

    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = new BetaRadials;
        radset_[itype]->build(nls[itype], itype);

        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    iter_build();
    set_rcut_max();
}

void RadialCollection::build(const int nfile, const std::string* const file, const char file_type)
{
#ifdef __DEBUG
    //assert(file_type == 'o' || file_type == 'p');
    assert(file_type == 'o'); // pseudopotential files are not read in this module
#endif

    cleanup();

    ntype_ = nfile;
    radset_ = new RadialSet*[ntype_];
    switch (file_type)
    {
    case 'o':
        for (int itype = 0; itype < ntype_; ++itype)
        {
            radset_[itype] = new AtomicRadials;
            radset_[itype]->build(file[itype], itype);
        }
        break;
    //case 'p':
    //    for (int itype = 0; itype < ntype_; ++itype)
    //    {
    //        radset_[itype] = new BetaRadials;
    //        radset_[itype]->build(file[itype], itype);
    //    }
    //    break;
    default:; /* not supposed to happen */
    }

    for (int itype = 0; itype < ntype_; ++itype)
    {
        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    iter_build();
    set_rcut_max();
}

void RadialCollection::build(const int ntype, const double* const charges, const int* const nmax, const std::string strategy)
{
    cleanup();
    ntype_ = ntype;
    radset_ = new RadialSet*[ntype_];

    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = new HydrogenRadials;
        radset_[itype]->build(itype, 
                              charges[itype], 
                              nmax[itype], 
                              10.0,             // rcut should be determined automatically, in principle...
                              0.01,
                              1e-6,
                              0,
                              strategy);

        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    // what are these two functions for? Do I need them?
    iter_build();
    set_rcut_max();
}

void RadialCollection::set_transformer(ModuleBase::SphericalBesselTransformer sbt, const int update)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_transformer(sbt, update);
    }
}

void RadialCollection::set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_grid(for_r_space, ngrid, grid, mode);
    }
    rcut_max_ = grid[ngrid - 1];
}

void RadialCollection::set_uniform_grid(const bool for_r_space,
                                        const int ngrid,
                                        const double cutoff,
                                        const char mode,
                                        const bool enable_fft)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_uniform_grid(for_r_space, ngrid, cutoff, mode, enable_fft);
    }
    rcut_max_ = cutoff;
}
