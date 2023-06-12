#include "module_basis/module_nao/radial_collection.h"

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/atomic_radials.h"
#include "module_basis/module_nao/beta_radials.h"

RadialCollection::RadialCollection(const RadialCollection& other)
{
    ntype_ = other.ntype_;
    lmax_ = other.lmax_;
    nchi_ = other.nchi_;
    nzeta_max_ = other.nzeta_max_;

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

    radset_ = new RadialSet*[ntype_];
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = rhs.radset_[itype]->clone();
    }

    iter_build();

    return *this;
}

double RadialCollection::rcut_max() const
{
    double rmax = 0.0;
    for (int itype = 0; itype < ntype_; ++itype)
    {
        rmax = std::max(rmax, radset_[itype]->rcut_max());
    }
    return rmax;
}

RadialCollection::~RadialCollection()
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        delete radset_[itype];
    }
    delete[] radset_;
    delete[] iter_; // iterator does not control memory; simply delete the pointer array
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
     *                                                                                      */
    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    iter_ = new const NumericalRadial*[nchi_];
    int i = 0;
    for (int itype = 0; itype < ntype_; ++itype)
    {
        for (const NumericalRadial* it = radset_[itype]->cbegin(); it != radset_[itype]->cend(); ++it)
        {
            iter_[i] = it;
            ++i;
        }
    }
}

void RadialCollection::build(const int nfile, const std::string* const file, const char file_type)
{
    assert(file_type == 'o' || file_type == 'p');

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
    case 'p':
        for (int itype = 0; itype < ntype_; ++itype)
        {
            radset_[itype] = new BetaRadials;
            radset_[itype]->build(file[itype], itype);
        }
        break;
    default:; /* not supposed to happen */
    }

    for (int itype = 0; itype < ntype_; ++itype)
    {
        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    iter_build();
}

void RadialCollection::set_transformer(ModuleBase::SphericalBesselTransformer* const sbt, const int update)
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
}
