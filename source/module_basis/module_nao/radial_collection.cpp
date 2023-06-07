#include "module_basis/module_nao/radial_collection.h"

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/atomic_radials.h"
#include "module_basis/module_nao/beta_radials.h"

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
}

void RadialCollection::cleanup()
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        delete radset_[itype];
    }
    delete[] radset_;
    radset_ = nullptr;
    ntype_ = 0;
    lmax_ = -1;
    nchi_ = 0;
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
    }
}

void RadialCollection::set_transformer(ModuleBase::SphericalBesselTransformer* sbt, int update)
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

void RadialCollection::set_uniform_grid(const bool for_r_space, const int ngrid, const double cutoff, const char mode)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_uniform_grid(for_r_space, ngrid, cutoff, mode);
    }
}
