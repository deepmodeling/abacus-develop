#include "module_basis/module_nao/radial_set.h"

RadialSet::~RadialSet()
{
    delete[] nzeta_;
    delete[] chi_;
    delete[] index_map_;
}

int RadialSet::index(int l, int izeta) const
{
    assert(l >= 0 && l <= lmax_);
    assert(izeta >= 0 && izeta < nzeta_[l]);
    return index_map_[l * nzeta_max_ + izeta];
}

const NumericalRadial& RadialSet::chi(int l, int izeta)
{
    int i = index_map_[l * nzeta_max_ + izeta];
    assert(i >= 0 && i < nchi_);
    return chi_[i];
}

void RadialSet::set_transformer(ModuleBase::SphericalBesselTransformer* sbt, int update)
{
    for (int i = 0; i < nchi_; i++)
        chi_[i].set_transformer(sbt, update);
}

void RadialSet::set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode)
{
    for (int i = 0; i < nchi_; i++)
        chi_[i].set_grid(for_r_space, ngrid, grid, mode);
}

void RadialSet::set_uniform_grid(const bool for_r_space, const int ngrid, const double cutoff, const char mode)
{
    for (int i = 0; i < nchi_; i++)
        chi_[i].set_uniform_grid(for_r_space, ngrid, cutoff, mode);
}

void RadialSet::cleanup()
{
    symbol_ = "";
    itype_ = 0;
    lmax_ = -1;

    delete[] nzeta_;
    nzeta_ = nullptr;
    nzeta_max_ = 0;
    nchi_ = 0;

    rcut_max_ = 0;

    delete[] chi_;
    chi_ = nullptr;

    delete[] index_map_;
    index_map_ = nullptr;
}
