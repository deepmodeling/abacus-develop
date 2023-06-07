#include "module_basis/module_nao/radial_set.h"

#include <algorithm>

RadialSet::~RadialSet()
{
    delete[] nzeta_;
    delete[] chi_;
    delete[] index_map_;
}

double RadialSet::rcut_max() const
{
    double rmax = 0.0;
    for (int i = 0; i != nchi_; ++i)
    {
        if (chi_[i].rcut() > rmax)
        {
            rmax = chi_[i].rcut();
        }
    }
    return rmax;
}

RadialSet::RadialSet(const RadialSet& other)
{
    symbol_ = other.symbol_;
    itype_ = other.itype_;
    lmax_ = other.lmax_;

    nzeta_ = nullptr;
    if (lmax_ >= 0)
    {
        nzeta_ = new int[lmax_ + 1];
        for (int l = 0; l <= lmax_; l++)
            nzeta_[l] = other.nzeta_[l];
    }
    nzeta_max_ = other.nzeta_max_;
    nchi_ = other.nchi_;

    chi_ = nullptr;
    if (nchi_ > 0)
    {
        chi_ = new NumericalRadial[nchi_];
        for (int i = 0; i < nchi_; i++)
        {
            chi_[i] = other.chi_[i]; // deep copy
        }
    }

    index_map_ = nullptr;
    if (lmax_ >= 0 && nzeta_max_ > 0)
    {
        index_map_ = new int[(lmax_ + 1) * nzeta_max_];
        for (int i = 0; i < (lmax_ + 1) * nzeta_max_; i++)
        {
            index_map_[i] = other.index_map_[i];
        }
    }
}

RadialSet& RadialSet::operator=(const RadialSet& rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    symbol_ = rhs.symbol_;
    itype_ = rhs.itype_;
    lmax_ = rhs.lmax_;

    delete[] nzeta_;
    nzeta_ = nullptr;
    if (lmax_ >= 0)
    {
        nzeta_ = new int[lmax_ + 1];
        for (int l = 0; l <= lmax_; l++)
            nzeta_[l] = rhs.nzeta_[l];
    }
    nzeta_max_ = rhs.nzeta_max_;
    nchi_ = rhs.nchi_;

    delete[] chi_;
    chi_ = nullptr;
    if (nchi_ > 0)
    {
        chi_ = new NumericalRadial[nchi_];
        for (int i = 0; i < nchi_; i++)
        {
            chi_[i] = rhs.chi_[i]; // deep copy
        }
    }

    delete[] index_map_;
    index_map_ = nullptr;
    if (lmax_ >= 0 && nzeta_max_ > 0)
    {
        index_map_ = new int[(lmax_ + 1) * nzeta_max_];
        for (int i = 0; i < (lmax_ + 1) * nzeta_max_; i++)
        {
            index_map_[i] = rhs.index_map_[i];
        }
    }
    return *this;
}

int RadialSet::index(int l, int izeta) const
{
    assert(l >= 0 && l <= lmax_);
    assert(izeta >= 0 && izeta < nzeta_[l]);
    return index_map_[l * nzeta_max_ + izeta];
}

void RadialSet::indexing()
{
    if (!nzeta_)
    {
        return;
    }

    assert(lmax_ >= 0);
    nzeta_max_ = *std::max_element(nzeta_, nzeta_ + lmax_ + 1);

    delete[] index_map_;
    index_map_ = new int[(lmax_ + 1) * nzeta_max_];
    int index_chi = 0;
    for (int l = 0; l <= lmax_; ++l)
    {
        for (int izeta = 0; izeta != nzeta_max_; ++izeta)
        {
            if (izeta >= nzeta_[l])
            {
                index_map_[l * nzeta_max_ + izeta] = -1; // -1 means no such orbital
            }
            else
            {
                index_map_[l * nzeta_max_ + izeta] = index_chi;
                ++index_chi;
            }
        }
    }
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

    delete[] chi_;
    chi_ = nullptr;

    delete[] index_map_;
    index_map_ = nullptr;
}
