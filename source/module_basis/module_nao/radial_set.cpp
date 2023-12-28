#include "module_basis/module_nao/radial_set.h"

#include <algorithm>
#include <cstring>
#include <memory>

#include "module_base/spherical_bessel_transformer.h"

RadialSet::~RadialSet()
{
    delete[] nzeta_;
    delete[] chi_;
    delete[] index_map_;
}

RadialSet::RadialSet(const RadialSet& other) :
    symbol_(other.symbol_),
    itype_(other.itype_),
    lmax_(other.lmax_),
    rcut_max_(other.rcut_max_),
    nzeta_(nullptr),
    nzeta_max_(other.nzeta_max_),
    nchi_(other.nchi_),
    chi_(nullptr),
    index_map_(nullptr)
{
    if (nchi_ == 0)
    {
        return;
    }

    nzeta_ = new int[lmax_ + 1];
    std::memcpy(nzeta_, other.nzeta_, (lmax_ + 1) * sizeof(int));

    index_map_ = new int[(lmax_ + 1) * nzeta_max_];
    std::memcpy(index_map_, other.index_map_, (lmax_ + 1) * nzeta_max_ * sizeof(int));

    chi_ = new NumericalRadial[nchi_];
    for (int i = 0; i < nchi_; i++)
    {
        chi_[i] = other.chi_[i]; // deep copy
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
    rcut_max_ = rhs.rcut_max_;
    nzeta_max_ = rhs.nzeta_max_;
    nchi_ = rhs.nchi_;

    delete[] nzeta_;
    delete[] chi_;
    delete[] index_map_;
    nzeta_ = nullptr;
    chi_ = nullptr;
    index_map_ = nullptr;

    if (nchi_ > 0)
    {
        nzeta_ = new int[lmax_ + 1];
        std::memcpy(nzeta_, rhs.nzeta_, (lmax_ + 1) * sizeof(int));

        index_map_ = new int[(lmax_ + 1) * nzeta_max_];
        std::memcpy(index_map_, rhs.index_map_, (lmax_ + 1) * nzeta_max_ * sizeof(int));

        chi_ = new NumericalRadial[nchi_];
        for (int i = 0; i < nchi_; i++)
        {
            chi_[i] = rhs.chi_[i]; // deep copy
        }
    }

    return *this;
}

void RadialSet::to_numerical_orbital(Numerical_Orbital& no, const int nk_legacy, const double lcao_dk) const
{
    delete[] no.chi();

    no.chi() = new Numerical_Orbital_Lm[nchi_];
    for (int i = 0; i < nchi_; i++)
    {
        chi_[i].to_numerical_orbital_lm(no.chi()[i], nk_legacy, lcao_dk);
    }

    no.set_orbital_info(itype_, symbol_, lmax_, nzeta_, nchi_);
}

void RadialSet::set_rcut_max()
{
    rcut_max_ = 0.0;
    for (int i = 0; i < nchi_; ++i)
    {
        rcut_max_ = std::max(rcut_max_, chi_[i].rcut());
    }
}

int RadialSet::index(const int l, const int izeta) const
{
#ifdef __DEBUG
    assert(l >= 0 && l <= lmax_);
    assert(izeta >= 0 && izeta < nzeta_[l]);
#endif
    return index_map_[l * nzeta_max_ + izeta];
}

void RadialSet::indexing()
{
    if (!nzeta_)
    {
        return;
    }

#ifdef __DEBUG
    assert(lmax_ >= 0);
#endif

    delete[] index_map_;
    index_map_ = new int[(lmax_ + 1) * nzeta_max_];
    int index_chi = 0;
    for (int l = 0; l <= lmax_; ++l)
    {
        for (int izeta = 0; izeta != nzeta_max_; ++izeta)
        {
            index_map_[l * nzeta_max_ + izeta] = izeta >= nzeta_[l] ? -1 : index_chi++;
        }
    }
}

const NumericalRadial& RadialSet::chi(const int l, const int izeta)
{
    int i = index_map_[l * nzeta_max_ + izeta];
#ifdef __DEBUG
    assert(i >= 0 && i < nchi_);
#endif
    return chi_[i];
}

void RadialSet::set_transformer(ModuleBase::SphericalBesselTransformer sbt, const int update)
{
    for (int i = 0; i < nchi_; i++)
    {
        chi_[i].set_transformer(sbt, update);
    }
}

void RadialSet::set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode)
{
    for (int i = 0; i < nchi_; i++)
    {
        chi_[i].set_grid(for_r_space, ngrid, grid, mode);
    }
    set_rcut_max();
}

void RadialSet::set_uniform_grid(const bool for_r_space,
                                 const int ngrid,
                                 const double cutoff,
                                 const char mode,
                                 const bool enable_fft)
{
    for (int i = 0; i < nchi_; i++)
    {
        chi_[i].set_uniform_grid(for_r_space, ngrid, cutoff, mode, enable_fft);
    }
    set_rcut_max();
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
