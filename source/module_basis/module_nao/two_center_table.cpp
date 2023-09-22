#include "module_basis/module_nao/two_center_table.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>

#include "module_base/constants.h"
#include "module_base/math_integral.h"
#include "module_base/cubic_spline.h"

void TwoCenterTable::build(const RadialCollection& bra,
                           const RadialCollection& ket,
                           const char op,
                           const int nr,
                           const double cutoff)
{
#ifdef __DEBUG
    assert(nr >= 3 && cutoff > 0.0);
#endif

    cleanup();

    op_ = op;
    nr_ = nr;
    rmax_ = cutoff;

    nchi_ket_.resize({ket.ntype(), ket.lmax() + 1});
    std::fill(nchi_ket_.data<int>(), nchi_ket_.data<int>() + nchi_ket_.NumElements(), 0);
    for (int itype = 0; itype < ket.ntype(); ++itype)
        for (int l = 0; l <= ket.lmax(itype); ++l)
            nchi_ket_.get_value<int>(itype, l) = ket.nzeta(itype, l);

    rgrid_ = new double[nr_];
    double dr = rmax_ / (nr_ - 1);
    std::for_each(rgrid_, rgrid_ + nr_, [this, dr](double& r) { r = (&r - rgrid_) * dr; });

    // index the table by generating a map from (itype1, l1, izeta1, itype2, l2, izeta2, l) to a row index
    index_map_.resize({bra.ntype(), bra.lmax() + 1, bra.nzeta_max(),
                       ket.ntype(), ket.lmax() + 1, ket.nzeta_max(), bra.lmax() + ket.lmax() + 1});
    std::fill(index_map_.data<int>(), index_map_.data<int>() + index_map_.NumElements(), -1);

    ntab_ = 0;
    two_center_loop(bra, ket, &TwoCenterTable::_indexing);

    table_.resize({ntab_, nr_});
    dtable_.resize({ntab_, nr_});
    two_center_loop(bra, ket, &TwoCenterTable::_tabulate);
}

const double* TwoCenterTable::table(const int itype1,
                                    const int l1,
                                    const int izeta1,
                                    const int itype2,
                                    const int l2,
                                    const int izeta2,
                                    const int l,
                                    const bool deriv) const
{
#ifdef __DEBUG
    assert(is_present(itype1, l1, izeta1, itype2, l2, izeta2, l));
#endif
    return deriv ? dtable_.inner_most_ptr<double>(index_map_.get_value<int>(itype1, l1, izeta1, itype2, l2, izeta2, l)):
                    table_.inner_most_ptr<double>(index_map_.get_value<int>(itype1, l1, izeta1, itype2, l2, izeta2, l));
}

void TwoCenterTable::lookup(const int itype1,
                            const int l1,
                            const int izeta1,
                            const int itype2,
                            const int l2,
                            const int izeta2,
                            const int l,
                            const double R,
                            double* val,
                            double* dval) const
{
#ifdef __DEBUG
    assert(R >= 0);
#endif

    if (R > rmax())
    {
        if (val) *val = 0.0;
        if (dval) *dval = 0.0;
        return;
    }

    const double*  tab = table(itype1, l1, izeta1, itype2, l2, izeta2, l, false);
    const double* dtab = table(itype1, l1, izeta1, itype2, l2, izeta2, l, true);
    ModuleBase::CubicSpline::eval(nr_, rgrid_, tab, dtab, 1, &R, val, dval);
}

int& TwoCenterTable::table_index(const NumericalRadial* it1, const NumericalRadial* it2, const int l)
{
    return index_map_.get_value<int>(it1->itype(), it1->l(), it1->izeta(), it2->itype(), it2->l(), it2->izeta(), l);
}

void TwoCenterTable::cleanup()
{
    op_ = '\0';
    nr_ = 0;
    delete[] rgrid_;
    rgrid_ = nullptr;

    table_.resize({0});
    dtable_.resize({0});
    index_map_.resize({0});
    nchi_ket_.resize({0});
}

bool TwoCenterTable::is_present(const int itype1,
                                const int l1,
                                const int izeta1,
                                const int itype2,
                                const int l2,
                                const int izeta2,
                                const int l) const
{
    // The given indices map to an entry in the table if they fall within the bounds of index_map_ and
    // the value of the entry in index_map_ is non-negative
    return itype1 >= 0 && itype1 < index_map_.shape().dim_size(0) && l1 >= 0 && l1 < index_map_.shape().dim_size(1)
           && izeta1 >= 0 && izeta1 < index_map_.shape().dim_size(2) && itype2 >= 0
           && itype2 < index_map_.shape().dim_size(3) && l2 >= 0 && l2 < index_map_.shape().dim_size(4) && izeta2 >= 0
           && izeta2 < index_map_.shape().dim_size(5) && l >= 0 && l <= index_map_.shape().dim_size(6)
           && index_map_.get_value<int>(itype1, l1, izeta1, itype2, l2, izeta2, l) >= 0;
}

double TwoCenterTable::dfact(int l) const
{
    double result = 1.0;
    for (int i = l; i > 1; i -= 2)
    {
        result *= i;
    }
    return result;
}

void TwoCenterTable::two_center_loop(const RadialCollection& bra, 
                                     const RadialCollection& ket, 
                                     looped_func f)
{
    for (int l = 0; l <= bra.lmax() + ket.lmax(); ++l)
        for (int l1 = 0; l1 <= bra.lmax(); ++l1)
            for (const NumericalRadial** it1 = bra.cbegin(l1); it1 != bra.cend(l1); ++it1)
                for (int l2 = std::abs(l1 - l); l2 <= std::min(ket.lmax(), l + l1); l2 += 2)
                    for (const NumericalRadial** it2 = ket.cbegin(l2); it2 != ket.cend(l2); ++it2)
                        (this->*f)(*it1, *it2, l);
}

void TwoCenterTable::_indexing(const NumericalRadial* it1, const NumericalRadial* it2, const int l)
{
    table_index(it1, it2, l) = ntab_++;
}

void TwoCenterTable::_tabulate(const NumericalRadial* it1, const NumericalRadial* it2, const int l)
{
    int itab = table_index(it1, it2, l);
    double* tab = table_.inner_most_ptr<double>(itab);
    it1->radtab(op_, *it2, l, tab, nr_, rmax_, false);

    // NOTE:
    // A radial table stores S(R)/R^l or T(R)/R^l instead of bare S/T.
    //
    // When calculating two-center integrals, we need the product between
    // S(or T) and Y. Since spherical harmonics are given as R^l * Y_lm,
    // it is more convenient to store S(R)/R^l or T(R)/R^l and their
    // derivatives instead of the bare ones.
    //
    // Note that the values at R=0 should be understood as, e.g.,
    //
    //                        S(R)
    //                  lim   ----
    //                  R->0    l
    //                         R
    //
    // which needs special treatment.
    //
    // See the developer's document for more details.
    double dr = rmax_ / (nr_ - 1);
    if ( l > 0 )
    {
        // divide S(R) by R^l (except the R=0 point)
        std::for_each(&tab[1], tab + nr_, [&](double& val) { val /= std::pow(dr * (&val - tab), l); });

        // special treatment for R=0
        int nk = it1->nk();
        const double* kgrid = it1->kgrid();

        double* fk = new double[nk];
        double* h = new double[nk];
        std::adjacent_difference(kgrid, kgrid + nk, h);

        int op_exp = l;
        switch (op_)
        {
        case 'S': op_exp += 2;
                  break;
        case 'T': op_exp += 4;
                  break;
        default: ; // currently not supposed to happen
        }

        for (int ik = 0; ik != nk; ++ik)
        {
            fk[ik] = it1->kvalue(ik) * it2->kvalue(ik) 
                    * std::pow(kgrid[ik], op_exp);
        }

        tab[0] = ModuleBase::Integral::simpson(nk, fk, &h[1]) 
                * ModuleBase::FOUR_PI / dfact(2 * l + 1);

        delete[] fk;
        delete[] h;
    }

    // The derivative table stores the derivative of S(R)/R^l or T(R)/R^l 
    // instead of bare dS(R)/dR or dT(R)/dR, which simplifies further calculation.
    //
    // The derivatives are computed from a cubic spline interpolation rather
    // than two spherical Bessel transforms. By doing so, we achieve a good
    // consistency between the table and its derivative during interpolation.
    using ModuleBase::CubicSpline;
    CubicSpline::build(nr_, rgrid_, table_.inner_most_ptr<double>(itab), dtable_.inner_most_ptr<double>(itab),
            CubicSpline::BoundaryCondition::first_deriv, CubicSpline::BoundaryCondition::first_deriv,
            0.0, 0.0);
}

