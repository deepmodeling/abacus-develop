#include "module_basis/module_nao/two_center_table.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>

#include "module_base/constants.h"
#include "module_base/math_integral.h"

void TwoCenterTable::build(const RadialCollection& bra,
                           const RadialCollection& ket,
                           const char op,
                           const int nr,
                           const double cutoff,
                           const bool with_deriv)
{
    assert(nr >= 4 && cutoff > 0.0); // nr >= 4 required for polynomial interpolation

    cleanup();

    op_ = op;
    with_deriv_ = with_deriv;

    nr_ = nr;
    rmax_ = cutoff;

    double* rgrid = new double[nr_];
    double dr = rmax_ / (nr_ - 1);
    std::for_each(rgrid, rgrid + nr_, [&rgrid, dr](double& r) { r = (&r - rgrid) * dr; });

    // index the table by generating a map from (itype1, l1, izeta1, itype2, l2, izeta2, l) to a row index
    index_map_.resize({bra.ntype(), bra.lmax() + 1, bra.nzeta_max(),
                       ket.ntype(), ket.lmax() + 1, ket.nzeta_max(), bra.lmax() + ket.lmax() + 1});
    std::fill(index_map_.data<int>(), index_map_.data<int>() + index_map_.NumElements(), -1);

    ntab_ = 0;
    two_center_loop(bra, ket, &TwoCenterTable::_indexing);

    table_.resize({ntab_, nr_});
    two_center_loop(bra, ket, &TwoCenterTable::_tabulate);

    if (with_deriv_)
    {
        // NOTE: for better performance, the loop of _tabulate_d should not be combined with _tabulate.
        // This is due to the caching mechanism of SphericalBesselTransformer.
        dtable_.resize({ntab_, nr_});
        two_center_loop(bra, ket, &TwoCenterTable::_tabulate_d);
    }
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
    assert(is_present(itype1, l1, izeta1, itype2, l2, izeta2, l));
    return deriv ? dtable_.inner_most_ptr<double>(index_map_.get_value<int>(itype1, l1, izeta1, itype2, l2, izeta2, l)):
                    table_.inner_most_ptr<double>(index_map_.get_value<int>(itype1, l1, izeta1, itype2, l2, izeta2, l));
}

double TwoCenterTable::lookup(const int itype1,
                              const int l1,
                              const int izeta1,
                              const int itype2,
                              const int l2,
                              const int izeta2,
                              const int l,
                              const double R,
                              const bool deriv) const
{
    assert(R >= 0);

    if (R > rmax())
    {
        return 0.0;
    }

    const double* tab = table(itype1, l1, izeta1, itype2, l2, izeta2, l, deriv);
    double dr = rmax_ / (nr_ - 1);

    // find the maximum ir satisfying ir * dr <= R
    int ir = static_cast<int>(R / dr);

    // adjust ir so that (ir-1, ir, ir+1, ir+2) all fall within the boundary
    ir += (ir == 0);
    ir -= (ir >= nr_ - 2) + (ir == nr_ -1);

    // apply Lagrange interpolation to data points with indices ir-1, ir, ir+1, ir+2
    double dx0 = R / dr - ir + 1.0;
    double dx1 = dx0 - 1.0;
    double dx2 = dx1 - 1.0;
    double dx3 = dx2 - 1.0;

    return dx1 * dx2 * (tab[ir+2] * dx0 - tab[ir-1] * dx3) / 6.0 +
           dx0 * dx3 * (tab[ir  ] * dx2 - tab[ir+1] * dx1) / 2.0;
}

int& TwoCenterTable::table_index(const NumericalRadial* it1, const NumericalRadial* it2, const int l)
{
    return index_map_.get_value<int>(it1->itype(), it1->l(), it1->izeta(), it2->itype(), it2->l(), it2->izeta(), l);
}

void TwoCenterTable::cleanup()
{
    op_ = '\0';
    with_deriv_ = false;
    nr_ = 0;

    table_.resize({0});
    dtable_.resize({0});
    index_map_.resize({0});
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
    double* tab = table_.inner_most_ptr<double>(table_index(it1, it2, l));
    it1->radtab(op_, *it2, l, tab, nr_, rmax_, false);

    // NOTE: 
    // A radial table stores S(R)/R^l or T(R)/R^l instead of bare S/T.
    //
    // When calculating two-center integrals, we need the product between
    // S(or T) and Y. Note that spherical harmonics are given as R^l * Y_lm,
    // thus the following limit (or the T(R) counterpart)
    //
    //                        S(R)
    //                  lim   ----
    //                  R->0    l
    //                         R
    //
    // is explicited required and should be stored in the table. In order
    // to maintain the consistency of the table, we choose to store S(R)/R^l
    // not only for R->0 but also for R>0.
    //
    // See the developer's document for more details.
    double dr = rmax_ / (nr_ - 1);
    if ( l > 0 )
    {
        //std::transform(&tab[1], tab + nr_, &rgrid[1], &tab[1],
        //        [&](double val, double r) { return val / std::pow(r, l); });
        std::for_each(&tab[1], tab + nr_, [&](double& val) { val /= std::pow(dr * (&val - tab), l); });

        // special treatment for R=0
        int nk = it1->nk();
        const double* kgrid = it1->ptr_kgrid();

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
            fk[ik] = it1->ptr_kvalue()[ik] * it2->ptr_kvalue()[ik] 
                    * std::pow(kgrid[ik], op_exp);
        }

        tab[0] = ModuleBase::Integral::simpson(nk, fk, &h[1]) 
                * ModuleBase::FOUR_PI / dfact(2 * l + 1);

        delete[] fk;
        delete[] h;
    }
}

void TwoCenterTable::_tabulate_d(const NumericalRadial* it1, const NumericalRadial* it2, const int l)
{
    // The R->0 issue does not occur in the derivative calculation
    // because the term which involves dS/dR or dT/dR vanishes at R=0.
    // Therefore, there's no need to find the limit of (dS/dR)/R^(l-1) at R=0,
    // and we choose to store dS/dR or dT/dR directly.
    //
    // See the developer's document for more details.
    it1->radtab(op_, *it2, l, dtable_.inner_most_ptr<double>(table_index(it1, it2, l)), nr_, rmax_, true);
}

