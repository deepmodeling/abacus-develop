#include "module_basis/module_nao/two_center_table.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>

TwoCenterTable::~TwoCenterTable() { delete[] rgrid_; }

void TwoCenterTable::build(const RadialCollection& bra,
                           const RadialCollection& ket,
                           const char op,
                           const int nr,
                           const double* const rgrid,
                           const bool deriv)
{
    assert(nr >= 4 && rgrid); // nr >= 4 required for polynomial interpolation
    cleanup();

    op_ = op;
    is_deriv_ = deriv;

    nr_ = nr;
    rgrid_ = new double[nr_];
    std::memcpy(rgrid_, rgrid, nr_ * sizeof(double));

    double tol = 4.0 * std::numeric_limits<double>::epsilon();
    double dr = rgrid[nr - 1] / (nr - 1);
    is_uniform_ = std::all_of(rgrid, rgrid + nr, [&](const double& r) { return std::abs(r - (&r - rgrid) * dr) < tol; });

    // index the table by generating a map from (itype1, l1, izeta1, itype2, l2, izeta2, l) to a row index
    index_map_.resize({bra.ntype(),
                       bra.lmax() + 1,
                       bra.nzeta_max(),
                       ket.ntype(),
                       ket.lmax() + 1,
                       ket.nzeta_max(),
                       bra.lmax() + ket.lmax() + 1});
    std::fill(index_map_.data<int>(), index_map_.data<int>() + index_map_.NumElements(), -1);

    ntab_ = 0;

    for (int l = 0; l <= bra.lmax() + ket.lmax(); ++l)
    {
        for (int l1 = 0; l1 <= bra.lmax(); ++l1)
        {
            for (const NumericalRadial** it1 = bra.cbegin(l1); it1 != bra.cend(l1); ++it1)
            {
                for (int l2 = std::abs(l1 - l); l2 <= std::min(ket.lmax(), l + l1); l2 += 2)
                {
                    for (const NumericalRadial** it2 = ket.cbegin(l2); it2 != ket.cend(l2); ++it2)
                    {
                        table_index(*it1, *it2, l) = ntab_;
                        ++ntab_;
                    }
                }
            }
        }
    }

    // irow is now the number of rows in the table
    table_.resize({ntab_, nr_});

    for (int l = 0; l <= bra.lmax() + ket.lmax(); ++l)
    {
        for (int l1 = 0; l1 <= bra.lmax(); ++l1)
        {
            for (const NumericalRadial** it1 = bra.cbegin(l1); it1 != bra.cend(l1); ++it1)
            {
                for (int l2 = std::abs(l1 - l); l2 <= std::min(ket.lmax(), l + l1); l2 += 2)
                {
                    for (const NumericalRadial** it2 = ket.cbegin(l2); it2 != ket.cend(l2); ++it2)
                    {
                        (*it1)->radtab(op,
                                       **it2,
                                       l,
                                       table_.inner_most_ptr<double>(table_index(*it1, *it2, l)),
                                       nr_,
                                       rgrid_,
                                       deriv);
                    }
                }
            }
        }
    }
}

const double* TwoCenterTable::ptr_table(const int itype1,
                                        const int l1,
                                        const int izeta1,
                                        const int itype2,
                                        const int l2,
                                        const int izeta2,
                                        const int l) const
{
    assert(is_present(itype1, l1, izeta1, itype2, l2, izeta2, l));
    return table_.inner_most_ptr<double>(index_map_.get_value<int>(itype1, l1, izeta1, itype2, l2, izeta2, l));
}

double TwoCenterTable::lookup(const int itype1,
                              const int l1,
                              const int izeta1,
                              const int itype2,
                              const int l2,
                              const int izeta2,
                              const int l,
                              const double R) const
{
    assert(R >= 0 && R <= rmax());
    const double* table = ptr_table(itype1, l1, izeta1, itype2, l2, izeta2, l);

    // find the maximum ir satisfying rgrid_[ir] <= R
    int ir = is_uniform_ ? static_cast<int>(R / rgrid_[1]) : std::upper_bound(rgrid_, rgrid_ + nr_, R) - rgrid_ - 1;

    // adjust ir so that (ir-1, ir, ir+1, ir+2) all fall within the boundary
    ir += (ir == 0);
    ir = ir >= nr_ - 2 ? nr_ - 3 : ir;

    // use Lagrange interpolation for data points with indices ir-1, ir, ir+1, ir+2
    if (is_uniform_)
    {
        // dxi = (R - rgrid_[ir-i-1]) / dr
        double dx0 = (R - rgrid_[ir - 1]) / rgrid_[1];
        double dx1 = dx0 - 1.0;
        double dx2 = dx1 - 1.0;
        double dx3 = dx2 - 1.0;

        // return - table[ir-1] * dx1 * dx2 * dx3 / 6.0
        //        + table[ir  ] * dx0 * dx2 * dx3 / 2.0
        //        - table[ir+1] * dx0 * dx1 * dx3 / 2.0
        //        + table[ir+2] * dx0 * dx1 * dx2 / 6.0;
        return dx1 * dx2 * (-table[ir - 1] * dx3 + table[ir + 2] * dx0) / 6.0
               + dx0 * dx3 * (table[ir] * dx2 - table[ir + 1] * dx1) / 2.0;
    }
    else
    {
        return table[ir - 1] * (R - rgrid_[ir]) * (R - rgrid_[ir + 1]) * (R - rgrid_[ir + 2]) /
            ((rgrid_[ir - 1] - rgrid_[ir]) * (rgrid_[ir - 1] - rgrid_[ir + 1]) * (rgrid_[ir - 1] - rgrid_[ir + 2]))
            + table[ir] * (R - rgrid_[ir - 1]) * (R - rgrid_[ir + 1]) * (R - rgrid_[ir + 2]) /
            ((rgrid_[ir] - rgrid_[ir - 1]) * (rgrid_[ir] - rgrid_[ir + 1]) * (rgrid_[ir] - rgrid_[ir + 2]))
            + table[ir + 1] * (R - rgrid_[ir - 1]) * (R - rgrid_[ir]) * (R - rgrid_[ir + 2]) /
            ((rgrid_[ir + 1] - rgrid_[ir - 1]) * (rgrid_[ir + 1] - rgrid_[ir]) * (rgrid_[ir + 1] - rgrid_[ir + 2]))
            + table[ir + 2] * (R - rgrid_[ir - 1]) * (R - rgrid_[ir]) * (R - rgrid_[ir + 1]) /
            ((rgrid_[ir + 2] - rgrid_[ir - 1]) * (rgrid_[ir + 2] - rgrid_[ir]) * (rgrid_[ir + 2] - rgrid_[ir + 1]));
    }
}

int& TwoCenterTable::table_index(const NumericalRadial* it1, const NumericalRadial* it2, const int l)
{
    return index_map_.get_value<int>(it1->itype(), it1->l(), it1->izeta(), it2->itype(), it2->l(), it2->izeta(), l);
}

void TwoCenterTable::cleanup()
{
    op_ = '\0';
    is_deriv_ = false;
    nr_ = 0;

    delete[] rgrid_;
    rgrid_ = nullptr;

    table_.resize({0});
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
