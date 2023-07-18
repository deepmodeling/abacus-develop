#include "module_basis/module_nao/two_center_table.h"

#include <algorithm>
#include <cstring>
#include <iostream>

TwoCenterTable::~TwoCenterTable() { delete[] rgrid_; }

void TwoCenterTable::build(const RadialCollection& bra,
                           const RadialCollection& ket,
                           const char op,
                           const int nr,
                           const double* const rgrid,
                           const bool deriv)
{
    assert(nr > 0 && rgrid);
    cleanup();

    op_ = op;
    is_deriv_ = deriv;

    nr_ = nr;
    rgrid_ = new double[nr_];
    std::memcpy(rgrid_, rgrid, nr_ * sizeof(double));

    // index the table by generating a map from (itype1, l1, izeta1, itype2, l2, izeta2, l) to a row index
    index_map_.resize({bra.ntype(), bra.lmax() + 1, bra.nzeta_max(),
                       ket.ntype(), ket.lmax() + 1, ket.nzeta_max(), bra.lmax() + ket.lmax() + 1});
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
