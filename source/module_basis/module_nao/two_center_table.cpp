#include "module_basis/module_nao/two_center_table.h"

#include <algorithm>
#include <cstring>
#include <iostream>

#include "module_base/spherical_bessel_transformer.h"

TwoCenterTable::~TwoCenterTable() {
    delete[] rgrid_;
#ifdef __CONTAINER

#else
    delete[] table_;
    delete[] index_map_;
#endif
}

void TwoCenterTable::build(const RadialCollection& bra, const RadialCollection& ket, const char op, const bool deriv, const int ntab, const double* const tabgrid) {
    // if ntab > 0, tabgrid must be provided; otherwise tabgrid must be nullptr
    assert( (ntab > 0 && tabgrid) || (ntab == 0 && !tabgrid) );

    cleanup();

    op_ = op;
    is_deriv_ = deriv;

    // determine the table grid
    nr_ = ntab > 0 ? ntab : (*bra.cbegin())->nr();
    rgrid_ = new double[nr_];

    if (tabgrid) {
        std::memcpy(rgrid_, tabgrid, nr_ * sizeof(double));
    } else {
        std::memcpy(rgrid_, (*bra.cbegin())->ptr_rgrid(), nr_ * sizeof(double));
    }

    // index the table by generating a map from (itype1, l1, izeta1, itype2, l2, izeta2, l) to a row index
#ifdef __CONTAINER
    index_map_.resize( {bra.ntype(), bra.lmax()+1, bra.nzeta_max(), 
                        ket.ntype(), ket.lmax()+1, ket.nzeta_max(), bra.lmax()+ket.lmax()+1} );
    std::fill(index_map_.data<int>(), index_map_.data<int>() + index_map_.NumElements(), -1);
#else
    ntype1_ = bra.ntype();
    lmax1_ = bra.lmax();
    nzeta_max1_ = bra.nzeta_max();
    ntype2_ = ket.ntype();
    lmax2_ = ket.lmax();
    nzeta_max2_ = ket.nzeta_max();
    lmax_ = lmax1_ + lmax2_;

    // C-style
    stride_[6] = 1;
    stride_[5] = stride_[6] * (lmax_ + 1);
    stride_[4] = stride_[5] * nzeta_max2_;
    stride_[3] = stride_[4] * (lmax2_ + 1);
    stride_[2] = stride_[3] * ntype2_;
    stride_[1] = stride_[2] * nzeta_max1_;
    stride_[0] = stride_[1] * (lmax1_ + 1);

    int sz_map = ntype1_ * (lmax1_+1) * nzeta_max1_ * ntype2_ * (lmax2_+1) * nzeta_max2_ * (lmax_+1);
    index_map_ = new int[sz_map];
    std::fill(index_map_, index_map_ + sz_map, -1);
#endif

    int irow = 0;
    for (const NumericalRadial** it1 = bra.cbegin(); it1 != bra.cend(); ++it1) {
        for (const NumericalRadial** it2 = ket.cbegin(); it2 != ket.cend(); ++it2) {
            // Gaunt coefficients are non-zero if |l1-l2| <= l <= l1+l2 && mod(l+l1+l2, 2) == 0
            for (int l = std::abs((*it1)->l() - (*it2)->l()); l <= (*it1)->l() + (*it2)->l(); l += 2) {
                table_index(*it1, *it2, l) = irow; 
                ++irow;
            }
        }
    }

    // irow is now the number of rows in the table
#ifdef __CONTAINER
    table_.resize({irow, nr_});
#else
    table_ = new double[irow * nr_];
#endif

    // fill in the table
    for (const NumericalRadial** it1 = bra.cbegin(); it1 != bra.cend(); ++it1) {
        for (const NumericalRadial** it2 = ket.cbegin(); it2 != ket.cend(); ++it2) {
            for (int l = std::abs((*it1)->l() - (*it2)->l()); l <= (*it1)->l() + (*it2)->l(); l += 2) {
#ifdef __CONTAINER
                (*it1)->radtab(op, **it2, l, table_.inner_most_ptr<double>(table_index(*it1, *it2, l)), deriv, nr_, rgrid_); 
#else
                (*it1)->radtab(op, **it2, l, &table_[table_index(*it1, *it2, l) * nr_], deriv, nr_, rgrid_);
#endif
            }
        }
    }
}

const double* TwoCenterTable::ptr_table(const int itype1, const int l1, const int izeta1, const int itype2, const int l2, const int izeta2, const int l) const {
    assert( is_present(itype1, l1, izeta1, itype2, l2, izeta2, l) );
#ifdef __CONTAINER
    return table_.inner_most_ptr<double>(index_map_.get_value<int>(itype1, l1, izeta1, itype2, l2, izeta2, l));
#else
    return &table_[index_map_[vectorized_index(itype1, l1, izeta1, itype2, l2, izeta2, l)] * nr_];
#endif
}

#ifdef __CONTAINER

#else
int TwoCenterTable::vectorized_index(const int itype1, const int l1, const int izeta1, const int itype2, const int l2, const int izeta2, const int l) const
{
    return itype1 * stride_[0] + l1 * stride_[1] + izeta1 * stride_[2] + 
        itype2 * stride_[3] + l2 * stride_[4] + izeta2 * stride_[5] + l * stride_[6];
}
#endif

int& TwoCenterTable::table_index(const NumericalRadial* it1, const NumericalRadial* it2, const int l) {
#ifdef __CONTAINER
    return index_map_.get_value<int>(it1->itype(), it1->l(), it1->izeta(), it2->itype(), it2->l(), it2->izeta(), l);
#else
    return index_map_[vectorized_index(it1->itype(), it1->l(), it1->izeta(), it2->itype(), it2->l(), it2->izeta(), l)];
#endif
}

void TwoCenterTable::cleanup() {
    op_ = '\0';
    is_deriv_ = false;
    nr_ = 0;

    delete[] rgrid_;
    rgrid_ = nullptr;

#ifdef __CONTAINER
    table_.resize({0});
    index_map_.resize({0});
#else
    delete[] table_;
    delete[] index_map_;
    table_ = nullptr;
    index_map_ = nullptr;

    ntype1_ = 0;
    lmax1_ = -1;
    nzeta_max1_ = 0;
    ntype2_ = 0;
    lmax2_ = -1;
    nzeta_max2_ = 0;
    lmax_ = -1;

    std::fill(stride_, stride_ + 7, 0);
#endif
}

bool TwoCenterTable::is_present(const int itype1, const int l1, const int izeta1, const int itype2, const int l2, const int izeta2, const int l) const {
    // The given indices map to an entry in the table if they fall within the bounds of index_map_ and 
    // the value of the entry in index_map_ is non-negative
#ifdef __CONTAINER
    return itype1 >= 0 && itype1 < index_map_.shape().dim_size(0) && 
           l1 >= 0 && l1 < index_map_.shape().dim_size(1) && 
           izeta1 >= 0 && izeta1 < index_map_.shape().dim_size(2) &&
           itype2 >= 0 && itype2 < index_map_.shape().dim_size(3) && 
           l2 >= 0 && l2 < index_map_.shape().dim_size(4) && 
           izeta2 >= 0 && izeta2 < index_map_.shape().dim_size(5) &&
           l >= 0 && l <= index_map_.shape().dim_size(6) &&
           index_map_.get_value<int>(itype1, l1, izeta1, itype2, l2, izeta2, l) >= 0;
#else
    return itype1 >= 0 && itype1 < ntype1_ && l1 >= 0 && l1 <= lmax1_ && izeta1 >= 0 && izeta1 < nzeta_max1_ && 
           itype2 >= 0 && itype2 < ntype2_ && l2 >= 0 && l2 <= lmax2_ && izeta2 >= 0 && izeta2 < nzeta_max2_ &&
           l >= 0 && l <= lmax_ && index_map_[vectorized_index(itype1, l1, izeta1, itype2, l2, izeta2, l)] >= 0;
#endif
}
