#include "module_basis/module_nao/real_gaunt_table.h"

#include <array>
#include <cassert>
#include <algorithm>

#include "module_base/constants.h"

void RealGauntTable::build(const int lmax)
{
    assert( lmax >= 0 );

    // do nothing if lmax < lmax_ -- or shall we shrink the map & tensor?
    if (lmax <= lmax_)
    {
        return;
    }

    // build the standard Gaunt table (with symmetry & selection rule considered)
    for (int l1 = 0; l1 <= 2 * lmax; ++l1)
    {
        for (int l2 = 0; l2 <= l1; ++l2)
        {
            for (int l3 = l1 - l2; l3 <= l2; l3 += 2)
            {
                for (int m3 = 0; m3 <= l3; ++m3)
                {
                    int m2_max = (l1 - m3) > l2 ? l2 : l1 - m3; // ensure min(m1) >= l1
                    for (int m2 = -l2; m2 <= m2_max; ++m2)
                    {
                        int m1 = -m2 - m3;
                        gaunt_table_[ std::array<int, 6>({l1, l2, l3, m1, m2, m3}) ] = gaunt(l1, l2, l3, m1, m2, m3);
                    }
                }
            }
        }
    }

    lmax_ = lmax;

    // build the real Gaunt table from tabulated standard Gaunt coefficients
    // this real Gaunt table is supposed to be used in two-center integrals, so the maximum
    // l of the third dimension is twice as large as the maximum l of the first two dimensions
    real_gaunt_table_.resize({ (lmax + 1) * (lmax + 1), (lmax + 1) * (lmax + 1), (2 * lmax + 1) * (2 * lmax + 1) });
    real_gaunt_table_.zero();

    for (int l1 = 0; l1 <= lmax; ++l1)
    {
        for (int m1 = -l1; m1 <= l1; ++m1)
        {
            int index1 = index_map(l1, m1);
            for (int l2 = 0; l2 <= lmax; ++l2)
            {
                for (int m2 = -l2; m2 <= l2; ++m2)
                {
                    int index2 = index_map(l2, m2);
                    for (int l3 = std::abs(l1 - l2); l3 <= l1 + l2; l3 += 2)
                    {
                        for (int m3 = -l3; m3 <= l3; ++m3)
                        {
                            int index3 = index_map(l3, m3);
                            real_gaunt_table_.get_value<double>(index1, index2, index3)
                                = real_gaunt_lookup(l1, l2, l3, m1, m2, m3);
                        }
                    }
                }
            }
        }
    }

}

const double& RealGauntTable::operator()(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const
{
    assert( is_valid_lm(l1, l2, l3, m1, m2, m3) );
    assert( l1 <= lmax_ && l2 <= lmax_ && l3 <= 2 * lmax_ );
    return real_gaunt_table_.get_value<double>(index_map(l1, m1), index_map(l2, m2), index_map(l3, m3));
}

double RealGauntTable::real_gaunt_lookup(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const
{
    // This function calculates and returns the Gaunt coefficients of real spherical harmonics
    // from tabulated standard Gaunt coefficients.

    assert( is_valid_lm(l1, l2, l3, m1, m2, m3) );
    assert( l1 <= lmax_ && l2 <= lmax_ && l3 <= 2 * lmax_ );

    if ( !gaunt_select_l(l1, l2, l3) || !real_gaunt_select_m(m1, m2, m3) )
    {
        return 0.0;
    }

    std::array<int, 3> m = {std::abs(m1), std::abs(m2), std::abs(m3)};
    int& m_absmax = *std::max_element(m.begin(), m.end());

    if ( m1 == 0 || m2 == 0 || m3 == 0 )
    {
        m_absmax = -m_absmax;
        return minus_1_pow(m_absmax) * gaunt_lookup(l1, l2, l3, m[0], m[1], m[2]);
    }
    else if ( m1 + m2 + m3 == 0 )
    {
        return ModuleBase::SQRT2 / 2.0 * minus_1_pow(m_absmax + 1) * gaunt_lookup(l1, l2, l3, m1, m2, m3);

    }
    else
    {
        m_absmax = -m_absmax;
        return ModuleBase::SQRT2 / 2.0 * minus_1_pow(m_absmax) * gaunt_lookup(l1, l2, l3, m[0], m[1], m[2]);
    }
}

double RealGauntTable::gaunt(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const
{
    // This function computes the Gaunt coefficients from the Wigner-3j expression

    assert( is_valid_lm(l1, l2, l3, m1, m2, m3) );
    if ( !gaunt_select_l(l1, l2, l3) || !gaunt_select_m(m1, m2, m3) )
    {
        return 0.0;
    }

    int g = (l1 + l2 + l3) / 2;
    double pref = std::sqrt( (2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / ModuleBase::FOUR_PI);
    double tri = std::sqrt( factorial(l1 + l2 - l3) * factorial(l2 + l3 - l1) * factorial(l3 + l1 - l2)
                            / factorial(l1+l2+l3+1) );

    // wigner3j(l1,l2,l3,0,0,0)
    double wigner1 = minus_1_pow(g) * tri * factorial(g) / factorial(g-l1) / factorial(g-l2) / factorial(g-l3);

    // wigner3j(l1,l2,l3,m1,m2,m3)
    int kmin = std::max(l2 - l3 - m1, l1 - l3 + m2);
    kmin = std::max(kmin, 0);

    int kmax = std::min(l1 - m1, l2 + m2);
    kmax = std::min(kmax, l1 + l2 - l3);

    double wigner2 = 0.0;
    for (int k = kmin; k <= kmax; ++k)
    {
        wigner2 += minus_1_pow(k) / factorial(k) / factorial(l1 - m1 - k) / factorial(l2 + m2 - k)
            / factorial(l3 - l2 + m1 + k) / factorial(l3 - l1 - m2 + k) / factorial(l1 + l2 - l3 - k);
    }

    wigner2 *= tri * minus_1_pow(l1 - l2 - m3) * std::sqrt(
            factorial(l1 + m1) * factorial(l1 - m1) *
            factorial(l2 + m2) * factorial(l2 - m2) *
            factorial(l3 + m3) * factorial(l3 - m3) );

    return pref * wigner1 * wigner2;
}

bool RealGauntTable::is_valid_lm(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const
{
    return std::abs(m1) <= l1 && std::abs(m2) <= l2 && std::abs(m3) <= l3;
}

bool RealGauntTable::gaunt_select_l(const int l1, const int l2, const int l3) const
{
    return l1 + l2 >= l3 && l1 + l3 >= l2 && l2 + l3 >= l1 && (l1 + l2 + l3) % 2 == 0;
}

bool RealGauntTable::real_gaunt_select_m(const int m1, const int m2, const int m3) const
{
    return  ( ( static_cast<int>(m1 < 0) + static_cast<int>(m2 < 0) + static_cast<int>(m3 < 0) ) % 2 == 0 ) && 
        ( std::abs(m1) + std::abs(m2) == std::abs(m3) || 
          std::abs(m2) + std::abs(m3) == std::abs(m1) || 
          std::abs(m3) + std::abs(m1) == std::abs(m2) );
}

double RealGauntTable::gaunt_lookup(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const
{
    assert( is_valid_lm(l1, l2, l3, m1, m2, m3) );
    assert( l1 <= 2 * lmax_ && l2 <= 2 * lmax_ && l3 <= 2 * lmax_ );

    return ( gaunt_select_l(l1, l2, l3) && gaunt_select_m(m1, m2, m3) ) ?
        gaunt_table_.at( gaunt_key(l1, l2, l3, m1, m2, m3) ) : 0.0;
}

std::array<int, 6> RealGauntTable::gaunt_key(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const
{
    assert( is_valid_lm(l1, l2, l3, m1, m2, m3) );

    std::array<int, 6> key{l1, l2, l3, m1, m2, m3};
    arrange(key[0], key[1], key[3], key[4]);
    arrange(key[0], key[2], key[3], key[5]);
    arrange(key[1], key[2], key[4], key[5]);
    if ( key[5] < 0 )
    {
        key[3] = -key[3];
        key[4] = -key[4];
        key[5] = -key[5];
    }
    return key;
}

void RealGauntTable::arrange(int& l1, int& l2, int& m1, int& m2) const
{
    if ( l1 < l2 )
    {
        std::swap(l1, l2);
        std::swap(m1, m2);
    }
}

double RealGauntTable::factorial(const int n) const
{
    assert( n >= 0 );
    double val = 1.0;
    for(int i = 2; i <= n; i++)
    {
        val *= static_cast<double>(i);
    }
    return val;
}

int RealGauntTable::index_map(int l, int m) const
{
    assert( std::abs(m) <= l );
    return l * l + l + m;
}
