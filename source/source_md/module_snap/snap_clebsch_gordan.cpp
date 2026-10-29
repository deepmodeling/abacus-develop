#include "snap_clebsch_gordan.h"

#include <cmath>

namespace ModuleMD
{
namespace SNAP
{

namespace
{
int half_sum(const int left, const int right)
{
    return (left + right) / 2;
}
}

ClebschGordanTable::ClebschGordanTable(const int twojmax) : twojmax_(twojmax < 0 ? 0 : twojmax)
{
    for (int j1 = 0; j1 <= twojmax_; ++j1)
    {
        for (int j2 = 0; j2 <= twojmax_; ++j2)
        {
            for (int j = std::abs(j1 - j2); j <= std::min(twojmax_, j1 + j2); j += 2)
            {
                for (int m1 = -j1; m1 <= j1; m1 += 2)
                {
                    for (int m2 = -j2; m2 <= j2; m2 += 2)
                    {
                        const int m = m1 + m2;
                        if (m < -j || m > j)
                        {
                            continue;
                        }
                        Entry entry;
                        entry.j1 = j1;
                        entry.m1 = m1;
                        entry.j2 = j2;
                        entry.m2 = m2;
                        entry.j = j;
                        entry.m = m;
                        entry.value = coefficient_uncached(j1, m1, j2, m2, j, m);
                        entries_.push_back(entry);
                    }
                }
            }
        }
    }
}

double ClebschGordanTable::factorial(const int value)
{
    if (value < 0)
    {
        return 0.0;
    }
    return std::tgamma(static_cast<double>(value) + 1.0);
}

bool ClebschGordanTable::allowed(const int j1,
                                 const int m1,
                                 const int j2,
                                 const int m2,
                                 const int j,
                                 const int m)
{
    if (j1 < 0 || j2 < 0 || j < 0 || m1 < -j1 || m1 > j1 || m2 < -j2 || m2 > j2 || m < -j || m > j)
    {
        return false;
    }
    if ((j1 - m1) % 2 != 0 || (j2 - m2) % 2 != 0 || (j - m) % 2 != 0)
    {
        return false;
    }
    if ((j1 + j2 + j) % 2 != 0 || j < std::abs(j1 - j2) || j > j1 + j2)
    {
        return false;
    }
    return m == m1 + m2;
}

double ClebschGordanTable::coefficient_uncached(const int j1,
                                                const int m1,
                                                const int j2,
                                                const int m2,
                                                const int j,
                                                const int m)
{
    if (!allowed(j1, m1, j2, m2, j, m))
    {
        return 0.0;
    }

    const double delta = std::sqrt(
        factorial((j1 + j2 - j) / 2) * factorial((j1 - j2 + j) / 2)
        * factorial((-j1 + j2 + j) / 2) / factorial((j1 + j2 + j) / 2 + 1));
    const double prefactor = std::sqrt(
        (j + 1.0) * factorial(half_sum(j1, m1)) * factorial(half_sum(j1, -m1))
        * factorial(half_sum(j2, m2)) * factorial(half_sum(j2, -m2))
        * factorial(half_sum(j, m)) * factorial(half_sum(j, -m)));

    double sum = 0.0;
    const int zmin = std::max(0, std::max(-(j - j2 + m1) / 2, -(j - j1 - m2) / 2));
    const int zmax = std::min((j1 + j2 - j) / 2,
                              std::min((j1 - m1) / 2, (j2 + m2) / 2));
    for (int z = zmin; z <= zmax; ++z)
    {
        const int a = z;
        const int b = (j1 + j2 - j) / 2 - z;
        const int c = (j1 - m1) / 2 - z;
        const int d = (j2 + m2) / 2 - z;
        const int e = (j - j2 + m1) / 2 + z;
        const int f = (j - j1 - m2) / 2 + z;
        if (a < 0 || b < 0 || c < 0 || d < 0 || e < 0 || f < 0)
        {
            continue;
        }
        const double denominator = factorial(a) * factorial(b) * factorial(c) * factorial(d) * factorial(e) * factorial(f);
        sum += (z % 2 == 0 ? 1.0 : -1.0) / denominator;
    }

    return delta * prefactor * sum;
}

const ClebschGordanTable::Entry* ClebschGordanTable::find(const int j1,
                                                          const int m1,
                                                          const int j2,
                                                          const int m2,
                                                          const int j,
                                                          const int m) const
{
    for (std::size_t index = 0; index < entries_.size(); ++index)
    {
        const Entry& entry = entries_[index];
        if (entry.j1 == j1 && entry.m1 == m1 && entry.j2 == j2 && entry.m2 == m2 && entry.j == j && entry.m == m)
        {
            return &entry;
        }
    }
    return nullptr;
}

double ClebschGordanTable::coefficient(const int j1,
                                       const int m1,
                                       const int j2,
                                       const int m2,
                                       const int j,
                                       const int m) const
{
    const Entry* entry = find(j1, m1, j2, m2, j, m);
    return entry == nullptr ? 0.0 : entry->value;
}

double ClebschGordanTable::h_coefficient(const int j1,
                                         const int m1,
                                         const int mp1,
                                         const int j2,
                                         const int m2,
                                         const int mp2,
                                         const int j,
                                         const int m,
                                         const int mp) const
{
    return coefficient(j1, m1, j2, m2, j, m) * coefficient(j1, mp1, j2, mp2, j, mp);
}

} // namespace SNAP
} // namespace ModuleMD
