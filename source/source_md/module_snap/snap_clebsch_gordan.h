#ifndef SNAP_CLEBSCH_GORDAN_H
#define SNAP_CLEBSCH_GORDAN_H

#include <vector>

namespace ModuleMD
{
namespace SNAP
{

/**
 * @brief Clebsch-Gordan coefficients indexed with doubled angular momenta.
 *
 * If the stored index is j, the physical angular momentum is j/2. The same
 * convention is used for m and m'.
 */
class ClebschGordanTable
{
  public:
    explicit ClebschGordanTable(int twojmax);

    double coefficient(int j1, int m1, int j2, int m2, int j, int m) const;

    double h_coefficient(int j1,
                         int m1,
                         int mp1,
                         int j2,
                         int m2,
                         int mp2,
                         int j,
                         int m,
                         int mp) const;

  private:
    struct Entry
    {
        int j1;
        int m1;
        int j2;
        int m2;
        int j;
        int m;
        double value;
    };

    static bool allowed(int j1, int m1, int j2, int m2, int j, int m);
    static double coefficient_uncached(int j1, int m1, int j2, int m2, int j, int m);
    static double factorial(int value);
    const Entry* find(int j1, int m1, int j2, int m2, int j, int m) const;

    int twojmax_;
    std::vector<Entry> entries_;
};

} // namespace SNAP
} // namespace ModuleMD

#endif // SNAP_CLEBSCH_GORDAN_H
