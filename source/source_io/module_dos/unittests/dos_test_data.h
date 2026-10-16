#ifndef DOS_TEST_DATA_H
#define DOS_TEST_DATA_H

// Small synthetic test data for module_dos unit tests:
// 4 k-points with equal weight (sum of wk = 2 for spin-unpolarized),
// 4 bands. Energies are in eV.

#include <vector>

namespace dos_test_data
{
const int nks = 4;
const int nkstot = 4;
const int nbands = 4;

// k-point weights (already doubled for spin degeneracy)
const std::vector<double> wk = {0.5, 0.5, 0.5, 0.5};

// band energies in eV, row-major [nks][nbands]
const std::vector<double> ekb_ev = {
    -5.0, 1.0, 3.0, 6.0,
    -4.5, 1.5, 3.5, 6.5,
    -4.0, 2.0, 4.0, 7.0,
    -3.5, 2.5, 4.5, 7.5,
};

// band occupation numbers (weight included), row-major [nks][nbands]
const std::vector<double> wg = {
    0.5, 0.0, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0,
};
} // namespace dos_test_data

#endif
