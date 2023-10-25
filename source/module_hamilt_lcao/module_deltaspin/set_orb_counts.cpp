#include "spin_constrain.h"

/// @brief  set orbital counts
template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::set_orb_counts(std::map<int, int> atomCounts_in,
                                                                          std::map<int, int> orbitalCounts_in)
{
    this->clear_atomCounts();
    this->clear_orbitalCounts();
    this->set_atomCounts(atomCounts_in);
    this->set_orbitalCounts(orbitalCounts_in);
}