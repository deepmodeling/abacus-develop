#ifndef STERNHEIMER_SIAB_DATA_H
#define STERNHEIMER_SIAB_DATA_H

#include <complex>
#include <string>
#include <vector>

namespace module_ri
{
namespace sternheimer_siab
{

struct PrimitiveBlock
{
    std::string element;
    int atom_index;
    int l;
    int m;
    int n_primitive;
    int offset;
};

struct ReferenceRow
{
    int occupied_state;
    int auxiliary_channel;
    double frequency_ha;
    double occupation;
    double frequency_weight;
    double norm;
    std::vector<std::complex<double>> q;
};

} // namespace sternheimer_siab
} // namespace module_ri

#endif
