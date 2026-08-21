#ifndef ABACUS_D3_DATA_H
#define ABACUS_D3_DATA_H

namespace vdw
{
namespace d3
{
namespace data
{

constexpr int max_element = 103;
constexpr int max_reference = 7;

int reference_count(int atomic_number);
double reference_cn(int atomic_number, int reference);
double reference_c6(int atomic_number_i, int reference_i, int atomic_number_j, int reference_j);
double covalent_radius(int atomic_number);
double r4r2(int atomic_number);
double vdw_radius(int atomic_number_i, int atomic_number_j);

} // namespace data
} // namespace d3
} // namespace vdw

#endif // ABACUS_D3_DATA_H
