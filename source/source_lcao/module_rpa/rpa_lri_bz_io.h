#ifndef RPA_LRI_BZ_IO_H
#define RPA_LRI_BZ_IO_H

#include <array>
#include <string>
#include <vector>

namespace RpaLriDetail
{
struct BzSamplePoint
{
    double weight = 0.0;
    std::array<double, 3> direct = {{0.0, 0.0, 0.0}};
    std::array<double, 3> cartesian = {{0.0, 0.0, 0.0}};
};

struct BzSampleData
{
    std::array<int, 3> mesh = {{0, 0, 0}};
    std::vector<BzSamplePoint> reduced;
    std::vector<std::array<double, 3>> full_cartesian;
    std::vector<int> full_map;
};

void write_bz_sample(const std::string& filename, const BzSampleData& data);
} // namespace RpaLriDetail

#endif
