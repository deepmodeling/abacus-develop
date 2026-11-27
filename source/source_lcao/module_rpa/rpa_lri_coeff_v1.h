#ifndef RPA_LRI_COEFF_V1_H
#define RPA_LRI_COEFF_V1_H

#include <cstdint>
#include <string>
#include <vector>

namespace RI
{
template <typename T>
class Tensor;
}

namespace RpaLriDetail
{
struct LriCoefV1Block
{
    std::int32_t ia1 = 0;
    std::int32_t ia2 = 0;
    std::int32_t r[3] = {0, 0, 0};
    double max_abs = 0.0;
    std::int64_t offset = 0;
    int nw1 = 0;
    int nw2 = 0;
    int naux = 0;
    const RI::Tensor<double>* tensor = nullptr;
};

void write_lri_coef_v1(const std::string& filename,
                       const std::int32_t natom,
                       std::vector<LriCoefV1Block>& blocks);
} // namespace RpaLriDetail

#endif
