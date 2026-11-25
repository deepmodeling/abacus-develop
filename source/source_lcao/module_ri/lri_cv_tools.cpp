#include "lri_cv_tools.h"

#include <complex>

template <typename T>
std::pair<LRI_CV_Tools::RIBlockMatrix<T>, std::array<LRI_CV_Tools::RIBlockMatrix<T>, 3>>
LRI_CV_Tools::cal_I_dI(const RIBlockMatrix<T>& metric,
                      const std::array<RI::Tensor<T>, 3>& derivative01,
                      const std::array<RI::Tensor<T>, 3>& derivative10,
                      const double threshold)
{
    const std::vector<std::size_t> sizes = {metric[0][0].shape[0], metric[1][1].shape[0]};
    std::array<RI::Tensor<T>, 3> derivatives;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        const RIBlockMatrix<T> blocks = {
            {RI::Tensor<T>({sizes[0], sizes[0]}), derivative01[axis]},
            {derivative10[axis], RI::Tensor<T>({sizes[1], sizes[1]})}};
        Inverse_Matrix<T> pack;
        pack.input(blocks);
        derivatives[axis] = pack.output();
    }
    Inverse_Matrix<T> inverse;
    inverse.input(metric);
    const auto responses = inverse.cal_inverse_with_derivatives(derivatives, threshold);
    std::array<RIBlockMatrix<T>, 3> block_responses;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        Inverse_Matrix<T> pack;
        pack.input(responses[axis]);
        block_responses[axis] = pack.output(sizes, sizes);
    }
    auto inverse_blocks = inverse.output(sizes, sizes);
    return std::make_pair(inverse_blocks, block_responses);
}

template std::pair<LRI_CV_Tools::RIBlockMatrix<double>,
                   std::array<LRI_CV_Tools::RIBlockMatrix<double>, 3>>
LRI_CV_Tools::cal_I_dI(
    const RIBlockMatrix<double>& metric,
    const std::array<RI::Tensor<double>, 3>& derivative01,
    const std::array<RI::Tensor<double>, 3>& derivative10,
    const double threshold);

template std::pair<LRI_CV_Tools::RIBlockMatrix<std::complex<double>>,
                   std::array<LRI_CV_Tools::RIBlockMatrix<std::complex<double>>, 3>>
LRI_CV_Tools::cal_I_dI(
    const RIBlockMatrix<std::complex<double>>& metric,
    const std::array<RI::Tensor<std::complex<double>>, 3>& derivative01,
    const std::array<RI::Tensor<std::complex<double>>, 3>& derivative10,
    const double threshold);
