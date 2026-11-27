#ifndef ABACUS_LR_WINDOW_RESPONSE_H
#define ABACUS_LR_WINDOW_RESPONSE_H
#include <ATen/core/tensor.h>
#include <string>
#include <vector>
namespace LR
{
template <typename T> struct GradientInputs;
// Owned corrections for one batch of amplitudes. Empty means a full window.
template <typename T>
struct WindowResponse
{
    // [state][spin * nk + k]; one spin block for spin-adapted amplitudes.
    std::vector<std::vector<container::Tensor>> density;
    // Direct sym(C_r epsilon_a W_ra C_a^T) only; the occupied contribution
    // is obtained from Hgs[T + Dwin + DZ]/2 in the existing Wc calculation.
    std::vector<std::vector<container::Tensor>> direct_edm;
};
template <typename T>
WindowResponse<T> build_window_response(const GradientInputs<T>& inputs,
    const container::Tensor& X, int nstates, int nactive, int ispin,
    const std::string& in_dir);
template <typename T>
WindowResponse<T> build_open_shell_window_response(const GradientInputs<T>& inputs,
    const container::Tensor& X, int nstates, const std::vector<int>& nactive);
}
#endif
