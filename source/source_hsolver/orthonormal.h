#ifndef HSOLVER_ORTHONORMAL_H
#define HSOLVER_ORTHONORMAL_H

#include "source_hsolver/linear_algebra.h"

#include <array>
#include <string>

namespace hsolver
{
enum class OrthMethod
{
    none,
    cholesky,
    lowdin,
    newton_schulz
};
OrthMethod parse_orth_method(const std::string& name);
const char* orth_method_name(OrthMethod method);

/** @brief Diagnostics for one collective correction; numerical failures are recoverable. */
struct OrthResult
{
    double before = 0.0;
    double after = 0.0;
    double seconds = 0.0;
    int passes = 0;
    int fallbacks = 0;
    bool skipped = false;
    bool invalid_input = false;
    OrthMethod actual = OrthMethod::none;
    std::string reason;
    // Failed methods, nonfinite Gram, rejected candidates, and last fallback methods.
    std::array<int, 8> events{};
    std::vector<double> norms;
};

/** @brief Maximum elementwise distance from the identity, including nonfinite detection. */
double orth_error(const std::vector<std::complex<double>>& gram, int bands);
/** @brief Construct a correction with bounded, explicitly reported fallback attempts. */
bool orth_transform(const std::vector<std::complex<double>>& gram,
                    int bands,
                    OrthMethod method,
                    std::vector<std::complex<double>>* transform,
                    OrthResult* result);

/** @brief Pool-local orthonormalization with FP64 products and native-precision updates. */
template <typename T, typename Device>
class Orthonormal
{
  private:
    const diag_comm_info comm_;
    LinearAlgebra<T, Device> algebra_;
    ct::Tensor candidate_;
    bool factor(const std::vector<std::complex<double>>& gram,
                int bands,
                OrthMethod method,
                std::vector<std::complex<double>>* transform,
                OrthResult* result);
    void rotate(const T* input, T* output, int ld, int dim, int bands, const std::vector<std::complex<double>>& transform);
    std::vector<std::complex<double>> gram(const T* input, int ld, int dim, int bands);

  public:
    explicit Orthonormal(const diag_comm_info& comm);
    /** @brief Preserve input on failed candidates; only valid rows may be modified. */
    OrthResult apply(T* input, int ld, int dim, int bands, OrthMethod method);
};
} // namespace hsolver
#endif
