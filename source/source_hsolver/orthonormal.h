#ifndef HSOLVER_ORTHONORMAL_H
#define HSOLVER_ORTHONORMAL_H

#include "source_hsolver/linear_algebra.h"

#include <array>
#include <limits>
#include <string>
#include <type_traits>

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

enum class OrthStatus
{
    accepted,
    unchanged,
    inspected, ///< Read-only Gram diagnostics, without enforcing the correction tolerance.
    disabled,
    failed
};

enum class OrthFailure
{
    none,
    nonfinite_gram,
    nonfinite_norm,
    nonpositive_norm,
    factorization_failed,
    invalid_candidate,
    no_improvement,
    tolerance_not_met
};
const char* orth_failure_name(OrthFailure failure);

/** @brief Acceptance tolerance for the native wavefunction precision. */
template <typename T>
constexpr double orth_tolerance()
{
    return std::is_same<T, std::complex<float>>::value ? 1e-6 : 1e-12;
}

/** @brief Collective outcome; failed states must not be used to construct a density. */
struct OrthResult
{
    double before = std::numeric_limits<double>::quiet_NaN();
    double after = std::numeric_limits<double>::quiet_NaN();
    bool gram_checked = false;
    double seconds = 0.0;
    int passes = 0;
    int fallbacks = 0;
    OrthStatus status = OrthStatus::failed;
    OrthFailure failure = OrthFailure::none;
    OrthMethod actual = OrthMethod::none;
    std::string reason;
    // Failed methods, nonfinite Gram, rejected candidates, and attempted fallback methods.
    std::array<int, 7> events{};
    std::vector<double> norms;
};

/** @brief Maximum elementwise distance from the identity, including nonfinite detection. */
double orth_error(const std::vector<std::complex<double>>& gram, int bands);
/** @brief Construct one correction; the caller validates candidates and chooses fallbacks. */
bool orth_transform(const std::vector<std::complex<double>>& gram,
                    int bands,
                    OrthMethod method,
                    std::vector<std::complex<double>>* transform);

/** @brief Pool-local orthonormalization with FP64 products and native-precision updates. */
template <typename T, typename Device>
class Orthonormal
{
  private:
    const diag_comm_info comm_;
    LinearAlgebra<T, Device> algebra_;
    ct::Tensor candidate_;
    bool factor(const std::vector<std::complex<double>>& gram, int bands, OrthMethod method, std::vector<std::complex<double>>* transform);
    bool try_candidate(T* input,
                       int ld,
                       int dim,
                       int bands,
                       const std::vector<std::complex<double>>& transform,
                       std::vector<std::complex<double>>* gram,
                       OrthResult* result);
    void correct(T* input, int ld, int dim, int bands, OrthMethod method, std::vector<std::complex<double>>* gram, OrthResult* result);
    void rotate(const T* input, T* output, int ld, int dim, int bands, const std::vector<std::complex<double>>& transform);
    std::vector<std::complex<double>> gram(const T* input, int ld, int dim, int bands);

  public:
    explicit Orthonormal(const diag_comm_info& comm);
    /** @brief Inspect norms and optionally overlaps without changing orbitals or enforcing a correction tolerance. */
    OrthResult inspect(const T* input, int ld, int dim, int bands, bool full_gram);
    /** @brief Preserve input on rejected candidates; a failed overall result is not usable. */
    OrthResult apply(T* input, int ld, int dim, int bands, OrthMethod method);
};
} // namespace hsolver
#endif
