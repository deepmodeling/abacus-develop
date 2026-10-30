#ifndef HSOLVER_FIXED_OVERLAP_H
#define HSOLVER_FIXED_OVERLAP_H

#include <complex>

namespace hsolver
{

/// Cache an orthonormal representation of fixed Hermitian subspaces on CUDA.
/// For S = U^H U, the cached matrices are U^-H H U^-1 and R U^-1.
/// All methods use the existing single-host-thread solver and BLAS context.
/// Inputs remain unchanged; allocation failure leaves an unavailable cache.
class FixedOverlapSubspace
{
  public:
    using Complex = std::complex<double>;

    FixedOverlapSubspace(int points, int bands, int projector_rows);
    ~FixedOverlapSubspace();
    FixedOverlapSubspace(const FixedOverlapSubspace&) = delete;
    FixedOverlapSubspace& operator=(const FixedOverlapSubspace&) = delete;

    bool available() const { return storage_ != nullptr; }

    /// Prepare a k point from compact column-major device matrices.
    /// H and S use their upper triangles, and S must be positive definite.
    void prepare(int point, const Complex* h, const Complex* s, const Complex* projectors);
    const Complex* hamiltonian(int point) const;
    const Complex* projectors(int point) const;

    /// Solve a field-shifted orthonormal problem and rotate its projectors.
    /// Eigenvalues are written to host memory; rotated projectors stay on CUDA.
    void diagonalize(const Complex* h, const Complex* projectors, Complex* rotated, double* eigenvalues);

  private:
    const int points_, bands_, rows_;
    Complex* storage_ = nullptr;
    Complex* h_ = nullptr;
    Complex* r_ = nullptr;
    Complex* cholesky_ = nullptr;
    Complex* vectors_ = nullptr;
    double* eigenvalues_ = nullptr;
    int* info_ = nullptr;
    int factor_work_ = 0;
    int eigen_work_ = 0;

    void whiten(Complex* h, Complex* projectors) const;
    void require_success(const char* operation) const;
};

} // namespace hsolver

#endif
