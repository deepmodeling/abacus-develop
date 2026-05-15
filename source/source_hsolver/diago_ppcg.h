#ifndef DIAGO_PPCG_H_
#define DIAGO_PPCG_H_

#include "source_base/macros.h"
#include "source_base/module_device/types.h"

#include <complex>
#include <functional>
#include <vector>

namespace hsolver
{

template <typename T = std::complex<double>, typename Device = base_device::DEVICE_CPU>
class DiagoPPCG
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    using HPsiFunc = std::function<void(T*, T*, const int, const int)>;

    explicit DiagoPPCG(const Real* precondition_in);

    void init_iter(const int nband, const int nband_l, const int nbasis, const int ndim);

    int diag(const HPsiFunc& hpsi_func,
             T* psi_in,
             Real* eigenvalue_in,
             const std::vector<double>& ethr_band);

  private:
    int n_band = 0;
    int n_band_l = 0;
    int n_basis = 0;
    int n_dim = 0;

    const Real* precondition = nullptr;

    std::vector<T> hpsi;
    std::vector<T> w;
    std::vector<T> hw;
    std::vector<T> p;
    std::vector<T> hp;
    std::vector<T> p_new;
    std::vector<T> hp_new;
    std::vector<T> hpsi_new;
    std::vector<T> work;
    std::vector<Real> eigen;
    std::vector<Real> err;

    T inner_product(const T* lhs, const T* rhs) const;
    Real vector_norm(const T* vec) const;
    void scale_vector(T* vec, const Real alpha) const;
    void axpy_vector(T* y, const T* x, const T alpha) const;
    void copy_vector(T* dst, const T* src) const;
    void zero_vector(T* vec) const;

    bool test_error(const std::vector<double>& ethr_band) const;
    void calc_hpsi(const HPsiFunc& hpsi_func, T* psi_in, std::vector<T>& hpsi_out) const;
    void modified_gram_schmidt(T* psi_in, std::vector<T>& hpsi_in) const;
    void rotate_block(T* block, const std::vector<T>& coeff, std::vector<T>& workspace) const;
    void rayleigh_ritz(T* psi_in, std::vector<T>& hpsi_in);
    void calc_preconditioned_residual(T* psi_in);
    void project_to_orthogonal_complement(T* psi_in, std::vector<T>& block) const;
    bool solve_small_problem(const int active_dim, T* hsmall, T* ssmall, T* coeff, Real* eval) const;
    void update_vectors_from_ppcg_subspace(T* psi_in);
};

} // namespace hsolver

#endif
