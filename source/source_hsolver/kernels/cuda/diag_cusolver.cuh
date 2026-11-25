#ifndef DIAG_CUSOLVER_CUH
#define DIAG_CUSOLVER_CUH

#include <complex>

#include <cusolverDn.h>

/// Dense generalized type-1 eigensolver on one CUDA device.
/// H and S are column-major with an authoritative upper triangle.
class Diag_Cusolver_gvd
{
  public:
    Diag_Cusolver_gvd();
    ~Diag_Cusolver_gvd();

    Diag_Cusolver_gvd(const Diag_Cusolver_gvd &) = delete;
    Diag_Cusolver_gvd &operator=(const Diag_Cusolver_gvd &) = delete;

    static bool use_subset(int n, int nbands);

    // Return false on a recoverable subset failure without updating W or V.
    bool solve_subset(
        int n,
        int nbands,
        const double *h,
        const double *s,
        double *w,
        double *v);
    // Solve the full spectrum and return the first nbands eigenpairs.
    void solve_full(
        int n,
        int nbands,
        const double *h,
        const double *s,
        double *w,
        double *v);
    void solve_full(
        int n,
        int nbands,
        const std::complex<double> *h,
        const std::complex<double> *s,
        double *w,
        std::complex<double> *v);

  private:
    cusolverDnHandle_t handle_ = nullptr;
};

#endif
