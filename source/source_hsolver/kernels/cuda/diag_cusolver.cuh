#ifndef DIAG_CUSOLVER_CUH
#define DIAG_CUSOLVER_CUH

#include <complex>
#include <string>

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
    // retain_vectors keeps the device eigenvectors alive for subsequent GPU DM.
    bool solve_subset(
        int n,
        int nbands,
        const double *h,
        const double *s,
        double *w,
        double *v,
        bool retain_vectors);
    // Solve the full spectrum and return the first nbands eigenpairs.
    void solve_full(
        int n,
        int nbands,
        const double *h,
        const double *s,
        double *w,
        double *v,
        bool retain_vectors);
    void solve_full(
        int n,
        int nbands,
        const std::complex<double> *h,
        const std::complex<double> *s,
        double *w,
        std::complex<double> *v);

    // Construct DM from the last successful real full or subset solve.
    // On failure, report the failing stage; the caller must terminate the solve.
    bool density_from_retained(int bands, const double* weights, double* density, std::string& error);

  private:
    void release_eigenvectors();
    double* eigenvectors_ = nullptr;
    double* weighted_vectors_ = nullptr;
    double* density_workspace_ = nullptr;
    double* occupations_ = nullptr;
    int rows_ = 0;
    int bands_ = 0;
    cusolverDnHandle_t handle_ = nullptr;
};

#endif
