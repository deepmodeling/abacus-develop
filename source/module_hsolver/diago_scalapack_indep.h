#ifndef DIAGO_SCALAPACK_H
#define DIAGO_SCALAPACK_H

#include <complex>
#include <utility>
#include <vector>

#include "diagh.h"

namespace hsolver
{
template <typename T = std::complex<double>>
class DiagoScalapack : public DiagH<T>
{
    private:
        using Real = typename GetTypeReal<T>::type;
        const int dim;
        const double diag_thr_in;
        const int iter_nmax;
    public:
        DiagoScalapack(const int& nband_in,
                        const int& nbasis_in,
                        const double& diag_thr_in,
                        const int& diag_nmax_in);
        ~DiagoScalapack();
        void diag(T* h_mat, T* s_mat, const int* const desc, T* psi, Real* eigenvalue_in) override;

    private:
        void pdsygvx_diag(const int *const desc,
                        const int ncol,
                        const int nrow,
                        const double *const h_mat,
                        const double *const s_mat,
                        double *const ekb,
                        double *wfc_2d);
        void pzhegvx_diag(const int *const desc,
                        const int ncol,
                        const int nrow,
                        const std::complex<double> *const h_mat,
                        const std::complex<double> *const s_mat,
                        double *const ekb,
                        std::complex<double> *wfc_2d);

        std::pair<int, std::vector<int>> pdsygvx_once(const int *const desc,
                                                    const int ncol,
                                                    const int nrow,
                                                    const double *const h_mat,
                                                    const double *const s_mat,
                                                    double *const ekb,
                                                    double *wfc_2d) const;
        std::pair<int, std::vector<int>> pzhegvx_once(const int *const desc,
                                                    const int ncol,
                                                    const int nrow,
                                                    const std::complex<double> *const h_mat,
                                                    const std::complex<double> *const s_mat,
                                                    double *const ekb,
                                                    std::complex<double> *wfc_2d) const;

        int degeneracy_max = 12; // For reorthogonalized memory. 12 followes siesta.

        void post_processing(const int info, const std::vector<int> &vec);
};

} // namespace hsolver

#endif