#ifndef MODULE_HAMILT_GENERAL_MODULE_XC_KERNELS_H_
#define MODULE_HAMILT_GENERAL_MODULE_XC_KERNELS_H_

#include <complex>
#include <module_psi/psi.h>
#include <module_base/macros.h>

namespace hamilt {

#if defined(__CUDA)
// Device arrays use the full FFT-grid layout or the compact plane-wave layout.
void xc_gradient_coefficients_gpu(int npw, int component, const int* grid_index,
                                 const double* gcar, const std::complex<double>* rhog,
                                 std::complex<double>* grid);
void xc_extract_real_gpu(int size, int stride, int offset, double scale,
                         const std::complex<double>* grid, double* result);
void xc_pack_component_gpu(int size, int component, const double* vectors,
                           std::complex<double>* grid);
void xc_accumulate_divergence_gpu(int npw, int nxyz, int component,
                                const int* grid_index, const double* gcar,
                                const std::complex<double>* grid,
                                std::complex<double>* reciprocal);

// Noncollinear fields have shape (4, size); channel gradients have shape (2, size, 3).
// Reduction storage contains four rows of xc_reduction_size_gpu(size) block sums.
int xc_reduction_size_gpu(int size);
void xc_pack_density_gpu(int size, const double* density, std::complex<double>* grid);
void xc_gather_density_gpu(int npw, int nxyz, const int* index,
                           const std::complex<double>* grid, std::complex<double>* reciprocal);
void xc_add_reciprocal_gpu(int npw, const std::complex<double>* source,
                          std::complex<double>* destination);
void xc_spin_gradient_gpu(int size, int spin, int component, double scale,
                         const std::complex<double>* grid, const double* rho, double* gradient);
void xc_pbe_spin_gpu(int size, const double* rho, const double* core, double* gradient_flux,
                    double* potential, double* sums);
void xc_pack_spin_flux_gpu(int size, int spin, int component, const double* rho,
                          const double* flux, std::complex<double>* grid);
void xc_apply_spin_divergence_gpu(int size, int spin, double scale,
                                const std::complex<double>* grid, const double* rho,
                                double* potential, double* magnetic_divergence, double* sums);
void xc_finish_spin_potential_gpu(int size, const double* rho, const double* magnetic_divergence,
                                double* potential, double* sums);
#endif

template <typename T, typename Device>
struct xc_functional_grad_wfc_op {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const int& ik,
        const int& pol,
        const int& npw,
        const int& npwx,
        const Real& tpiba,
        const Real * gcar,
        const Real * kvec_c,
        const T * rhog,
        T* porter);
    
    void operator()(
        const int& ipol,
        const int& nrxx,
        const T * porter,
        T* grad);
};

} // namespace hamilt
#endif // MODULE_HAMILT_GENERAL_MODULE_XC_KERNELS_H_
