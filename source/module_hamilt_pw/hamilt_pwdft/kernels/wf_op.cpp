#include "module_hamilt_pw/hamilt_pwdft/kernels/wf_op.h"
#include "module_base/libm/libm.h"

namespace hamilt{

template <typename FPTYPE>
struct cal_sk_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(const psi::DEVICE_CPU* ctx,
                    const int& ik,
                    const int& ntype,
                    const int& nx,
                    const int& ny,
                    const int& nz,
                    const int& rho_nx,
                    const int& rho_ny,
                    const int& rho_nz,
                    const int& npw,
                    const int& npwx,
                    const int& fftny,
                    const int& eigts1_nc,
                    const int& eigts2_nc,
                    const int& eigts3_nc,
                    const int* atom_na,
                    const int* igl2isz,
                    const int* is2fftixy,
                    const FPTYPE& TWO_PI,
                    const FPTYPE* kvec_c,
                    const FPTYPE* atom_tau,
                    std::complex<FPTYPE>* eigts1,
                    std::complex<FPTYPE>* eigts2,
                    std::complex<FPTYPE>* eigts3,
                    std::complex<FPTYPE>* sk)
    {
#ifdef _OPENMP
#pragma omp parallel
{
#endif
        int iat = 0;
        for (int it = 0; it < ntype; it++) {
            for (int ia = 0; ia < atom_na[it]; ia++) {
                FPTYPE arg = 0.0;
                for (int ii = 0; ii < 3; ii++) {
                    arg += kvec_c[ik * 3 + ii] * atom_tau[iat * 3 + ii];
                }
                arg *= TWO_PI;
                FPTYPE sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                const std::complex<FPTYPE> kphase = std::complex<FPTYPE>(cosp, -sinp);
#ifdef _OPENMP
#pragma omp for
#endif
                for (int igl = 0; igl < npw; ++igl) {
                    const int isz = igl2isz[ik * npwx + igl];
                    int iz = isz % nz;
                    const int is = isz / nz;
                    const int ixy = is2fftixy[is];
                    int ix = ixy / fftny;
                    int iy = ixy % fftny;
                    if (ix >= int(nx / 2) + 1)
                        ix -= nx;
                    if (iy >= int(ny / 2) + 1)
                        iy -= ny;
                    if (iz >= int(nz / 2) + 1)
                        iz -= nz;
                    ix += rho_nx;
                    iy += rho_ny;
                    iz += rho_nz;
                    sk[iat * npw + igl] = kphase * eigts1[iat * eigts1_nc + ix] * eigts2[iat * eigts2_nc + iy]
                                          * eigts3[iat * eigts3_nc + iz];
                }
                iat++;
            }
        }
#ifdef _OPENMP
}
#endif
    }
};

template struct cal_sk_op<float, psi::DEVICE_CPU>;
template struct cal_sk_op<double, psi::DEVICE_CPU>;

}  // namespace hamilt

