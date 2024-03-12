#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hamilt{

template <typename FPTYPE>
struct cal_vkb1_nl_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
        const psi::DEVICE_CPU *ctx,
        const int &nkb,
        const int &npwx,
        const int &npwk_max,
        const int &vkb_nc,
        const int &nbasis,
        const int &ik,
        const int &ipol,
        const std::complex<FPTYPE> &NEG_IMAG_UNIT,
        const std::complex<FPTYPE> *vkb,
        const FPTYPE *gcar,
        std::complex<FPTYPE> *vkb1)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int i = 0; i < nkb; i++) {
            for (int ig = 0; ig < nbasis; ig++) {
                std::complex<FPTYPE> *pvkb1 = vkb1 + i * npwx;
                const std::complex<FPTYPE> *pvkb = vkb + i * vkb_nc;
                pvkb1[ig] = pvkb[ig] * NEG_IMAG_UNIT * gcar[(ik * npwk_max + ig) * 3 + ipol];
            }
        }
    }
};

template <typename FPTYPE>
struct cal_force_nl_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(const psi::DEVICE_CPU* ctx,
                    const bool& nondiagonal,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force)
    {
#ifdef _OPENMP
#pragma omp parallel
{
#endif
        int iat0 = 0;
        int sum0 = 0;
        for (int it = 0; it < ntype; it++) {
            const int Nprojs = atom_nh[it];
#ifdef _OPENMP
#pragma omp for collapse(2)
#endif
            for (int ia = 0; ia < atom_na[it]; ia++) {
                for (int ib = 0; ib < nbands_occ; ib++) {
                    FPTYPE local_force[3] = {0, 0, 0};
                    FPTYPE fac = d_wg[ik * wg_nc + ib] * 2.0 * tpiba;
                    FPTYPE ekb_now = d_ekb[ik * wg_nc + ib];
                    int iat = iat0 + ia;
                    int sum = sum0 + ia * Nprojs;
                    for (int ip = 0; ip < Nprojs; ip++) {
                        // Effective values of the D-eS coefficients
                        FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip]
                                    - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip];
                        const int inkb = sum + ip;
                        //out<<"\n ps = "<<ps;

                        for (int ipol=0; ipol<3; ipol++) {
                            const FPTYPE dbb = (conj(dbecp[ipol * nbands * nkb + ib * nkb + inkb]) * becp[ib * nkb + inkb]).real();
                            local_force[ipol] -= ps * fac * dbb;
                            //cf[iat*3+ipol] += ps * fac * dbb;
                        }
                        if (nondiagonal)
                        {
                            for (int ip2=0; ip2<Nprojs; ip2++) {
                                if ( ip != ip2 ) {
                                    const int jnkb = sum + ip2;
                                    FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2]
                                                - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip2];

                                    for (int ipol = 0; ipol < 3; ipol++) {
                                        const FPTYPE dbb = ( conj( dbecp[ipol * nbands * nkb + ib * nkb + inkb] ) * becp[ib * nkb + jnkb] ).real();
                                        local_force[ipol] -= ps * fac * dbb;
                                    }
                                }
                            }
                        }
                    }
#ifdef _OPENMP
                    if (omp_get_num_threads() > 1)
                    {
                        for (int ipol=0; ipol<3; ipol++) {
                            #pragma omp atomic
                            force[iat * forcenl_nc + ipol] += local_force[ipol];
                        }
                    }
                    else
#endif
                    {
                        for (int ipol=0; ipol<3; ipol++) {
                            force[iat * forcenl_nc + ipol] += local_force[ipol];
                        }
                    }
                }
            } // end ia
            iat0 += atom_na[it];
            sum0 += atom_na[it] * Nprojs;
        } //end it
#ifdef _OPENMP
}
#endif
    }
};

template struct cal_vkb1_nl_op<float, psi::DEVICE_CPU>;
template struct cal_force_nl_op<float, psi::DEVICE_CPU>;

template struct cal_vkb1_nl_op<double, psi::DEVICE_CPU>;
template struct cal_force_nl_op<double, psi::DEVICE_CPU>;

}  // namespace hamilt

