#include "module_hamilt_pw/hamilt_pwdft/kernels/vnl_op.h"

namespace hamilt {

template <typename FPTYPE>
static inline FPTYPE _polynomial_interpolation(
        const FPTYPE *table,
        const int &dim1,
        const int &dim2,
        const int &tab_2,
        const int &tab_3,
        const int &table_length,
        const FPTYPE &table_interval,
        const FPTYPE &x)
{
    const FPTYPE position = x / table_interval;
    const int iq = static_cast<int>(position);

    const FPTYPE x0 = position - static_cast<FPTYPE>(iq);
    const FPTYPE x1 = 1.0 - x0;
    const FPTYPE x2 = 2.0 - x0;
    const FPTYPE x3 = 3.0 - x0;
    const FPTYPE y =
            table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0] * x1 * x2 * x3 / 6.0 +
            table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 1] * x0 * x2 * x3 / 2.0 -
            table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 2] * x1 * x0 * x3 / 2.0 +
            table[(dim1 * tab_2 + dim2) * tab_3 + iq + 0 + 3] * x1 * x2 * x0 / 6.0 ;

//	ModuleBase::timer::tick("PolyInt","Poly_Interpo_2");
    return y;
}

template <typename FPTYPE>
struct cal_vnl_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
        const psi::DEVICE_CPU *ctx,
        const int &ntype,
        const int &npw,
        const int &npwx,
        const int &nhm,
        const int &NQX,
        const int &tab_2,
        const int &tab_3,
        const int * atom_na,
        const int * atom_nb,
        const int * atom_nh,
        const FPTYPE &DQ,
        const FPTYPE &tpiba,
        const std::complex<FPTYPE> &NEG_IMAG_UNIT,
        const FPTYPE *gk,
        const FPTYPE *ylm,
        const FPTYPE *indv,
        const FPTYPE *nhtol,
        const FPTYPE *nhtolm,
        const FPTYPE *tab,
        FPTYPE *vkb1,
        const std::complex<FPTYPE> *sk,
        std::complex<FPTYPE> *vkb_in)
    {
        const int imag_pow_period = 4;
        // result table of pow(0-1i, int)
        static const std::complex<FPTYPE> pref_tab[imag_pow_period] = {{1, 0}, {0, -1}, {-1, 0}, {0, 1}};
#ifdef _OPENMP
#pragma omp parallel
{
#endif
        int jkb = 0, iat = 0;
        FPTYPE vq = 0.0;
        for (int it = 0; it < ntype; it++) {
            // calculate beta in G-space using an interpolation table
            const int nh = atom_nh[it];
            const int nbeta = atom_nb[it];

            for (int nb = 0; nb < nbeta; nb++) {
#ifdef _OPENMP
                #pragma omp for
#endif
                for (int ig = 0; ig < npw; ig++) {
                    const FPTYPE gnorm = sqrt(gk[ig * 3 + 0] * gk[ig * 3 + 0] + gk[ig * 3 + 1] * gk[ig * 3 + 1] +
                                              gk[ig * 3 + 2] * gk[ig * 3 + 2]) * tpiba;

                    vq = _polynomial_interpolation(
                            tab, it, nb, tab_2, tab_3, NQX, DQ, gnorm);

                    // add spherical harmonic part
                    for (int ih = 0; ih < nh; ih++) {
                        if (nb == indv[it * nhm + ih]) {
                            const int lm = static_cast<int>(nhtolm[it * nhm + ih]);
                            vkb1[ih * npw + ig] = ylm[lm * npw + ig] * vq;
                        }
                    } // end ih
                }
            } // end nbeta

            // vkb1 contains all betas including angular part for type nt
            // now add the structure factor and factor (-i)^l
            for (int ia = 0; ia < atom_na[it]; ia++) {
                for (int ih = 0; ih < nh; ih++) {
                    // std::complex<FPTYPE> pref = pow(NEG_IMAG_UNIT, nhtol[it * nhm + ih]);    //?
                    std::complex<FPTYPE> pref = pref_tab[int(nhtol[it * nhm + ih]) % imag_pow_period];
                    std::complex<FPTYPE> *pvkb = vkb_in + jkb * npwx;
#ifdef _OPENMP
                    #pragma omp for
#endif
                    for (int ig = 0; ig < npw; ig++) {
                        pvkb[ig] = vkb1[ih * npw + ig] * sk[iat * npw + ig] * pref;
                    }
                    ++jkb;
                } // end ih
                iat++;
            } // end ia
        } // enddo
#ifdef _OPENMP
}
#endif
    }
};

template struct cal_vnl_op<float, psi::DEVICE_CPU>;
template struct cal_vnl_op<double, psi::DEVICE_CPU>;

}  // namespace hamilt

