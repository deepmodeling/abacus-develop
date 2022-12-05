#include "src_pw/include/force_multi_device.h"
#include <iostream>

namespace src_pw{

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
        for (int i = 0; i < nkb; i++) {
            std::complex<FPTYPE> *pvkb1 = vkb1 + i * npwx;
            const std::complex<FPTYPE> *pvkb = vkb + i * vkb_nc;
            for (int ig = 0; ig < nbasis; ig++) {
                pvkb1[ig] = pvkb[ig] * NEG_IMAG_UNIT * gcar[(ik * npwk_max + ig) * 3 + ipol];
            }
        }
    }
};

template <typename FPTYPE>
struct cal_force_nl_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
        const psi::DEVICE_CPU *ctx,
        const bool &multi_proj,
        const int &nbands_occ,
        const int &wg_nc,
        const int &ntype,
        const int &spin,
        const int &deeq_2,
        const int &deeq_3,
        const int &deeq_4,
        const int &forcenl_nc,
        const int &nbands,
        const int &ik,
        const int &nkb,
        const int *atom_nh,
        const int *atom_na,
        const FPTYPE &tpiba,
        const FPTYPE *d_wg,
        const FPTYPE *deeq,
        const std::complex<FPTYPE> *becp,
        const std::complex<FPTYPE> *dbecp,
        FPTYPE *force)
    {
        for (int ib = 0; ib < nbands_occ; ib++) {
            FPTYPE fac = d_wg[ik * wg_nc + ib] * 2.0 * tpiba;
            int iat = 0;
            int sum = 0;
            for (int it = 0; it < ntype; it++) {
                const int Nprojs = atom_nh[it];
                for (int ia = 0; ia < atom_na[it]; ia++) {
                    for (int ip = 0; ip < Nprojs; ip++) {
                        // FPTYPE ps = GlobalC::ppcell.deeq[GlobalV::CURRENT_SPIN, iat, ip, ip];
                        FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip];
                        const int inkb = sum + ip;
                        //out<<"\n ps = "<<ps;

                        for (int ipol=0; ipol<3; ipol++) {
                            const FPTYPE dbb = (conj(dbecp[ipol * nbands * nkb + ib * nkb + inkb]) * becp[ib * nkb + inkb]).real();
                            force[iat * forcenl_nc + ipol] -= ps * fac * dbb;
                            //cf[iat*3+ipol] += ps * fac * dbb;
                        }
                    }

                    if(multi_proj) {
                        for (int ip = 0; ip < Nprojs; ip++) {
                            const int inkb = sum + ip;
                            //for (int ip2=0; ip2<Nprojs; ip2++)
                            for (int ip2=0; ip2<Nprojs; ip2++) {
                                if ( ip != ip2 ) {
                                    const int jnkb = sum + ip2;
                                    FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] ;

                                    for (int ipol = 0; ipol < 3; ipol++) {
                                        const FPTYPE dbb = ( conj( dbecp[ipol * nbands * nkb + ib * nkb + inkb] ) * becp[ib * nkb + jnkb] ).real();
                                        force[iat * forcenl_nc + ipol] -= ps * fac * dbb;
                                    }
                                }
                            }
                        }
                    }
                    ++iat;
                    sum+=Nprojs;
                }
            } //end it
        }
    }
};

template struct cal_vkb1_nl_op<double, psi::DEVICE_CPU>;
template struct cal_force_nl_op<double, psi::DEVICE_CPU>;

}  // namespace src_pw

