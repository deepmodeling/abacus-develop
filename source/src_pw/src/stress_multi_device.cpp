#include "src_pw/include/stress_multi_device.h"
#include <iostream>
#include <iomanip>

namespace src_pw{

    template <typename FPTYPE>
    struct cal_dbecp_noevc_nl_op<FPTYPE, psi::DEVICE_CPU> {
        void operator()(
                const psi::DEVICE_CPU *ctx,
                const int &ipol,
                const int &jpol,
                const int &nkb,
                const int &npw,
                const int &npwx,
                const int &ik,
                const FPTYPE &tpiba,
                const FPTYPE *gcar,
                const FPTYPE *kvec_c,
                std::complex<FPTYPE> *vkbi,
                std::complex<FPTYPE> *vkbj,
                std::complex<FPTYPE> *vkb,
                std::complex<FPTYPE> *vkb1,
                std::complex<FPTYPE> *vkb2,
                std::complex<FPTYPE> *dbecp_noevc)
        {
            for (int i = 0; i < nkb; i++)
            {
                const std::complex<FPTYPE>* pvkb0i = vkbi + i * npwx;
                const std::complex<FPTYPE>* pvkb0j = vkbj + i * npwx;
                std::complex<FPTYPE>* pvkb = nullptr;
                std::complex<FPTYPE>* pdbecp_noevc = dbecp_noevc + i * npwx;

                // third term of dbecp_noevc
                //std::complex<FPTYPE>* pvkb = &vkb2(i,0);
                //std::complex<FPTYPE>* pdbecp_noevc = &dbecp_noevc(i, 0);
                FPTYPE qvec[3] = {0, 0, 0};
                for (int ig = 0; ig < npw; ig++)
                {
                    pvkb = vkb1 + i * npwx;
                    qvec[ipol] = gcar[(ik * npwx + ig) * 3 + ipol] + kvec_c[ik * 3 + ipol];
                    qvec[jpol] = gcar[(ik * npwx + ig) * 3 + jpol] + kvec_c[ik * 3 + jpol];
                    pvkb[ig] += 0.5 * qvec[ipol] * pvkb0j[ig] +
                                0.5 * qvec[jpol] * pvkb0i[ig];

                    pdbecp_noevc[ig] -= 2.0 * pvkb[ig];

                    if (ipol == jpol) {
                        pvkb = vkb + i * npwx;
                        pdbecp_noevc[ig] -= pvkb[ig];
                    }
                    pvkb = vkb2 + i * npwx;
                    for (int ii = 0; ii < 3; ii++) {
                        qvec[ii] = gcar[(ik * npwx + ig) * 3 + ii] + kvec_c[ik * 3 + ii];
                    }
                    FPTYPE qvec_norm2 = qvec[0] * qvec[0] + qvec[1] * qvec[1] + qvec[2] * qvec[2];
                    FPTYPE qm1 = qvec_norm2 > 1e-16 ? 1.0 / sqrt(qvec_norm2) : 0;
                    pdbecp_noevc[ig] -= 2.0 * pvkb[ig] * qvec[ipol] *
                                        qvec[jpol] * qm1 *	tpiba;
                } // end ig
            }//end nkb
        }
    };

    template <typename FPTYPE>
    struct cal_stress_nl_op<FPTYPE, psi::DEVICE_CPU> {
        void operator()(
                const psi::DEVICE_CPU *ctx,
                const bool &multi_proj,
                const int &ipol,
                const int &jpol,
                const int &nkb,
                const int &nbands_occ,
                const int &ntype,
                const int &spin,
                const int &wg_nc,
                const int &ik,
                const int &deeq_2,
                const int &deeq_3,
                const int &deeq_4,
                const int *atom_nh,
                const int *atom_na,
                const FPTYPE *d_wg,
                const FPTYPE *deeq,
                const std::complex<FPTYPE> *becp,
                const std::complex<FPTYPE> *dbecp,
                FPTYPE *stress)
        {
            for (int ib = 0; ib < nbands_occ; ib++)
            {
                FPTYPE fac = d_wg[ik * wg_nc + ib] * 1.0;
                int iat = 0, sum = 0;
                for (int it = 0; it < ntype; it++)
                {
                    const int Nprojs = atom_nh[it];
                    for (int ia=0; ia < atom_na[it]; ia++)
                    {
                        for (int ip1=0; ip1<Nprojs; ip1++)
                        {
                            for(int ip2=0; ip2<Nprojs; ip2++)
                            {
                                if(!multi_proj && ip1 != ip2) {
                                    continue;
                                }
                                FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip1) * deeq_4 + ip2];
                                const int inkb1 = sum + ip1;
                                const int inkb2 = sum + ip2;
                                //out<<"\n ps = "<<ps;


                                const FPTYPE dbb = ( conj( dbecp[ ib * nkb + inkb1] ) * becp[ ib * nkb + inkb2] ).real();
                                stress[ipol * 3 + jpol] -= ps * fac * dbb;
                            }
                        }//end ip
                        ++iat;
                        sum+=Nprojs;
                    }//ia
                } //end it
            }
        }
    };

    template struct cal_dbecp_noevc_nl_op<double, psi::DEVICE_CPU>;
    template struct cal_stress_nl_op<double, psi::DEVICE_CPU>;

}  // namespace src_pw

