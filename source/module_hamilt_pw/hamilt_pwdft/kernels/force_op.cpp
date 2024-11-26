#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hamilt
{

template <typename FPTYPE>
struct cal_vkb1_nl_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nkb,
                    const int& npwx,
                    const int& vkb_nc,
                    const int& nbasis,
                    const int& ipol,
                    const std::complex<FPTYPE>& NEG_IMAG_UNIT,
                    const std::complex<FPTYPE>* vkb,
                    const FPTYPE* gcar,
                    std::complex<FPTYPE>* vkb1)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int i = 0; i < nkb; i++)
        {
            for (int ig = 0; ig < nbasis; ig++)
            {
                std::complex<FPTYPE>* pvkb1 = vkb1 + i * npwx;
                const std::complex<FPTYPE>* pvkb = vkb + i * vkb_nc;
                pvkb1[ig] = pvkb[ig] * NEG_IMAG_UNIT * gcar[ig * 3 + ipol];
            }
        }
    }
};

template <typename FPTYPE>
struct cal_force_nl_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const bool& nondiagonal,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
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
            for (int it = 0; it < ntype; it++)
            {
                const int nproj = atom_nh[it];
#ifdef _OPENMP
#pragma omp for collapse(2)
#endif
                for (int ia = 0; ia < atom_na[it]; ia++)
                {
                    for (int ib = 0; ib < nbands_occ; ib++)
                    {
                        FPTYPE local_force[3] = {0, 0, 0};
                        FPTYPE fac;
                        if(occ)
                        {
                            fac = d_wg[ib] * 2.0 * tpiba;
                        }
                        else
                        {
                            fac = d_wg[0] * 2.0 * tpiba;
                        }
                        FPTYPE ekb_now = 0.0;
                        if (d_ekb != nullptr)
                        {
                            ekb_now = d_ekb[ib];
                        }
                        int iat = iat0 + ia;
                        int sum = sum0 + ia * nproj;
                        for (int ip = 0; ip < nproj; ip++)
                        {
                            FPTYPE ps_qq = 0;
                            if(ekb_now != 0)
                            {
                                ps_qq = - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip];
                            }
                            // Effective values of the D-eS coefficients
                            FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip] + ps_qq;
                            const int inkb = sum + ip;
                            // out<<"\n ps = "<<ps;

                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                const FPTYPE dbb
                                    = (conj(dbecp[ipol * nbands * nkb + ib * nkb + inkb]) * becp[ib * nkb + inkb])
                                          .real();
                                local_force[ipol] -= ps * fac * dbb;
                                // cf[iat*3+ipol] += ps * fac * dbb;
                            }
                            if (nondiagonal)
                            {
                                for (int ip2 = 0; ip2 < nproj; ip2++)
                                {
                                    if (ip != ip2)
                                    {
                                        const int jnkb = sum + ip2;
                                        FPTYPE ps_qq = 0;
                                        if (ekb_now != 0)
                                        {
                                            ps_qq = -ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip2];
                                        }
                                        FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;
                                        for (int ipol = 0; ipol < 3; ipol++)
                                        {
                                            const FPTYPE dbb = (conj(dbecp[ipol * nbands * nkb + ib * nkb + inkb])
                                                                * becp[ib * nkb + jnkb])
                                                                   .real();
                                            local_force[ipol] -= ps * fac * dbb;
                                        }
                                    }
                                }
                            }
                        }
#ifdef _OPENMP
                        if (omp_get_num_threads() > 1)
                        {
                            for (int ipol = 0; ipol < 3; ipol++)
                            {
#pragma omp atomic
                                force[iat * forcenl_nc + ipol] += local_force[ipol];
                            }
                        }
                        else
#endif
                        {
                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                force[iat * forcenl_nc + ipol] += local_force[ipol];
                            }
                        }
                    }
                } // end ia
                iat0 += atom_na[it];
                sum0 += atom_na[it] * nproj;
            } // end it
#ifdef _OPENMP
        }
#endif
    };

    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const std::complex<FPTYPE>* deeq_nc,
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
            for (int it = 0; it < ntype; it++)
            {
                const int nprojs = atom_nh[it];
#ifdef _OPENMP
#pragma omp for collapse(2)
#endif
                for (int ia = 0; ia < atom_na[it]; ia++)
                {
                    for (int ib = 0; ib < nbands_occ; ib++)
                    {
                        const int ib2 = ib*2;
                        FPTYPE local_force[3] = {0, 0, 0};
                        FPTYPE fac;
                        if(occ)
                        {
                            fac = d_wg[ib] * 2.0 * tpiba;
                        }
                        else
                        {
                            fac = d_wg[0] * 2.0 * tpiba;
                        }
                        FPTYPE ekb_now = 0.0;
                        if (d_ekb != nullptr)
                        {
                            ekb_now = d_ekb[ib];
                        }
                        int iat = iat0 + ia;
                        int sum = sum0 + ia * nprojs;
                        for (int ip = 0; ip < nprojs; ip++)
                        {
                            const int inkb = sum + ip;
                            // out<<"\n ps = "<<ps;
                            for (int ip2 = 0; ip2 < nprojs; ip2++)
                            {
                                // Effective values of the D-eS coefficients
                                std::complex<FPTYPE> ps_qq = 0;
                                if(ekb_now)
                                {
                                    ps_qq = - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip2];
                                }
                                const int jnkb = sum + ip2;
                                std::complex<FPTYPE> ps0 = deeq_nc[((0 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;
                                std::complex<FPTYPE> ps1 = deeq_nc[((1 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2];
                                std::complex<FPTYPE> ps2 = deeq_nc[((2 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2];
                                std::complex<FPTYPE> ps3 = deeq_nc[((3 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;
                                

                                for (int ipol = 0; ipol < 3; ipol++)
                                {
                                    const int index0 = ipol * nbands * 2 * nkb + ib2 * nkb + inkb;
                                    const int index1 = ib2 * nkb + jnkb;
                                    const std::complex<FPTYPE> dbb0 = conj(dbecp[index0]) * becp[index1];
                                    const std::complex<FPTYPE> dbb1 = conj(dbecp[index0]) * becp[index1 + nkb];
                                    const std::complex<FPTYPE> dbb2 = conj(dbecp[index0 + nkb]) * becp[index1];
                                    const std::complex<FPTYPE> dbb3 = conj(dbecp[index0 + nkb]) * becp[index1 + nkb];

                                    local_force[ipol] -= fac * (ps0 * dbb0 + ps1 * dbb1 + ps2 * dbb2 + ps3 * dbb3).real();
                                }
                            }
                        }
#ifdef _OPENMP
                        if (omp_get_num_threads() > 1)
                        {
                            for (int ipol = 0; ipol < 3; ++ipol)
                            {
#pragma omp atomic
                                force[iat * forcenl_nc + ipol] += local_force[ipol];
                            }
                        }
                        else
#endif
                        {
                            for (int ipol = 0; ipol < 3; ++ipol)
                            {
                                force[iat * forcenl_nc + ipol] += local_force[ipol];
                            }
                        }
                    }
                } // end ia
                iat0 += atom_na[it];
                sum0 += atom_na[it] * nprojs;
            } // end it
#ifdef _OPENMP
        }
#endif
    }
};

template struct cal_vkb1_nl_op<float, base_device::DEVICE_CPU>;
template struct cal_force_nl_op<float, base_device::DEVICE_CPU>;

template struct cal_vkb1_nl_op<double, base_device::DEVICE_CPU>;
template struct cal_force_nl_op<double, base_device::DEVICE_CPU>;

} // namespace hamilt
