#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"

#include "module_base/math_polyint.h"
#include "module_base/memory.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/vnl_op.h"
#include "vnl_tools.hpp"
#include <vector>
#include <iomanip>
#include "module_base/libm/libm.h"
namespace hamilt
{
template <typename FPTYPE>
struct cal_dbecp_noevc_nl_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& npw,
                    const int& npwx,
                    const int& ik,
                    const FPTYPE& tpiba,
                    const FPTYPE* gcar,
                    const FPTYPE* kvec_c,
                    std::complex<FPTYPE>* vkbi,
                    std::complex<FPTYPE>* vkbj,
                    std::complex<FPTYPE>* vkb,
                    std::complex<FPTYPE>* vkb1,
                    std::complex<FPTYPE>* vkb2,
                    std::complex<FPTYPE>* dbecp_noevc)
    {
        // npwx >= npw
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int i = 0; i < nkb; i++)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                auto pvkb0i = vkbi + i * npwx;
                auto pvkb0j = vkbj + i * npwx;
                auto pdbecp_noevc = dbecp_noevc + i * npwx;
                FPTYPE qvec[3];
                for (int ii = 0; ii < 3; ii++)
                {
                    qvec[ii] = gcar[(ik * npwx + ig) * 3 + ii] + kvec_c[ik * 3 + ii];
                }
                auto pvkb1 = vkb1 + i * npwx;
                pvkb1[ig] += static_cast<FPTYPE>(0.5) * qvec[ipol] * pvkb0j[ig]
                             + static_cast<FPTYPE>(0.5) * qvec[jpol] * pvkb0i[ig];
                pdbecp_noevc[ig] -= static_cast<FPTYPE>(2.0) * pvkb1[ig];

                if (ipol == jpol)
                {
                    auto pvkb = vkb + i * npwx;
                    pdbecp_noevc[ig] -= pvkb[ig];
                }
                auto pvkb2 = vkb2 + i * npwx;

                FPTYPE qvec_norm2 = qvec[0] * qvec[0] + qvec[1] * qvec[1] + qvec[2] * qvec[2];
                FPTYPE qm1 = qvec_norm2 > 1e-16 ? 1.0 / sqrt(qvec_norm2) : 0;
                pdbecp_noevc[ig] -= static_cast<FPTYPE>(2.0) * pvkb2[ig] * qvec[ipol] * qvec[jpol] * qm1 * tpiba;
            } // end ig
        }     // end nkb
    }
};

template <typename FPTYPE>
struct cal_stress_nl_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const bool& nondiagonal,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* stress)
    {
        FPTYPE local_stress = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+ : local_stress)
        {
#endif
            int iat = 0;
            int sum = 0;
            for (int it = 0; it < ntype; it++)
            {
                const int nproj = atom_nh[it];
#ifdef _OPENMP
#pragma omp for
#endif
                for (int ib = 0; ib < nbands_occ; ib++)
                {
                    FPTYPE ekb_now = 0.0;
                    if(d_ekb != nullptr)
                    {
                        ekb_now = d_ekb[ib];
                    }
                    for (int ia = 0; ia < atom_na[it]; ia++)
                    {
                        for (int ip1 = 0; ip1 < nproj; ip1++)
                        {
                            for (int ip2 = 0; ip2 < nproj; ip2++)
                            {
                                if (!nondiagonal && ip1 != ip2)
                                {
                                    continue;
                                }

                                FPTYPE fac;
                                if (occ)
                                {
                                    fac = d_wg[ib];
                                }
                                else
                                {
                                    fac = d_wg[0];
                                }
                                FPTYPE ps_qq = 0;
                                if(ekb_now != 0)
                                {
                                    ps_qq = - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip1 * deeq_4 + ip2];
                                }
                                FPTYPE ps = deeq[((spin * deeq_2 + iat + ia) * deeq_3 + ip1) * deeq_4 + ip2] + ps_qq;
                                const int inkb1 = sum + ia * nproj + ip1;
                                const int inkb2 = sum + ia * nproj + ip2;
                                // out<<"\n ps = "<<ps;

                                const FPTYPE dbb = (conj(dbecp[ib * nkb + inkb1]) * becp[ib * nkb + inkb2]).real();
                                local_stress -= ps * fac * dbb;
                            }
                        } // end ip
                    }     // ia
                }
                sum += atom_na[it] * nproj;
                iat += atom_na[it];
            } // end it
#ifdef _OPENMP
        }
#endif
        stress[ipol * 3 + jpol] += local_stress;
    };

    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const std::complex<FPTYPE>* deeq_nc,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* stress)
    {
        FPTYPE local_stress = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+ : local_stress)
        {
#endif
            int iat = 0, sum = 0;
            for (int it = 0; it < ntype; it++)
            {
                const int nproj = atom_nh[it];
#ifdef _OPENMP
#pragma omp for
#endif
                for (int ib = 0; ib < nbands_occ; ib++)
                {
                    FPTYPE ekb_now = 0.0;
                    if(d_ekb != nullptr)
                    {
                        ekb_now = d_ekb[ib];
                    }
                    for (int ia = 0; ia < atom_na[it]; ia++)
                    {
                        for (int ip1 = 0; ip1 < nproj; ip1++)
                        {
                            for (int ip2 = 0; ip2 < nproj; ip2++)
                            {
                                const int ib2 = ib*2;
                                FPTYPE fac;
                                if (occ)
                                {
                                    fac = d_wg[ib];
                                }
                                else
                                {
                                    fac = d_wg[0];
                                }
                                std::complex<FPTYPE> ps_qq = 0;
                                if(ekb_now != 0)
                                {
                                    ps_qq = - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip1 * deeq_4 + ip2];
                                }
                                std::complex<FPTYPE> ps0 = deeq_nc[((iat + ia) * deeq_3 + ip1) * deeq_4 + ip2] + ps_qq;
                                std::complex<FPTYPE> ps1 = deeq_nc[((1 * deeq_2 + iat + ia) * deeq_3 + ip1) * deeq_4 + ip2];
                                std::complex<FPTYPE> ps2 = deeq_nc[((2 * deeq_2 + iat + ia) * deeq_3 + ip1) * deeq_4 + ip2];
                                std::complex<FPTYPE> ps3 = deeq_nc[((3 * deeq_2 + iat + ia) * deeq_3 + ip1) * deeq_4 + ip2] + ps_qq;
                                const int inkb1 = sum + ia * nproj + ip1;
                                const int inkb2 = sum + ia * nproj + ip2;
                                // out<<"\n ps = "<<ps;

                                const std::complex<FPTYPE> dbb0 = conj(dbecp[ib2 * nkb + inkb1]) * becp[ib2 * nkb + inkb2];
                                const std::complex<FPTYPE> dbb1 = conj(dbecp[ib2 * nkb + inkb1]) * becp[(ib2+1) * nkb + inkb2];
                                const std::complex<FPTYPE> dbb2 = conj(dbecp[(ib2+1) * nkb + inkb1]) * becp[ib2 * nkb + inkb2];
                                const std::complex<FPTYPE> dbb3 = conj(dbecp[(ib2+1) * nkb + inkb1]) * becp[(ib2+1) * nkb + inkb2];
                                local_stress -= fac * (ps0 * dbb0 + ps1 * dbb1 + ps2 * dbb2 + ps3 * dbb3).real();
                            }
                        } // end ip
                    }     // ia
                }
                sum += atom_na[it] * nproj;
                iat += atom_na[it];
            } // end it
#ifdef _OPENMP
        }
#endif
        stress[ipol * 3 + jpol] += local_stress;
    }
};

template <typename T, typename Device>
void cal_stress_mgga_op<T, Device>::operator()(const int& spin,
                                               const int& nrxx,
                                               const Real& w1,
                                               const T* gradwfc,
                                               Real* crosstaus)
{
    for (int ir = 0; ir < nrxx; ir++)
    {
        int ipol = 0;
        for (int ix = 0; ix < 3; ix++)
        {
            for (int iy = 0; iy < ix + 1; iy++)
            {
                crosstaus[spin * nrxx * 6 + ipol * nrxx + ir]
                    += 2.0 * w1
                       * (gradwfc[ix * nrxx + ir].real() * gradwfc[iy * nrxx + ir].real()
                          + gradwfc[ix * nrxx + ir].imag() * gradwfc[iy * nrxx + ir].imag());
                ipol += 1;
            }
        }
    }
}

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vkb_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int nh,
                    const int npw,
                    const int* indexes,
                    const FPTYPE* vqs_in,
                    const FPTYPE* ylms_in,
                    const std::complex<FPTYPE>* sk_in,
                    const std::complex<FPTYPE>* pref_in,
                    std::complex<FPTYPE>* vkbs_out)
    {
        // loop over all beta functions
        for (int ih = 0; ih < nh; ih++)
        {
            std::complex<FPTYPE>* vkb_ptr = vkbs_out + ih * npw;
            const FPTYPE* ylm_ptr = ylms_in + indexes[ih * 4] * npw;
            const FPTYPE* vq_ptr = vqs_in + indexes[ih * 4 + 1] * npw;
            // loop over all G-vectors
            for (int ig = 0; ig < npw; ig++)
            {
                vkb_ptr[ig] = ylm_ptr[ig] * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
            }
        }
    }
};

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vkb_deri_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int nh,
                    const int npw,
                    const int ipol,
                    const int jpol,
                    const int* indexes,
                    const FPTYPE* vqs_in,
                    const FPTYPE* vqs_deri_in,
                    const FPTYPE* ylms_in,
                    const FPTYPE* ylms_deri_in,
                    const std::complex<FPTYPE>* sk_in,
                    const std::complex<FPTYPE>* pref_in,
                    const FPTYPE* gk_in,
                    std::complex<FPTYPE>* vkbs_out)
    {
        int ih = 0;
        // loop over all beta functions
        for (int ih = 0; ih < nh; ih++)
        {
            // move ptrs
            std::complex<FPTYPE>* vkb_ptr = vkbs_out + ih * npw;
            const FPTYPE* ylm_ptr = ylms_in + indexes[ih * 4] * npw;
            const FPTYPE* vq_ptr = vqs_in + indexes[ih * 4 + 1] * npw;
            // set vkb to zero
            for (int ig = 0; ig < npw; ig++)
            {
                vkb_ptr[ig] = std::complex<FPTYPE>(0.0, 0.0);
            }
            // first term: ylm * vq * sk * pref
            // loop over all G-vectors
            if (ipol == jpol)
            {
                for (int ig = 0; ig < npw; ig++)
                {
                    vkb_ptr[ig] -= ylm_ptr[ig] * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
                }
            }
            // second term: ylm_deri * vq_deri * sk * pref
            //  loop over all G-vectors
            const FPTYPE* ylm_deri_ptr1 = ylms_deri_in + indexes[ih * 4 + 2] * npw;
            const FPTYPE* ylm_deri_ptr2 = ylms_deri_in + indexes[ih * 4 + 3] * npw;
            const FPTYPE* vq_deri_ptr = vqs_deri_in + indexes[ih * 4 + 1] * npw;
            const FPTYPE* gkn = &gk_in[4 * npw];
            for (int ig = 0; ig < npw; ig++)
            {
                vkb_ptr[ig] -= (gk_in[ig * 3 + ipol] * ylm_deri_ptr2[ig] + gk_in[ig * 3 + jpol] * ylm_deri_ptr1[ig])
                               * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
            }
            // third term: ylm * vq_deri * sk * pref
            //  loop over all G-vectors
            for (int ig = 0; ig < npw; ig++)
            {
                FPTYPE two = 2.0;
                vkb_ptr[ig] -= two * ylm_ptr[ig] * vq_deri_ptr[ig] * sk_in[ig] * pref_in[ih] * gk_in[ig * 3 + ipol]
                               * gk_in[ig * 3 + jpol] * gkn[ig];
            }
        }
    }
};

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vq_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const FPTYPE* tab,
                    int it,
                    const FPTYPE* gk,
                    int npw,
                    const int tab_2,
                    const int tab_3,
                    const FPTYPE table_interval,
                    const int nbeta,
                    FPTYPE* vq)
    {
        for (int nb = 0; nb < nbeta; nb++)
        {
            FPTYPE* vq_ptr = &vq[nb * npw];
            const FPTYPE* gnorm = &gk[3 * npw];
            for (int ig = 0; ig < npw; ig++)
            {
                vq_ptr[ig] = _polynomial_interpolation<FPTYPE>(tab, it, nb, tab_2, tab_3, table_interval, gnorm[ig]);
            }
        }
    }
};

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vq_deri_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const FPTYPE* tab,
                    int it,
                    const FPTYPE* gk,
                    int npw,
                    const int tab_2,
                    const int tab_3,
                    const FPTYPE table_interval,
                    const int nbeta,
                    FPTYPE* vq)
    {
        for (int nb = 0; nb < nbeta; nb++)
        {
            const FPTYPE* gnorm = &gk[3 * npw];
            FPTYPE* vq_ptr = &vq[nb * npw];
            for (int ig = 0; ig < npw; ig++)
            {
                vq_ptr[ig] = _polynomial_interpolation_nl<FPTYPE>(tab, it, nb, tab_2, tab_3, table_interval, gnorm[ig]);
            }
        }
        return;
    }
};


template <typename FPTYPE>
void Simpson_Integral
(
    const int mesh,
    FPTYPE * func,
    const FPTYPE * rab,
    FPTYPE &asum
)
{
    assert(mesh&1);

    asum = 0.00;
	const size_t end = mesh-2;
    for( size_t i=1; i!=end; i+=2 )
    {
		const double f1 = func[i]*rab[i];
		asum += f1 + f1 + func[i+1]*rab[i+1];
    }
	const double f1 = func[mesh-2]*rab[mesh-2];
	asum += f1+f1;
	asum += asum;
	asum += func[0]*rab[0] + func[mesh-1]*rab[mesh-1];
	asum /= 3.0;
    return;
}// end subroutine simpson

template <typename FPTYPE>
struct cal_stress_drhoc_aux_op<FPTYPE, base_device::DEVICE_CPU> {
    void operator()(const FPTYPE* r,
                    const FPTYPE* rhoc,
                    const FPTYPE* gx_arr,
                    const FPTYPE* rab,
                    FPTYPE* drhocg,
                    const int mesh,
                    const int igl0,
                    const int ngg,
                    const double omega,
                    int type) {
    const double FOUR_PI = 4.0 * 3.14159265358979323846;
    // printf("%d,%d,%lf\n",ngg,mesh,omega);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif

#ifdef _OPENMP
#pragma omp for
#endif
        for(int igl = 0;igl< ngg;igl++)
        {
            FPTYPE rhocg1 = 0;
            //FPTYPE *aux = new FPTYPE[mesh];
            std::vector<FPTYPE> aux(mesh);
            for( int ir = 0;ir< mesh; ir++)
            {
                if(type ==0 ){
                    aux [ir] = r [ir] * rhoc [ir] * (r [ir] * cos (gx_arr[igl] * r [ir] ) / gx_arr[igl] - sin (gx_arr[igl] * r [ir] ) / pow(gx_arr[igl],2));
                } else if(type == 1) {
                    aux [ir] = ir!=0 ? std::sin(gx_arr[igl] * r[ir]) / (gx_arr[igl] * r[ir]) : 1.0;
                    aux [ir] = r[ir] * r[ir] * rhoc [ir] * aux [ir];
                } else if(type == 2) {
                    aux [ir] = r[ir] < 1.0e-8 ? rhoc [ir] : rhoc [ir] * sin(gx_arr[igl] * r[ir]) / (gx_arr[igl] * r[ir]);
                } else if(type == 3) {
                    FPTYPE sinp, cosp;
                    sinp = std::sin(gx_arr[igl] * r[ir]);
                    cosp = std::cos(gx_arr[igl] * r[ir]);
                    aux[ir] = rhoc [ir] *  (r [ir] * cosp / gx_arr[igl] - sinp / pow(gx_arr[igl],2));
                }
            }//ir
            Simpson_Integral<FPTYPE>(mesh, aux.data(), rab, rhocg1);
            if(type ==0 ) { drhocg [igl] = FOUR_PI / omega * rhocg1;
            } else if(type == 1) { drhocg [igl] = FOUR_PI * rhocg1 / omega;
            } else if(type == 2) { drhocg [igl] = rhocg1;
            } else if(type == 3) {
                rhocg1 *= FOUR_PI / omega / 2.0 / gx_arr[igl];
                FPTYPE g2a = (gx_arr[igl]*gx_arr[igl]) / 4.0;
                rhocg1 += FOUR_PI / omega * gx_arr[ngg] * ModuleBase::libm::exp ( - g2a) * (g2a + 1) / pow(gx_arr[igl]*gx_arr[igl] , 2);
                drhocg [igl] = rhocg1;
            }
        }
#ifdef _OPENMP
        }
#endif
    }
};


template <typename FPTYPE>
struct cal_force_npw_op<FPTYPE, base_device::DEVICE_CPU> {
    void operator()(const std::complex<FPTYPE> *psiv,
                    const FPTYPE* gv_x, const FPTYPE* gv_y, const FPTYPE* gv_z,
                    const FPTYPE* rhocgigg_vec,
                    FPTYPE* force,
                    const FPTYPE pos_x, const FPTYPE pos_y, const FPTYPE pos_z,
                    const int npw,
                    const FPTYPE omega, const FPTYPE tpiba) {
        const double TWO_PI = 2.0 * 3.14159265358979323846;
        // printf("%d,%d,%lf\n",ngg,mesh,omega);

        // printf("iminininininin\n\n\n\n");
#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (int ig = 0; ig < npw; ig++)
        {
            const std::complex<FPTYPE> psiv_conj = conj(psiv[ig]);

            const FPTYPE arg = TWO_PI * (gv_x[ig] * pos_x + gv_y[ig] * pos_y + gv_z[ig] * pos_z);
            FPTYPE sinp, cosp;
            ModuleBase::libm::sincos(arg, &sinp, &cosp);
            const std::complex<FPTYPE> expiarg = std::complex<FPTYPE>(sinp, cosp);

            const std::complex<FPTYPE> tmp_var = psiv_conj * expiarg * tpiba * omega * rhocgigg_vec[ig];

            const std::complex<FPTYPE> ipol0 = tmp_var * gv_x[ig];
            force[0] += ipol0.real();

            const std::complex<FPTYPE> ipol1 = tmp_var * gv_y[ig];
            force[1] += ipol1.real();

            const std::complex<FPTYPE> ipol2 = tmp_var * gv_z[ig];
            force[2] += ipol2.real();
        }
    }
};

template <typename FPTYPE>
struct cal_multi_dot_op<FPTYPE, base_device::DEVICE_CPU> {
    FPTYPE operator()(const int& npw,
                      const FPTYPE& fac,
                      const FPTYPE* gk1,
                      const FPTYPE* gk2,
                      const FPTYPE* d_kfac,
                      const std::complex<FPTYPE>* psi)
    {
        FPTYPE sum = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
        for (int i = 0; i < npw; i++)
        {
            sum += fac * gk1[i] * gk2[i] * d_kfac[i] * std::norm(psi[i]);
        }
        return sum;
    }
};

// // cpu version first, gpu version later
// template <typename FPTYPE>
// struct prepare_vkb_deri_ptr_op<FPTYPE, base_device::DEVICE_CPU>{
//     void operator()(
//         const base_device::DEVICE_CPU* ctx,
//         int nbeta, double* nhtol, int nhtol_nc, int npw, int it,
//         int ipol, int jpol,
//         std::complex<FPTYPE>*vkb_out, std::complex<FPTYPE>** vkb_ptrs,
//         FPTYPE* ylm_in, FPTYPE** ylm_ptrs,
//         FPTYPE* ylm_deri_in, FPTYPE** ylm_deri_ptr1s, FPTYPE** ylm_deri_ptr2s,
//         FPTYPE* vq_in, FPTYPE** vq_ptrs,
//         FPTYPE* vq_deri_in, FPTYPE** vq_deri_ptrs
//     ){
//         // int ih=0;
//         // int x1 = (nlpp->lmaxkb + 1) * (nlpp->lmaxkb + 1);
//         // for(int nb=0;nb<nbeta;nb++)
//         // {
//         //     int l = nhtol[it*nhtol_nc+ih];
//         //     for(int m=0;m<2*l+1;m++)
//         //     {
//         //         int lm = l*l + m;
//         //         vkb_ptrs[ih] = &vkb_out[ih * npw];
//         //         ylm_ptrs[ih] = &ylm_in[lm * npw];
//         //         vq_ptrs[ih] = &vq_in[nb * npw];

//         //         ylm_deri_ptr1s[ih] = &ylm_deri_in[(ipol * x1 + lm) * npw];
//         //         ylm_deri_ptr2s[ih] = &ylm_deri_in[(jpol * x1 + lm) * npw];
//         //         vq_deri_ptrs[ih] = &vq_deri_in[nb * npw];

//         //         ih++;

//         //     }
//         // }
//     }
// };

template <>
void pointer_array_malloc<base_device::DEVICE_CPU>::operator()(void** ptr, const int n)
{
    return;
}

template <>
void synchronize_ptrs<base_device::DEVICE_CPU>::operator()(void** ptr_out, const void** ptr_in, const int size)
{
    return;
}

template struct pointer_array_malloc<base_device::DEVICE_CPU>;
template struct synchronize_ptrs<base_device::DEVICE_CPU>;

template struct cal_stress_mgga_op<std::complex<float>, base_device::DEVICE_CPU>;
template struct cal_stress_mgga_op<std::complex<double>, base_device::DEVICE_CPU>;

template struct cal_dbecp_noevc_nl_op<float, base_device::DEVICE_CPU>;
template struct cal_stress_nl_op<float, base_device::DEVICE_CPU>;

template struct cal_dbecp_noevc_nl_op<double, base_device::DEVICE_CPU>;
template struct cal_stress_nl_op<double, base_device::DEVICE_CPU>;

template struct cal_vkb_op<float, base_device::DEVICE_CPU>;
template struct cal_vkb_op<double, base_device::DEVICE_CPU>;

template struct cal_vkb_deri_op<float, base_device::DEVICE_CPU>;
template struct cal_vkb_deri_op<double, base_device::DEVICE_CPU>;

template struct cal_vq_op<float, base_device::DEVICE_CPU>;
template struct cal_vq_op<double, base_device::DEVICE_CPU>;

template struct cal_vq_deri_op<float, base_device::DEVICE_CPU>;
template struct cal_vq_deri_op<double, base_device::DEVICE_CPU>;

template struct cal_stress_drhoc_aux_op<float, base_device::DEVICE_CPU>;
template struct cal_stress_drhoc_aux_op<double, base_device::DEVICE_CPU>;

template struct cal_force_npw_op<float, base_device::DEVICE_CPU>;
template struct cal_force_npw_op<double, base_device::DEVICE_CPU>;

template struct cal_multi_dot_op<float, base_device::DEVICE_CPU>;
template struct cal_multi_dot_op<double, base_device::DEVICE_CPU>;


// template struct prepare_vkb_deri_ptr_op<float, base_device::DEVICE_CPU>;
// template struct prepare_vkb_deri_ptr_op<double, base_device::DEVICE_CPU>;
} // namespace hamilt
