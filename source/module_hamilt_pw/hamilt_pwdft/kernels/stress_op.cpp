#include "module_hamilt_pw/hamilt_pwdft/kernels/stress_op.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include "module_base/memory.h"
#include "module_base/math_polyint.h"

#include <iomanip>

namespace hamilt{

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
                for (int ii = 0; ii < 3; ii++) {
                    qvec[ii] = gcar[(ik * npwx + ig) * 3 + ii] + kvec_c[ik * 3 + ii];
                }
                auto pvkb1 = vkb1 + i * npwx;
                pvkb1[ig] += static_cast<FPTYPE>(0.5) * qvec[ipol] * pvkb0j[ig] +
                             static_cast<FPTYPE>(0.5) * qvec[jpol] * pvkb0i[ig];
                pdbecp_noevc[ig] -= static_cast<FPTYPE>(2.0) * pvkb1[ig];

                if (ipol == jpol) {
                    auto pvkb = vkb + i * npwx;
                    pdbecp_noevc[ig] -= pvkb[ig];
                }
                auto pvkb2 = vkb2 + i * npwx;
                
                FPTYPE qvec_norm2 = qvec[0] * qvec[0] + qvec[1] * qvec[1] + qvec[2] * qvec[2];
                FPTYPE qm1 = qvec_norm2 > 1e-16 ? 1.0 / sqrt(qvec_norm2) : 0;
                pdbecp_noevc[ig] -= static_cast<FPTYPE>(2.0) * pvkb2[ig] * qvec[ipol] *
                                        qvec[jpol] * qm1 *	tpiba;
            } // end ig
        }//end nkb
    }
};

template <typename FPTYPE>
struct cal_stress_nl_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(const psi::DEVICE_CPU* ctx,
                    const bool& nondiagonal,
                    const int& ipol,
                    const int& jpol,
                    const int& nkb,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& wg_nc,
                    const int& ik,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE* d_wg,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* stress)
    {
        FPTYPE local_stress = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:local_stress)
{
#endif
        int iat = 0, sum = 0;
        for (int it = 0; it < ntype; it++)
        {
            const int Nprojs = atom_nh[it];
#ifdef _OPENMP
#pragma omp for collapse(4)
#endif
            for (int ib = 0; ib < nbands_occ; ib++)
            {
                for (int ia=0; ia < atom_na[it]; ia++)
                {
                    for (int ip1=0; ip1<Nprojs; ip1++)
                    {
                        for(int ip2=0; ip2<Nprojs; ip2++)
                        {
                            if (!nondiagonal && ip1 != ip2)
                            {
                                continue;
                            }
                            FPTYPE fac = d_wg[ik * wg_nc + ib] * 1.0;
                            FPTYPE ekb_now = d_ekb[ik * wg_nc + ib];
                            FPTYPE ps = deeq[((spin * deeq_2 + iat + ia) * deeq_3 + ip1) * deeq_4 + ip2]
                                        - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip1 * deeq_4 + ip2];
                            const int inkb1 = sum + ia * Nprojs + ip1;
                            const int inkb2 = sum + ia * Nprojs + ip2;
                            //out<<"\n ps = "<<ps;


                            const FPTYPE dbb = ( conj( dbecp[ ib * nkb + inkb1] ) * becp[ ib * nkb + inkb2] ).real();
                            local_stress -= ps * fac * dbb;
                        }
                    }//end ip
                }//ia
            } 
            sum += atom_na[it] * Nprojs;
            iat += atom_na[it];
        } //end it
#ifdef _OPENMP
}
#endif
        stress[ipol * 3 + jpol] += local_stress;
    }
};

template <typename T, typename Device>
void cal_stress_mgga_op<T, Device>::operator()(
    const int& spin,
    const int& nrxx,
    const Real& w1,
    const T * gradwfc,
    Real * crosstaus)
{
    for (int ir = 0; ir < nrxx; ir++) {
        int ipol = 0;
        for (int ix = 0; ix < 3; ix++) {
            for (int iy = 0; iy < ix + 1; iy++) {
                crosstaus[spin * nrxx * 6 + ipol * nrxx + ir] 
                    += 2.0 * w1 
                    * (gradwfc[ix*nrxx + ir].real() * gradwfc[iy*nrxx + ir].real()
                    +  gradwfc[ix*nrxx + ir].imag() * gradwfc[iy*nrxx + ir].imag());
                ipol += 1;
            }
        }
    }
}


// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vkb_op<FPTYPE, psi::DEVICE_CPU>{
    void operator()(
        const psi::DEVICE_CPU* ctx,
        int nh,
        int npw,
        const FPTYPE** vqs_in,
        const FPTYPE** ylms_in,
        const std::complex<FPTYPE>* sk_in,
        const std::complex<FPTYPE>* pref_in,
        std::complex<FPTYPE>** vkbs_out
    ){
        // loop over all beta functions
        for(int ih=0;ih<nh;ih++)
        {
            std::complex<FPTYPE>* vkb_ptr = vkbs_out[ih];
            const FPTYPE* ylm_ptr = ylms_in[ih];
            const FPTYPE* vq_ptr = vqs_in[ih];
            // loop over all G-vectors
            for(int ig=0;ig<npw;ig++)
            {
                vkb_ptr[ig] = ylm_ptr[ig] * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
            }
        }
    }
};


// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vkb_deri_op<FPTYPE, psi::DEVICE_CPU>{
    void operator()(
        const psi::DEVICE_CPU* ctx,
        int nh,
        int npw,
        int ipol,
        int jpol,
        const FPTYPE** vqs_in, const FPTYPE** vqs_deri_in,
        const FPTYPE** ylms_in, const FPTYPE** ylms_deri_in1,const FPTYPE** ylms_deri_in2,
        const std::complex<FPTYPE>* sk_in,
        const std::complex<FPTYPE>* pref_in,
        const FPTYPE* gk_in,
        std::complex<FPTYPE>** vkbs_out
    ){
        int x1 = (GlobalC::ppcell.lmaxkb + 1) * (GlobalC::ppcell.lmaxkb + 1);
        int ih=0;
        // loop over all beta functions
        for(int ih=0;ih<nh;ih++)
        {
            std::complex<FPTYPE>* vkb_ptr = vkbs_out[ih];
            const FPTYPE* ylm_ptr = ylms_in[ih];
            const FPTYPE* vq_ptr = vqs_in[ih];
            // set vkb to zero
            for(int ig=0;ig<npw;ig++)
            {
                vkb_ptr[ig] = std::complex<FPTYPE>(0.0, 0.0);
            }
            // first term: ylm * vq * sk * pref
            // loop over all G-vectors
            if(ipol == jpol)
            {
                for(int ig=0;ig<npw;ig++)
                {
                    vkb_ptr[ig] -= ylm_ptr[ig] * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
                }
            }
            //second term: ylm_deri * vq_deri * sk * pref
            // loop over all G-vectors
            const FPTYPE* ylm_deri_ptr1 = &ylms_deri_in1[ih];
            const FPTYPE* ylm_deri_ptr2 = &ylms_deri_in2[ih];
            const FPTYPE* vq_deri_ptr = &vqs_deri_in[ih];
            const FPTYPE* gkn = &gk_in[4 * npw];
            for(int ig=0;ig<npw;ig++)
            {
                vkb_ptr[ig] -= (gk_in[ig*3+ipol] * ylm_deri_ptr2[ig] + gk_in[ig*3+jpol] * ylm_deri_ptr1[ig]) 
                                * vq_ptr[ig] * sk_in[ig] * pref_in[ih];
            }
            //third term: ylm * vq_deri * sk * pref
            // loop over all G-vectors
            for(int ig=0;ig<npw;ig++)
            {
                vkb_ptr[ig] -= 2.0 * ylm_ptr[ig] * vq_deri_ptr[ig] * sk_in[ig] * pref_in[ih]
                            * gk_in[ig*3+ipol] * gk_in[ig*3+jpol] * gkn[ig];
            }  
        }
    }
};

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vq_op<FPTYPE, psi::DEVICE_CPU>{
    void operator()(
        const psi::DEVICE_CPU* ctx,
        const FPTYPE* tab,
        int it, const FPTYPE* gk, int npw,
        const int tab_2,const int tab_3, const FPTYPE table_interval, 
        const int nbeta, FPTYPE* vq
    ){
        for (int nb = 0; nb < nbeta; nb++)
        {
            FPTYPE* vq_ptr = &vq[nb * npw];
            const FPTYPE* gnorm = &gk[3 * npw];
            for (int ig = 0; ig < npw; ig++)
            {
                vq_ptr[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(GlobalC::ppcell.tab,
                                                                        it,
                                                                        nb,
                                                                        GlobalV::NQX,
                                                                        GlobalV::DQ,
                                                                        gnorm[ig]);
            }
        }
    }
};


template <typename FPTYPE, typename Device>
FPTYPE Polynomial_Interpolation_nl
(
    const ModuleBase::realArray &table,
    const int &dim1,
    const int &dim2,
    const FPTYPE &table_interval,
    const FPTYPE &x                             // input value
)
{

	assert(table_interval>0.0);
	const FPTYPE position = x  / table_interval;
	const int iq = static_cast<int>(position);

	const FPTYPE x0 = position - static_cast<FPTYPE>(iq);
	const FPTYPE x1 = 1.0 - x0;
	const FPTYPE x2 = 2.0 - x0;
	const FPTYPE x3 = 3.0 - x0;
	const FPTYPE y=
			( table(dim1, dim2, iq)   * (-x2*x3-x1*x3-x1*x2) / 6.0 +
			table(dim1, dim2, iq+1) * (+x2*x3-x0*x3-x0*x2) / 2.0 -
			table(dim1, dim2, iq+2) * (+x1*x3-x0*x3-x0*x1) / 2.0 +
			table(dim1, dim2, iq+3) * (+x1*x2-x0*x2-x0*x1) / 6.0 )/table_interval ;


	return y;
}

// cpu version first, gpu version later
template <typename FPTYPE>
struct cal_vq_deri_op<FPTYPE, psi::DEVICE_CPU>{
    void operator()(
        const psi::DEVICE_CPU* ctx,
        const FPTYPE* tab,
        int it, const FPTYPE* gk, int npw,
        const int tab_2,const int tab_3, const FPTYPE &table_interval, 
        const int nbeta, FPTYPE* vq
    ){
        for (int nb = 0; nb < nbeta; nb++)
        {
            const FPTYPE* gnorm = &gk[3 * npw];
            FPTYPE* vq_ptr = &vq[nb * npw];
            for (int ig = 0; ig < npw; ig++)
            {
                vq_ptr[ig] = Polynomial_Interpolation_nl<FPTYPE, psi::DEVICE_CPU>(
                            GlobalC::ppcell.tab, 
                            it, 
                            nb, 
                            GlobalV::DQ, 
                            gnorm[ig] );
            }
        }
        return ;
    }
};




template struct cal_stress_mgga_op<std::complex<float>,  psi::DEVICE_CPU>;
template struct cal_stress_mgga_op<std::complex<double>, psi::DEVICE_CPU>;

template struct cal_dbecp_noevc_nl_op<float, psi::DEVICE_CPU>;
template struct cal_stress_nl_op<float, psi::DEVICE_CPU>;

template struct cal_dbecp_noevc_nl_op<double, psi::DEVICE_CPU>;
template struct cal_stress_nl_op<double, psi::DEVICE_CPU>;

template struct cal_vkb_op<float, psi::DEVICE_CPU>;
template struct cal_vkb_op<double, psi::DEVICE_CPU>;

template struct cal_vq_op<float, psi::DEVICE_CPU>;
template struct cal_vq_op<double, psi::DEVICE_CPU>;


template struct cal_vq_deri_op<float, psi::DEVICE_CPU>;
template struct cal_vq_deri_op<double, psi::DEVICE_CPU>;
}  // namespace hamilt

