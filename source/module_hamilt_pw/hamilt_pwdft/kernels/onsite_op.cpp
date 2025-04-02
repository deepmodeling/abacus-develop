#include "module_hamilt_pw/hamilt_pwdft/kernels/onsite_op.h"

namespace hamilt
{

template <typename FPTYPE>
struct onsite_ps_op<FPTYPE, base_device::DEVICE_CPU>
{
    // kernel for DeltaSpin calculation
    void operator()(const base_device::DEVICE_CPU* /*dev*/,
                    const int& npm,
                    const int npol,
                    const int* ip_iat,
                    const int& tnp,
                    const std::complex<FPTYPE>* lambda_array,
                    std::complex<FPTYPE>* ps,
                    const std::complex<FPTYPE>* becp)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int ib = 0; ib < npm / npol; ib++)
        {
            for (int ip = 0; ip < tnp; ip++)
            {
                int ib2 = ib * npol;
                int iat = ip_iat[ip];
                const int psind = ip * npm + ib2;
                const int becpind = ib2 * tnp + ip;
                ps[psind] += lambda_array[iat * 4] * becp[becpind] 
                            + lambda_array[iat * 4 + 2] * becp[becpind + tnp];
                ps[psind + 1] += lambda_array[iat * 4 + 1] * becp[becpind] 
                            + lambda_array[iat * 4 + 3] * becp[becpind + tnp];
            } // end ip
        } // end ib
    };

    // kernel for DFT+U calculation
    void operator()(const base_device::DEVICE_CPU* dev,
      const int& npm,
      const int npol,
      const int* orb_l_iat,
      const int* ip_iat,
      const int* ip_m,
      const int* vu_begin_iat,
      const int& tnp,
      const std::complex<FPTYPE>* vu,
      std::complex<FPTYPE>* ps,
      const std::complex<FPTYPE>* becp)
  {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int ib = 0; ib < npm / npol; ib++)
        {
            for (int ip = 0; ip < tnp; ip++)
            {
                int m1 = ip_m[ip];
                if(m1 < 0) continue;
                int ib2 = ib * npol;
                int iat = ip_iat[ip];
                const std::complex<FPTYPE>* vu_iat = vu + vu_begin_iat[iat];
                int orb_l = orb_l_iat[iat];
                int tlp1 = 2 * orb_l + 1;
                int tlp1_2 = tlp1 * tlp1;
                int ip2_begin = ip - m1;
                int ip2_end = ip - m1 + tlp1;
                const int psind = ip * npm + ib2;
                for(int ip2 = ip2_begin;ip2<ip2_end;ip2++)
                {
                    const int becpind = ib2 * tnp + ip2;
                    int m2 = ip_m[ip2];
                    const int index_mm = m1 * tlp1 + m2;
                    ps[psind] += vu_iat[index_mm] * becp[becpind]
                                + vu_iat[index_mm + tlp1_2 * 2] * becp[becpind + tnp];
                    ps[psind + 1] += vu_iat[index_mm + tlp1_2 * 1] * becp[becpind]
                                + vu_iat[index_mm + tlp1_2 * 3] * becp[becpind + tnp];
                }
            } // end ip
        } // end ib
  }
};

template struct onsite_ps_op<float, base_device::DEVICE_CPU>;
template struct onsite_ps_op<double, base_device::DEVICE_CPU>;

} // namespace hamilt