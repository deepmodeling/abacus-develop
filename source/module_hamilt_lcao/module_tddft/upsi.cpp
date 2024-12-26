#include "upsi.h"

#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
#ifdef __MPI
void upsi(const Parallel_Orbitals* pv,
          const int nband,
          const int nlocal,
          const std::complex<double>* U_operator,
          const std::complex<double>* psi_k_laststep,
          std::complex<double>* psi_k,
          const int print_matrix)
{

    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             U_operator,
                             1,
                             1,
                             pv->desc,
                             psi_k_laststep,
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             psi_k,
                             1,
                             1,
                             pv->desc_wfc);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = psi_k[i * pv->ncol + j].real();
                bb = psi_k[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                    aa = 0.0;
                if (std::abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k_laststep:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = psi_k_laststep[i * pv->ncol + j].real();
                bb = psi_k_laststep[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                    aa = 0.0;
                if (std::abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
}

void upsi_tensor(const Parallel_Orbitals* pv,
                 const int nband,
                 const int nlocal,
                 const container::Tensor& U_operator,
                 const container::Tensor& psi_k_laststep,
                 container::Tensor& psi_k,
                 const int print_matrix)
{

    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             U_operator.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc,
                             psi_k_laststep.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             psi_k.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = psi_k.data<std::complex<double>>()[i * pv->ncol + j].real();
                bb = psi_k.data<std::complex<double>>()[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                    aa = 0.0;
                if (std::abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " psi_k_laststep:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa, bb;
                aa = psi_k_laststep.data<std::complex<double>>()[i * pv->ncol + j].real();
                bb = psi_k_laststep.data<std::complex<double>>()[i * pv->ncol + j].imag();
                if (std::abs(aa) < 1e-8)
                    aa = 0.0;
                if (std::abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
}

#endif
} // namespace module_tddft
