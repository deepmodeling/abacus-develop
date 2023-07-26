#include "middle_hamilt.h"

#include <complex>
#include <iostream>

#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"

namespace module_tddft
{
#ifdef __MPI

void half_Hmatrix(const Parallel_Orbitals* pv,
                  const int nband,
                  const int nlocal,
                  std::complex<double>* Htmp,
                  std::complex<double>* Stmp,
                  const std::complex<double>* H_laststep,
                  const std::complex<double>* S_laststep,
                  const int print_matrix)
{
    if (print_matrix)
    {
        GlobalV::ofs_running << std::setprecision(10);
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H(t+dt) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << Htmp[i * pv->ncol + j].real() << "+" << Htmp[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H(t):" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << H_laststep[i * pv->ncol + j].real() << "+"
                                     << H_laststep[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    std::complex<double> alpha = {0.5, 0.0};
    std::complex<double> beta = {0.5, 0.0};
    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, H_laststep, 1, 1, pv->desc, beta, Htmp, 1, 1, pv->desc);
    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, S_laststep, 1, 1, pv->desc, beta, Stmp, 1, 1, pv->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H (t+dt/2) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << Htmp[i * pv->ncol + j].real() << "+" << Htmp[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
}
#endif
} // namespace module_tddft