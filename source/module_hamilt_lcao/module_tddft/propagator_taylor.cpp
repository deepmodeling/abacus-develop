#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/module_container/ATen/kernels/lapack.h"
#include "module_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "module_base/module_device/memory_op.h"              // memory operations
#include "module_base/scalapack_connector.h"
#include "module_parameter/parameter.h"
#include "propagator.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
#ifdef __MPI
void Propagator::compute_propagator_taylor(const int nlocal,
                                           const std::complex<double>* Stmp,
                                           const std::complex<double>* Htmp,
                                           std::complex<double>* U_operator,
                                           std::ofstream& ofs_running,
                                           const int print_matrix,
                                           const int tag) const
{
    assert(this->ParaV->nloc > 0);

    ModuleBase::GlobalFunc::ZEROS(U_operator, this->ParaV->nloc);
    std::complex<double>* A_matrix = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(A_matrix, this->ParaV->nloc);
    std::complex<double>* rank0 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(rank0, this->ParaV->nloc);
    std::complex<double>* rank2 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(rank2, this->ParaV->nloc);
    std::complex<double>* rank3 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(rank3, this->ParaV->nloc);
    std::complex<double>* rank4 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(rank4, this->ParaV->nloc);
    std::complex<double>* tmp1 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, this->ParaV->nloc);
    std::complex<double>* tmp2 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(tmp2, this->ParaV->nloc);
    std::complex<double>* Sinv = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Sinv, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Stmp, 1, Sinv, 1);

    if (print_matrix)
    {
        ofs_running << std::endl;
        ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Stmp[in + j].real() << "+" << Stmp[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Htmp[in + j].real() << "+" << Htmp[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    // set rank0
    int info = 0;
    int naroc[2] = {0, 0}; // maximum number of row or column

    for (int iprow = 0; iprow < this->ParaV->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < this->ParaV->dim1; ++ipcol)
        {
            if (iprow == ParaV->coord[0] && ipcol == ParaV->coord[1])
            {
                naroc[0] = this->ParaV->nrow;
                naroc[1] = this->ParaV->ncol;
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, this->ParaV->nb, this->ParaV->dim1, ipcol);
                    if (igcol >= nlocal)
                    {
                        continue;
                    }
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, this->ParaV->nb, this->ParaV->dim0, iprow);
                        if (igrow >= nlocal)
                        {
                            continue;
                        }
                        if (igcol == igrow)
                        {
                            rank0[j * naroc[0] + i] = {1.0, 0.0};
                        }
                        else
                        {
                            rank0[j * naroc[0] + i] = {0.0, 0.0};
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow

    std::complex<double> beta = {0.0, -0.5 * this->dt / tag}; // for ETRS tag=2 , for taylor tag=1

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // invert Stmp
    int* ipiv = new int[this->ParaV->nloc];
    // (3.1) compute ipiv
    ScalapackConnector::getrf(nlocal, nlocal, Sinv, 1, 1, this->ParaV->desc, ipiv, &info);
    int lwork = -1;
    int liwotk = -1;
    std::vector<std::complex<double>> work(1, 0);
    std::vector<int> iwork(1, 0);
    // (3.2) compute work
    ScalapackConnector::getri(nlocal,
                              Sinv,
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv,
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwotk,
                              &info);
    lwork = work[0].real();
    work.resize(lwork, 0);
    liwotk = iwork[0];
    iwork.resize(liwotk, 0);
    ScalapackConnector::getri(nlocal,
                              Sinv,
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv,
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwotk,
                              &info);
    assert(0 == info);

    //  A_matrix = - idt S^-1 H ;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             beta,
                             Sinv,
                             1,
                             1,
                             this->ParaV->desc,
                             Htmp,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc);

    //  rank2 = A^2 ;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             rank2,
                             1,
                             1,
                             this->ParaV->desc);

    //  rank3 = A^3 ;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc,
                             rank2,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             rank3,
                             1,
                             1,
                             this->ParaV->desc);

    //  rank4 = A^4 ;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc,
                             rank3,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             rank4,
                             1,
                             1,
                             this->ParaV->desc);

    std::complex<double> p1 = {1.0, 0.0};
    std::complex<double> p2 = {1.0 / 2.0, 0.0};
    std::complex<double> p3 = {1.0 / 6.0, 0.0};
    std::complex<double> p4 = {1.0 / 24.0, 0.0};

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              p1,
                              rank0,
                              1,
                              1,
                              this->ParaV->desc,
                              p1,
                              U_operator,
                              1,
                              1,
                              this->ParaV->desc);

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              p2,
                              rank2,
                              1,
                              1,
                              this->ParaV->desc,
                              p1,
                              U_operator,
                              1,
                              1,
                              this->ParaV->desc);

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              p3,
                              rank3,
                              1,
                              1,
                              this->ParaV->desc,
                              p1,
                              U_operator,
                              1,
                              1,
                              this->ParaV->desc);

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              p4,
                              rank4,
                              1,
                              1,
                              this->ParaV->desc,
                              p1,
                              U_operator,
                              1,
                              1,
                              this->ParaV->desc);

    if (print_matrix)
    {
        ofs_running << " A_matrix:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << A_matrix[in + j].real() << "+" << A_matrix[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa = U_operator[in + j].real();
                double bb = U_operator[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
    }
    delete[] A_matrix;
    delete[] rank0;
    delete[] rank2;
    delete[] rank3;
    delete[] rank4;
    delete[] tmp1;
    delete[] tmp2;
    delete[] Sinv;
    delete[] ipiv;
}
#endif // __MPI
} // namespace module_tddft
