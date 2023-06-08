#include "propagator.h"

#include <complex>
#include <iostream>

#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_io/input.h"

namespace module_tddft
{
Propagator::~Propagator()
{
}

#ifdef __MPI

inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

void Propagator::compute_propagator(const int nlocal,
                                    const std::complex<double>* Stmp,
                                    const std::complex<double>* Htmp,
                                    const std::complex<double>* H_laststep,
                                    std::complex<double>* U_operator,
                                    const int print_matrix) const
{
    int tag;
    switch (ptype)
    {
    case 0:
        compute_propagator_cn2(nlocal, Stmp, Htmp, U_operator, print_matrix);
        break;

    case 1:
        tag = 1;
        compute_propagator_taylor(nlocal, Stmp, Htmp, U_operator, print_matrix, tag);
        break;

    case 2:
        compute_propagator_etrs(nlocal, Stmp, Htmp, H_laststep, U_operator, print_matrix);

        break;

    default:
        std::cout << "method of propagator is wrong" << std::endl;
        break;
    }
}

void Propagator::compute_propagator_cn2(const int nlocal,
                                        const std::complex<double>* Stmp,
                                        const std::complex<double>* Htmp,
                                        std::complex<double>* U_operator,
                                        const int print_matrix) const
{
    // (1) copy Htmp to Numerator & Denominator
    std::complex<double>* Numerator = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Numerator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Numerator, 1);

    std::complex<double>* Denominator = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Denominator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Denominator, 1);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Stmp[i * this->ParaV->ncol + j].real() << "+"
                                     << Stmp[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Numerator[i * this->ParaV->ncol + j].real() << "+"
                                     << Numerator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute Numerator & Denominator by GEADD
    // Numerator = Stmp - i*para * Htmp;     beta1 = - para = -0.25 * INPUT.mdp.md_dt
    // Denominator = Stmp + i*para * Htmp;   beta2 = para = 0.25 * INPUT.mdp.md_dt
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta1 = {0.0, -0.25 * INPUT.mdp.md_dt};
    std::complex<double> beta2 = {0.0, 0.25 * INPUT.mdp.md_dt};

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp,
                              1,
                              1,
                              this->ParaV->desc,
                              beta1,
                              Numerator,
                              1,
                              1,
                              this->ParaV->desc);
    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp,
                              1,
                              1,
                              this->ParaV->desc,
                              beta2,
                              Denominator,
                              1,
                              1,
                              this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " beta=" << beta1 << std::endl;
        GlobalV::ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Denominator[i * this->ParaV->ncol + j].real() << "+"
                                     << Denominator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) Next, invert Denominator
    int* ipiv = new int[this->ParaV->nloc];
    int info = 0;
    // (3.1) compute ipiv
    ScalapackConnector::getrf(nlocal, nlocal, Denominator, 1, 1, this->ParaV->desc, ipiv, &info);
    int lwork = -1;
    int liwotk = -1;
    std::vector<std::complex<double>> work(1, 0);
    std::vector<int> iwork(1, 0);
    // (3.2) compute work
    ScalapackConnector::getri(nlocal,
                              Denominator,
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
    // (3.3) compute inverse matrix of Denominator
    ScalapackConnector::getri(nlocal,
                              Denominator,
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

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // (4) U_operator = Denominator * Numerator;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             Denominator,
                             1,
                             1,
                             this->ParaV->desc,
                             Numerator,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " fenmu^-1:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Denominator[i * this->ParaV->ncol + j].real() << "+"
                                     << Denominator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Numerator[i * this->ParaV->ncol + j].real() << "+"
                                     << Numerator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa, bb;
                aa = U_operator[i * this->ParaV->ncol + j].real();
                bb = U_operator[i * this->ParaV->ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
    }

    delete[] Numerator;
    delete[] Denominator;
    delete[] ipiv;
}

void Propagator::compute_propagator_taylor(const int nlocal,
                                           const std::complex<double>* Stmp,
                                           const std::complex<double>* Htmp,
                                           std::complex<double>* U_operator,
                                           const int print_matrix,
                                           const int tag) const
{
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
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Stmp[i * this->ParaV->ncol + j].real() << "+"
                                     << Stmp[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Htmp[i * this->ParaV->ncol + j].real() << "+"
                                     << Htmp[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    // set rank0
    int info;
    int myid;
    MPI_Comm_rank(this->ParaV->comm_2D, &myid);
    int naroc[2]; // maximum number of row or column

    for (int iprow = 0; iprow < this->ParaV->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < this->ParaV->dim1; ++ipcol)
        {
            const int coord[2] = {iprow, ipcol};
            int src_rank;
            info = MPI_Cart_rank(this->ParaV->comm_2D, coord, &src_rank);
            if (myid == src_rank)
            {
                naroc[0] = this->ParaV->nrow;
                naroc[1] = this->ParaV->ncol;
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, this->ParaV->nb, this->ParaV->dim1, ipcol);
                    if (igcol >= nlocal)
                        continue;
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, this->ParaV->nb, this->ParaV->dim0, iprow);
                        if (igrow >= nlocal)
                            continue;
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
    }     // loop iprow

    std::complex<double> beta = {0.0, -0.5 * INPUT.mdp.md_dt / tag}; // for ETRS tag=2 , for taylor tag=1

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
        GlobalV::ofs_running << " A_matrix:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << A_matrix[i * this->ParaV->ncol + j].real() << "+"
                                     << A_matrix[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa, bb;
                aa = U_operator[i * this->ParaV->ncol + j].real();
                bb = U_operator[i * this->ParaV->ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
    }

    delete[] rank0;
    delete[] rank2;
    delete[] rank3;
    delete[] rank4;
    delete[] tmp1;
    delete[] tmp2;
    delete[] ipiv;
}

void Propagator::compute_propagator_etrs(const int nlocal,
                                         const std::complex<double>* Stmp,
                                         const std::complex<double>* Htmp,
                                         const std::complex<double>* H_laststep,
                                         std::complex<double>* U_operator,
                                         const int print_matrix) const
{
    std::complex<double>* U1 = new std::complex<double>[this->ParaV->nloc];
    std::complex<double>* U2 = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(U1, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(U2, this->ParaV->nloc);
    int tag = 2;
    compute_propagator_taylor(nlocal, Stmp, Htmp, U1, print_matrix, tag);
    compute_propagator_taylor(nlocal, Stmp, H_laststep, U2, print_matrix, tag);
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             U1,
                             1,
                             1,
                             this->ParaV->desc,
                             U2,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc);
}

#endif
} // namespace module_tddft