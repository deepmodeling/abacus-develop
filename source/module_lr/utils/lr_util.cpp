#include "module_base/constants.h"
#include "lr_util.h"
#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
namespace LR_Util
{
    /// =================PHYSICS====================
    int cal_nocc(int nelec) { return nelec / ModuleBase::DEGSPIN + nelec % static_cast<int>(ModuleBase::DEGSPIN); }

    std::pair<ModuleBase::matrix, std::vector<std::pair<int, int>>>
        set_ix_map_diagonal(bool mode, int nocc, int nvirt)
    {
        int npairs = nocc * nvirt;
        ModuleBase::matrix ioiv2ix(nocc, nvirt, true);
        std::vector<std::pair<int, int>> ix2ioiv(npairs);
        int io = nocc - 1, iv = 0;    //start：leftup
        if (mode == 0)  // leftdown->rightup
        {
            for (int ix = 0;ix < npairs - 1;++ix)
            {
                // 1. set value
                ioiv2ix(io, iv) = ix;
                ix2ioiv[ix] = std::make_pair(io, iv);
                // 2. move
                if (io == nocc - 1 || iv == nvirt - 1)    // rightup bound
                {
                    int io_next = std::max(nocc - iv - 1 - (nocc - io), 0);
                    iv -= (io - io_next) - 1;
                    io = io_next;
                }
                else { ++io;++iv; }//move rightup
            }
        }
        else    //rightup->leftdown
        {
            for (int ix = 0;ix < npairs - 1;++ix)
            {
                // 1. set value
                ioiv2ix(io, iv) = ix;
                ix2ioiv[ix] = std::make_pair(io, iv);
                // 2. move
                if (io == 0 || iv == 0)    // leftdown bound
                {
                    int iv_next = std::min(nocc - io + iv, nvirt - 1);
                    io += (iv_next - iv) - 1;
                    iv = iv_next;
                }
                else { --iv;--io; }//move leftdown
            }
        }
        //final set: rightdown
        assert(io == 0);
        assert(iv == nvirt - 1);
        ioiv2ix(io, iv) = npairs - 1;
        ix2ioiv[npairs - 1] = std::make_pair(io, iv);
        return std::make_pair(std::move(ioiv2ix), std::move(ix2ioiv));
    }

    /// =================ALGORITHM====================

#ifdef __MPI
    template<>
    void matsym<double>(const double* in, const int n, const Parallel_2D& pmat, double* out)
    {
        std::copy(in, in + pmat.get_local_size(), out);
        const double alpha = 0.5, beta = 0.5;
        const int i1 = 1;
        pdtran_(&n, &n, &alpha, in, &i1, &i1, pmat.desc, &beta, out, &i1, &i1, pmat.desc);
    }
    template<>
    void matsym<double>(double* inout, const int n, const Parallel_2D& pmat)
    {
        std::vector<double> tmp(pmat.get_local_size());
        std::copy(inout, inout + pmat.get_local_size(), tmp.begin());
        const double alpha = 0.5, beta = 0.5;
        const int i1 = 1;
        pdtran_(&n, &n, &alpha, tmp.data(), &i1, &i1, pmat.desc, &beta, inout, &i1, &i1, pmat.desc);
    }
    template<>
    void matsym<std::complex<double>>(const std::complex<double>* in, const int n, const Parallel_2D& pmat, std::complex<double>* out)
    {
        std::copy(in, in + pmat.get_local_size(), out);
        const std::complex<double> alpha(0.5, 0.0), beta(0.5, 0.0);
        const int i1 = 1;
        pztranc_(&n, &n, &alpha, in, &i1, &i1, pmat.desc, &beta, out, &i1, &i1, pmat.desc);
    }
    template<>
    void matsym<std::complex<double>>(std::complex<double>* inout, const int n, const Parallel_2D& pmat)
    {
        std::vector<std::complex<double>> tmp(pmat.get_local_size());
        std::copy(inout, inout + pmat.get_local_size(), tmp.begin());
        const std::complex<double> alpha(0.5, 0.0), beta(0.5, 0.0);
        const int i1 = 1;
        pztranc_(&n, &n, &alpha, tmp.data(), &i1, &i1, pmat.desc, &beta, inout, &i1, &i1, pmat.desc);
    }
#endif

    // for the first matrix in the commutator
    void setup_2d_division(Parallel_2D& pv, int nb, int gr, int gc)
    {
        ModuleBase::TITLE("LR_Util", "setup_2d_division");
#ifdef __MPI
        pv.init(gr, gc, nb, MPI_COMM_WORLD);
#else
        pv.set_serial(gr, gc);
#endif
    }

#ifdef __MPI
    // for the other matrices in the commutator other than the first one
    void setup_2d_division(Parallel_2D& pv, int nb, int gr, int gc, const int& blacs_ctxt_in)
    {
        ModuleBase::TITLE("LR_Util", "setup_2d_division");
        pv.set(gr, gc, nb, blacs_ctxt_in);
    }
#endif

    void diag_lapack(const int& n, double* mat, double* eig)
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack<double>");
        int info = 0;
        char jobz = 'V', uplo = 'U';
        double work_tmp;
        const int minus_one = -1;
        dsyev_(&jobz, &uplo, &n, mat, &n, eig, &work_tmp, &minus_one, &info);		// get best lwork
        const int lwork = work_tmp;
        double* work2 = new double[lwork];
        dsyev_(&jobz, &uplo, &n, mat, &n, eig, work2, &lwork, &info);
        if (info) { std::cout << "ERROR: Lapack solver, info=" << info << std::endl; }
        delete[] work2;
    }

    void diag_lapack(const int& n, std::complex<double>* mat, double* eig)
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack<complex<double>>");
        int lwork = 2 * n;
        std::complex<double>* work2 = new std::complex<double>[lwork];
        double* rwork = new double[3 * n - 2];
        int info = 0;
        char jobz = 'V', uplo = 'U';
        zheev_(&jobz, &uplo, &n, mat, &n, eig, work2, &lwork, rwork, &info);
        if (info) { std::cout << "ERROR: Lapack solver, info=" << info << std::endl; }
        delete[] rwork;
        delete[] work2;
    }

    void diag_lapack_nh(const int& n, double* mat, std::complex<double>* eig)
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack_nh<double>");
        int info = 0;
        char jobvl = 'N', jobvr = 'V';  //calculate right eigenvectors
        double work_tmp;
        const int minus_one = -1;
        std::vector<double> eig_real(n);
        std::vector<double> eig_imag(n);
        const int ldvl = 1, ldvr = n;
        std::vector<double> vl(ldvl * n), vr(ldvr * n);
        dgeev_(&jobvl, &jobvr, &n, mat, &n, eig_real.data(), eig_imag.data(),
            vl.data(), &ldvl, vr.data(), &ldvr, &work_tmp, &minus_one /*lwork*/, &info);		// get best lwork
        const int lwork = work_tmp;
        std::vector<double> work2(lwork);
        dgeev_(&jobvl, &jobvr, &n, mat, &n, eig_real.data(), eig_imag.data(),
            vl.data(), &ldvl, vr.data(), &ldvr, work2.data(), &lwork, &info);
        if (info) { std::cout << "ERROR: Lapack solver dgeev, info=" << info << std::endl; }
        for (int i = 0;i < n;++i) { eig[i] = std::complex<double>(eig_real[i], eig_imag[i]); }
    }

    void diag_lapack_nh(const int& n, std::complex<double>* mat, std::complex<double>* eig)
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack_nh<complex<double>>");
        int lwork = 2 * n;
        std::vector<std::complex<double>> work2(lwork);
        std::vector<double> rwork(3 * n - 2);
        int info = 0;
        char jobvl = 'N', jobvr = 'V';
        const int ldvl = 1, ldvr = n;
        std::vector<std::complex<double>> vl(ldvl * n), vr(ldvr * n);
        zgeev_(&jobvl, &jobvr, &n, mat, &n, eig,
            vl.data(), &ldvl, vr.data(), &ldvr, work2.data(), &lwork, rwork.data(), &info);
        if (info) { std::cout << "ERROR: Lapack solver zgeev, info=" << info << std::endl; }
    }

#ifdef __MPI
    void diag_scalapack(const int& n, double* mat, double* eigval, double* eigvec, const int(&desc)[9])
    {
        ModuleBase::TITLE("LR_Util", "diag_scalapack<double>");
        const char jobz = 'V', uplo = 'U';
        const int minus_one = -1;
        const int i1 = 1;
        int info = 0;
        double lwork_tmp = 0.0;
        pdsyev_(&jobz, &uplo, &n,
            mat, &i1, &i1, desc,
            eigval, eigvec, &i1, &i1, desc,
            &lwork_tmp, &minus_one, &info);		// get the optimal size of work into lwork
        const int lwork = lwork_tmp;
        std::vector<double> work(lwork);
        pdsyev_(&jobz, &uplo, &n,
            mat, &i1, &i1, desc,
            eigval, eigvec, &i1, &i1, desc,
            work.data(), &lwork, &info);
        if (info) { std::cout << "ERROR: Scalapack solver, info=" << info << std::endl; }
    }
    void diag_scalapack(const int& n, std::complex<double>* mat, double* eigval, std::complex<double>* eigvec, const int(&desc)[9])
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack<complex<double>>");
        const char jobz = 'V', uplo = 'U';
        const int minus_one = -1;
        const int i1 = 1;
        int info = 0;
        std::vector<std::complex<double>> work(1, 0.0);
        std::vector<double>rwork(1, 0.0);
        // pzheev_(&jobz, &uplo, &n,
        //     mat, &i1, &i1, desc,
        //     eigval, eigvec, &i1, &i1, desc,
        //     work.data(), &minus_one, rwork.data(), &minus_one, &info);   // get the optimal workspace size
        /// try pzheevd
        // int liwork = 0;
        // pzheevd_(&jobz, &uplo, &n,
        //     mat, &i1, &i1, desc,
        //     eigval, eigvec, &i1, &i1, desc,
        //     &lwork_tmp, &minus_one, &lrwork_tmp, &minus_one, &liwork, &minus_one, &info);   // get the optimal workspace size

        // try pzheevx
        const char range = 'A';
        const double zero = 0.0;
        double abstol = 0.0;
        int nz = n;
        std::vector<int> iwork(1, 0);
        std::vector<int> ifail(n, 0);
        std::vector<int> iclustr(2 * GlobalV::DSIZE);
        std::vector<double> gap(GlobalV::DSIZE);
        pzheevx_(&jobz, &range, &uplo, &n,
            mat, &i1, &i1, desc,
            &zero, &zero, &i1, &i1, &zero,
            &nz, &nz, eigval, &zero,
            eigvec, &i1, &i1, desc,
            work.data(), &minus_one, rwork.data(), &minus_one, iwork.data(), &minus_one,
            ifail.data(), iclustr.data(), gap.data(), &info);

        const int lwork = work.at(0).real();
        work.resize(lwork);
        const int lrwork = rwork.at(0);
        rwork.resize(lrwork);
        const int liwork = iwork.at(0);
        iwork.resize(liwork);
        // std::cout << "pzheevx: query result: lwork=" << work.at(0) << ", lrwork=" << rwork.at(0) << ", liwork=" << iwork.at(0) << std::endl;

        // pzheev_(&jobz, &uplo, &n,
        //     mat, &i1, &i1, desc,
        //     eigval, eigvec, &i1, &i1, desc,
        //     work.data(), &lwork, rwork.data(), &lrwork, &info);
        // std::vector<int> iwork(liwork);
        // pzheevd_(&jobz, &uplo, &n,
        //     mat, &i1, &i1, desc,
        //     eigval, eigvec, &i1, &i1, desc,
        //     work.data(), &lwork, rwork.data(), &lrwork, iwork.data(), &liwork, &info);
        pzheevx_(&jobz, &range, &uplo, &n,
            mat, &i1, &i1, desc,
            &zero, &zero, &i1, &i1, &zero,
            &nz, &nz, eigval, &zero,
            eigvec, &i1, &i1, desc,
            work.data(), &lwork, rwork.data(), &lrwork, iwork.data(), &liwork,
            ifail.data(), iclustr.data(), gap.data(), &info);
        if (info) { std::cout << "ERROR: Scalapack solver, info=" << info << std::endl; }
    }

    void diag_scalapack(const int& n, std::complex<double>* hmat, std::complex<double>* const smat, double* eigval, std::complex<double>* eigvec, const int(&desc)[9])
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack<complex<double>>");
        const char jobz = 'V', uplo = 'U', range = 'A';
        int minus_one = -1;
        const int i1 = 1;
        const double zero = 0.0;
        int info = 0;
        double abstol = 0.0;
        int nz = n;
        std::vector<std::complex<double>> work(1, 0.0);
        std::vector<double>rwork(1, 0.0);
        std::vector<int> iwork(1, 0);
        std::vector<int> ifail(n, 0);
        std::vector<int> iclustr(2 * GlobalV::DSIZE);
        std::vector<double> gap(GlobalV::DSIZE);
        pzhegvx_(&i1, &jobz, &range, &uplo, &n,
            hmat, &i1, &i1, desc, smat, &i1, &i1, desc,
            &zero, &zero, &i1, &i1, &zero,
            &nz, &nz, eigval, &zero,
            eigvec, &i1, &i1, desc,
            work.data(), &minus_one, rwork.data(), &minus_one, iwork.data(), &minus_one,
            ifail.data(), iclustr.data(), gap.data(), &info);

        int lwork = work.at(0).real();
        work.resize(lwork);
        int lrwork = rwork.at(0);
        rwork.resize(lrwork);
        int liwork = iwork.at(0);
        iwork.resize(liwork);
        // std::cout << "pzhegvx: query result: lwork=" << work.at(0) << ", lrwork=" << rwork.at(0) << ", liwork=" << iwork.at(0) << std::endl;
        pzhegvx_(&i1, &jobz, &range, &uplo, &n,
            hmat, &i1, &i1, desc,
            smat, &i1, &i1, desc,
            &zero, &zero, &i1, &i1, &zero,
            &nz, &nz, eigval, &zero,
            eigvec, &i1, &i1, desc,
            work.data(), &lwork, rwork.data(), &lrwork, iwork.data(), &liwork,
            ifail.data(), iclustr.data(), gap.data(), &info);
        if (info) { std::cout << "ERROR: Scalapack solver, info=" << info << std::endl; }
    }
#endif

    std::string tolower(const std::string& str)
    {
        std::string str_lower = str;
        std::transform(str_lower.begin(), str_lower.end(), str_lower.begin(), ::tolower);
        return str_lower;
    }
    std::string toupper(const std::string& str)
    {
        std::string str_upper = str;
        std::transform(str_upper.begin(), str_upper.end(), str_upper.begin(), ::toupper);
        return str_upper;
    }
}