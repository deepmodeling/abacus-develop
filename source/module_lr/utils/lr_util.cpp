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
        for (int i = 0;i < pmat.get_local_size();++i) { out[i] = in[i]; }
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
        for (int i = 0;i < pmat.get_local_size();++i) { out[i] = in[i]; }
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
    container::Tensor mat2ten_double(ModuleBase::matrix& m)
    {
        container::Tensor t(DAT::DT_DOUBLE, DEV::CpuDevice, { m.nr, m.nc });
        for (int i = 0;i < t.NumElements();++i) {t.data<double>()[i] = m.c[i];
}
        return t;
    }
    std::vector<container::Tensor> mat2ten_double(std::vector<ModuleBase::matrix>& m)
    {
        std::vector<container::Tensor> t;
        for (int i = 0;i < m.size();++i) { t.push_back(mat2ten_double(m[i]));
}
        return t;
    }
    ModuleBase::matrix ten2mat_double(container::Tensor& t)
    {
        ModuleBase::matrix m(t.shape().dims()[0], t.shape().dims()[1]);
        for (int i = 0;i < t.NumElements();++i) {m.c[i] = t.data<double>()[i];
}
        return m;
    }
    std::vector<ModuleBase::matrix> ten2mat_double(std::vector<container::Tensor>& t)
    {
        std::vector<ModuleBase::matrix> m;
        for (int i = 0;i < t.size();++i) { m.push_back(ten2mat_double(t[i]));
}
        return m;
    }
    container::Tensor mat2ten_complex(ModuleBase::ComplexMatrix& m)
    {
        container::Tensor t(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { m.nr, m.nc });
        for (int i = 0;i < t.NumElements();++i) {t.data<std::complex<double>>()[i] = m.c[i];
}
        return t;
    }
    std::vector<container::Tensor> mat2ten_complex(std::vector<ModuleBase::ComplexMatrix>& m)
    {
        std::vector<container::Tensor> t;
        for (int i = 0;i < m.size();++i) { t.push_back(mat2ten_complex(m[i]));
}
        return t;
    }
    ModuleBase::ComplexMatrix ten2mat_complex(container::Tensor& t)
    {
        ModuleBase::ComplexMatrix m(t.shape().dims()[0], t.shape().dims()[1]);
        for (int i = 0;i < t.NumElements();++i) {m.c[i] = t.data<std::complex<double>>()[i];
}
        return m;
    }
    std::vector<ModuleBase::ComplexMatrix> ten2mat_complex(std::vector<container::Tensor>& t)
    {
        std::vector<ModuleBase::ComplexMatrix> m;
        for (int i = 0;i < t.size();++i) { m.push_back(ten2mat_complex(t[i]));
}
        return m;
    }

    ModuleBase::matrix vec2mat(const std::vector<double>& v, const int nr, const int nc)
    {
        assert(v.size() == nr * nc);
        ModuleBase::matrix m(nr, nc, false);
        for (int i = 0;i < v.size();++i) { m.c[i] = v[i];
}
        return m;
    }
    ModuleBase::ComplexMatrix vec2mat(const std::vector<std::complex<double>>& v, const int nr, const int nc)
    {
        assert(v.size() == nr * nc);
        ModuleBase::ComplexMatrix m(nr, nc, false);
        for (int i = 0;i < v.size();++i) { m.c[i] = v[i];
}
        return m;
    }
    std::vector<ModuleBase::matrix> vec2mat(const std::vector<std::vector<double>>& v, const int nr, const int nc)
    {
        std::vector<ModuleBase::matrix> m(v.size());
        for (int i = 0;i < v.size();++i) { m[i] = vec2mat(v[i], nr, nc);
}
        return m;
    }
    std::vector<ModuleBase::ComplexMatrix> vec2mat(const std::vector<std::vector<std::complex<double>>>& v, const int nr, const int nc)
    {
        std::vector<ModuleBase::ComplexMatrix> m(v.size());
        for (int i = 0;i < v.size();++i) { m[i] = vec2mat(v[i], nr, nc);
}
        return m;
    }

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
#ifdef USE_LIBXC
    void grad(const double* rhor,
        ModuleBase::Vector3<double>* gdr,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba)
    {
        std::vector<std::complex<double>> rhog(rho_basis.npw);
        rho_basis.real2recip(rhor, rhog.data());
        XC_Functional::grad_rho(rhog.data(), gdr, &rho_basis, tpiba);
    }
    void laplace(const double* rhor, double* lapn,
        const ModulePW::PW_Basis& rho_basis,
        const double& tpiba2)
    {
        ModuleBase::GlobalFunc::ZEROS(lapn, rho_basis.nrxx);
        std::vector<std::complex<double>> rhog(rho_basis.npw);
        std::vector<double> tmp_rhor(rho_basis.nrxx);
        rho_basis.real2recip(rhor, rhog.data());
        for (int i = 0;i < 3;++i)
        {
            for (int ig = 0; ig < rho_basis.npw; ig++) {
                rhog[ig] *= pow(rho_basis.gcar[ig][i], 2);
}
            rho_basis.recip2real(rhog.data(), tmp_rhor.data());
            for (int ir = 0; ir < rho_basis.nrxx; ir++) {
                lapn[ir] -= tmp_rhor[ir] * tpiba2;
}
        }
    }
#endif
}
