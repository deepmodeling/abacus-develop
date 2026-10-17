#include "source_base/constants.h"
#include "lr_util.h"
#include "source_base/parallel_reduce.h"
#include "source_cell/unitcell.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/module_external/scalapack_connector.h"
namespace LR_Util
{
    int spin_population_difference(const ModuleBase::matrix& wg, const int kpar, const int rank_in_pool)
    {   // only for nspin=2
        const int nk = wg.nr / 2;
        auto occ_sum_k = [&wg, nk](const int is, const int ib) -> double
        {
            double occupation = 0.0;
            for (int ik = 0; ik < nk; ++ik)
            {
                const int spin_k = is * nk + ik;
                occupation += wg(spin_k, ib);
            }
            return occupation;
        };
        // Sum the occupations of each channel FIRST and round once, instead of rounding band by band
        // and summing the differences. A half-occupied degenerate frontier pair (OH's 2-Pi doublet
        // smears its odd electron as 0.5/0.5 over the two pi_down orbitals) otherwise makes the answer
        // a coin flip: the stored values are 0.5000000052 and 0.4999999947, so one rounds up and one
        // down, and which way they land is pure noise.
        double up = 0.0;
        double dn = 0.0;
        for (int ib = 0; ib < wg.nc; ++ib)
        {
            up += occ_sum_k(0, ib);
            dn += occ_sum_k(1, ib);
        }
        // wg is replicated within a pool, but each pool holds different k points.
        if (kpar > 1)
        {
            if (rank_in_pool != 0)
            {
                up = 0.0;
                dn = 0.0;
            }
            Parallel_Reduce::reduce_all(up);
            Parallel_Reduce::reduce_all(dn);
        }
        return static_cast<int>(std::lround(up) - std::lround(dn));
    }

    std::string force_window_warning(const int basis_size, const int gs_bands,
        const int skipped_occupied, const std::vector<int>& x_virtual,
        const std::vector<int>& z_virtual)
    {
        assert(x_virtual.size() == z_virtual.size());
        const bool occupied_truncated = skipped_occupied > 0;
        const bool bands_truncated = gs_bands < basis_size;
        bool virtual_truncated = false;
        for (size_t spin = 0; spin < x_virtual.size(); ++spin)
        {
            virtual_truncated = virtual_truncated || x_virtual[spin] < z_virtual[spin];
        }
        if (!occupied_truncated && !bands_truncated && !virtual_truncated)
        {
            return "";
        }
        std::string message = " WARNING: truncated LR force window:";
        if (occupied_truncated)
        {
            message += " occupied X window is truncated;";
        }
        if (virtual_truncated)
        {
            message += " virtual X window is truncated;";
        }
        if (bands_truncated)
        {
            message += " ground-state band space is incomplete;";
        }
        message += " zero-padding X into the available Z virtual space is retained, but does not"
            " restore all active/inactive orbital-response terms. The force may differ from finite"
            " differences of the selected truncated excitation energy. For reliable forces and"
            " relaxation, use all occupied orbitals, nbands = NLOCAL (";
        message += std::to_string(basis_size);
        message += ") and the full spin-resolved virtual window (nvirt = nbands - nocc)."
            " Check force/window convergence before relaxation.";
        return message;
    }

    void check_force_pp(const UnitCell& cell, const bool cal_force, const std::string& calculation)
    {
        const bool relaxation = calculation == "relax" || calculation == "cell-relax";
        if (!cal_force && !relaxation)
        {
            return;
        }
        for (int type = 0; type < cell.ntype; ++type)
        {
            if (cell.atoms[type].ncpp.nlcc)
            {
                throw std::invalid_argument(
                    "LR excited-state forces and relaxation do not support NLCC pseudopotentials: element "
                    + cell.atoms[type].label + ", file " + cell.pseudo_fn.at(type)
                    + ". Core-density response-force and XC-kernel derivative terms are not implemented."
                      " LR spectra remain available with cal_force=0.");
            }
        }
    }

    const std::vector<double>& prepare_xc_sigma(const std::vector<double>& sigma,
                                               const bool is_hse06,
                                               std::vector<double>& hse_buffer)
    {
        if (!is_hse06)
        {
            return sigma;
        }
        hse_buffer.resize(sigma.size());
        for (std::size_t index = 0; index < sigma.size(); ++index)
        {
            hse_buffer[index] = std::max(sigma[index], 1e-6);
        }
        return hse_buffer;
    }

    /// =================PHYSICS====================
    int cal_nocc(int nelec) { return nelec / ModuleBase::DEGSPIN + nelec % static_cast<int>(ModuleBase::DEGSPIN); }

    int cal_nocc(double nelec, int nspin, int nupdown)
    {
        const double polarization = nspin == 2 ? std::abs(nupdown) : 0.0;
        if (!std::isfinite(nelec) || nelec <= 0.0 || polarization > nelec)
        {
            throw std::invalid_argument("LR: invalid electron number or spin population");
        }
        // Use ceil to include a partially occupied frontier orbital, without promoting
        // a numerically noisy integer population to the next orbital.
        const double population = (nelec + polarization) / 2.0;
        const double nearest_integer = std::round(population);
        const double stable_population = std::abs(population - nearest_integer) < 1e-8
                                             ? nearest_integer : population;
        return static_cast<int>(std::ceil(stable_population));
    }

    int cal_nocc_window(int requested_nocc, int nocc_max)
    {
        if (nocc_max <= 0)
        {
            throw std::invalid_argument("LR: no occupied orbitals");
        }
        return requested_nocc <= 0 ? nocc_max : std::min(requested_nocc, nocc_max);
    }

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
        ScalapackConnector::tranu(n,
            n,
            alpha,
            in,
            i1,
            i1,
            pmat.desc,
            beta,
            out,
            i1,
            i1,
            pmat.desc);
    }
    template<>
    void matsym<double>(double* inout, const int n, const Parallel_2D& pmat)
    {
        std::vector<double> tmp(pmat.get_local_size());
        std::copy(inout, inout + pmat.get_local_size(), tmp.begin());
        const double alpha = 0.5, beta = 0.5;
        const int i1 = 1;
        ScalapackConnector::tranu(n,
            n,
            alpha,
            tmp.data(),
            i1,
            i1,
            pmat.desc,
            beta,
            inout,
            i1,
            i1,
            pmat.desc);
    }
    template<>
    void matsym<std::complex<double>>(const std::complex<double>* in, const int n, const Parallel_2D& pmat, std::complex<double>* out)
    {
        std::copy(in, in + pmat.get_local_size(), out);
        const std::complex<double> alpha(0.5, 0.0), beta(0.5, 0.0);
        const int i1 = 1;
        ScalapackConnector::tranc(n,
            n,
            alpha,
            in,
            i1,
            i1,
            pmat.desc,
            beta,
            out,
            i1,
            i1,
            pmat.desc);
    }
    template<>
    void matsym<std::complex<double>>(std::complex<double>* inout, const int n, const Parallel_2D& pmat)
    {
        std::vector<std::complex<double>> tmp(pmat.get_local_size());
        std::copy(inout, inout + pmat.get_local_size(), tmp.begin());
        const std::complex<double> alpha(0.5, 0.0), beta(0.5, 0.0);
        const int i1 = 1;
        ScalapackConnector::tranc(n,
            n,
            alpha,
            tmp.data(),
            i1,
            i1,
            pmat.desc,
            beta,
            inout,
            i1,
            i1,
            pmat.desc);
    }

    template<>
    void mattrans<double>(const double* in, const int n, const Parallel_2D& pmat, double* out)
    {
        std::copy(in, in + pmat.get_local_size(), out);
        const double alpha = 1.0;
        const double beta = 0.0;
        const int i1 = 1;
        ScalapackConnector::tranu(n,
            n,
            alpha,
            in,
            i1,
            i1,
            pmat.desc,
            beta,
            out,
            i1,
            i1,
            pmat.desc);
    }
    template<>
    void mattrans<double>(double* inout, const int n, const Parallel_2D& pmat)
    {
        std::vector<double> tmp(pmat.get_local_size());
        std::copy(inout, inout + pmat.get_local_size(), tmp.begin());
        const double alpha = 1.0;
        const double beta = 0.0;
        const int i1 = 1;
        ScalapackConnector::tranu(n,
            n,
            alpha,
            tmp.data(),
            i1,
            i1,
            pmat.desc,
            beta,
            inout,
            i1,
            i1,
            pmat.desc);
    }
    template<>
    void mattrans<std::complex<double>>(const std::complex<double>* in, const int n, const Parallel_2D& pmat, std::complex<double>* out)
    {
        std::copy(in, in + pmat.get_local_size(), out);
        const std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);
        const int i1 = 1;
        ScalapackConnector::tranc(n,
            n,
            alpha,
            in,
            i1,
            i1,
            pmat.desc,
            beta,
            out,
            i1,
            i1,
            pmat.desc);
    }
    template<>
    void mattrans<std::complex<double>>(std::complex<double>* inout, const int n, const Parallel_2D& pmat)
    {
        std::vector<std::complex<double>> tmp(pmat.get_local_size());
        std::copy(inout, inout + pmat.get_local_size(), tmp.begin());
        const std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);
        const int i1 = 1;
        ScalapackConnector::tranc(n,
            n,
            alpha,
            tmp.data(),
            i1,
            i1,
            pmat.desc,
            beta,
            inout,
            i1,
            i1,
            pmat.desc);
    }


    template<>
    void matantisym<double>(double* inout, const int n, const Parallel_2D& pmat)
    {
        std::vector<double> tmp(pmat.get_local_size());
        std::copy(inout, inout + pmat.get_local_size(), tmp.begin());
        const double alpha = -0.5;
        const double beta = 0.5;
        const int i1 = 1;
        ScalapackConnector::tranu(n,
            n,
            alpha,
            tmp.data(),
            i1,
            i1,
            pmat.desc,
            beta,
            inout,
            i1,
            i1,
            pmat.desc);
    }
    template<>
    void matantisym<std::complex<double>>(std::complex<double>* inout, const int n, const Parallel_2D& pmat)
    {
        std::vector<std::complex<double>> tmp(pmat.get_local_size());
        std::copy(inout, inout + pmat.get_local_size(), tmp.begin());
        const std::complex<double> alpha(-0.5, 0.0);
        const std::complex<double> beta(0.5, 0.0);
        const int i1 = 1;
        ScalapackConnector::tranc(n,
            n,
            alpha,
            tmp.data(),
            i1,
            i1,
            pmat.desc,
            beta,
            inout,
            i1,
            i1,
            pmat.desc);
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

    template<>
    void diag_lapack<double>(const int& n, double* mat, double* eig)
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack<double>");
        int info = 0;
        char jobz = 'V', uplo = 'U';
        double work_tmp = 0.0;
        const int minus_one = -1;
        container::lapackConnector::heev(jobz, uplo, n, mat, n, eig, &work_tmp, minus_one, info);		// get best lwork
        const int lwork = work_tmp;
        std::vector<double> work2(lwork);
        container::lapackConnector::heev(jobz, uplo, n, mat, n, eig, work2.data(), lwork, info);
        if (info) { std::cout << "ERROR: Lapack solver, info=" << info << std::endl; }
    }
    template<>
    void diag_lapack<std::complex<double>>(const int& n, std::complex<double>* mat, double* eig)
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack <std::complex<double>>");
        int lwork = 2 * n;
        std::vector<std::complex<double>> work2(lwork);
        std::vector<double> rwork(3 * n - 2);
        int info = 0;
        char jobz = 'V', uplo = 'U';
        container::lapackConnector::heev(jobz, uplo, n, mat, n, eig, work2.data(), lwork, rwork.data(), info);
        if (info) { std::cout << "ERROR: Lapack solver, info=" << info << std::endl; }
    }
    template<>
    void diag_lapack_nh<double>(const int& n, double* mat, std::complex<double>* eig)
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack_nh<double>");
        int info = 0;
        char jobvl = 'N', jobvr = 'V';  //calculate right eigenvectors
        double work_tmp = 0.0;
        const int minus_one = -1;
        std::vector<double> eig_real(n);
        std::vector<double> eig_imag(n);
        const int ldvl = 1, ldvr = n;
        std::vector<double> vl(ldvl * n), vr(ldvr * n);
        container::lapackConnector::geev(jobvl, jobvr, n, mat, n, eig_real.data(), eig_imag.data(), vl.data(), ldvl, vr.data(), ldvr, &work_tmp, minus_one /*lwork*/, info);		// get best lwork
        const int lwork = work_tmp;
        std::vector<double> work2(lwork);
        container::lapackConnector::geev(jobvl, jobvr, n, mat, n, eig_real.data(), eig_imag.data(), vl.data(), ldvl, vr.data(), ldvr, work2.data(), lwork, info);
        if (info) { std::cout << "ERROR: Lapack solver dgeev, info=" << info << std::endl; }
        for (int i = 0;i < n;++i) { eig[i] = std::complex<double>(eig_real[i], eig_imag[i]); }
    }
    template<>
    void diag_lapack_nh<std::complex<double>>(const int& n, std::complex<double>* mat, std::complex<double>* eig)
    {
        ModuleBase::TITLE("LR_Util", "diag_lapack_nh <std::complex<double>>");
        int lwork = 2 * n;
        std::vector<std::complex<double>> work2(lwork);
        std::vector<double> rwork(3 * n - 2);
        int info = 0;
        char jobvl = 'N', jobvr = 'V';
        const int ldvl = 1, ldvr = n;
        std::vector<std::complex<double>> vl(ldvl * n), vr(ldvr * n);
        container::lapackConnector::geev(jobvl, jobvr, n, mat, n, eig, vl.data(), ldvl, vr.data(), ldvr, work2.data(), lwork, rwork.data(), info);
        if (info) { std::cout << "ERROR: Lapack solver zgeev, info=" << info << std::endl; }
    }

    template<>
    int lapack_linear_solver<double>(const double* A, double* x, const double* b, const int n, const int nrhs)
    {
        ModuleBase::TITLE("LR_Util", "lapack_linear_solver<double>");
        // 1. copy A to a mutable array
        std::vector<double> A_copy(A, A + n * n);
        // copy b to x
        std::copy(b, b + n * nrhs, x);
        // 2. LU decomposition: A->LU
        std::vector<int> ipiv(n);   // pivot indices
        int info = 0;
        container::lapackConnector::getrf(n, n, A_copy.data(), n, ipiv.data(), info);
        if (info)
        {
            std::cout << "ERROR: Lapack solver dgetrf, info=" << info << std::endl;
        }
        // 3. Solve Ax=b
        const char trans = 'N';
        container::lapackConnector::getrs(trans, n, nrhs, A_copy.data(), n, ipiv.data(), x, n, info);
        if (info)
        {
            std::cout << "ERROR: Lapack solver dgetrs, info=" << info << std::endl;
        }
        return info;
    }

    template<>
    int lapack_linear_solver<std::complex<double>>(const std::complex<double>* A, std::complex<double>* x, const std::complex<double>* b, const int n, const int nrhs)
    {
        ModuleBase::TITLE("LR_Util", "lapack_linear_solver<complex<double>>");
        // 1. copy A to a mutable array
        std::vector<std::complex<double>> A_copy(A, A + n * n);
        // copy b to x
        std::copy(b, b + n * nrhs, x);
        // 2. LU decomposition: A->LU
        std::vector<int> ipiv(n);   // pivot indices, for
        int info = 0;
        container::lapackConnector::getrf(n, n, A_copy.data(), n, ipiv.data(), info);
        if (info)
        {
            std::cout << "ERROR: Lapack solver zgetrf, info=" << info << std::endl;
        }
        // 3. Solve Ax=b
        const char trans = 'N';
        container::lapackConnector::getrs(trans, n, nrhs, A_copy.data(), n, ipiv.data(), x, n, info);
        if (info)
        {
            std::cout << "ERROR: Lapack solver zgetrs, info=" << info << std::endl;
        }
        return info;
    }

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
