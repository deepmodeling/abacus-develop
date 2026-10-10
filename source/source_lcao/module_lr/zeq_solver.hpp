#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_ZEQ_SOLVER_HPP
#define ABACUS_SOURCE_LCAO_MODULE_LR_ZEQ_SOLVER_HPP
#include "gradient_inputs.h"
#include <fstream>
#include "zeq_solver.h"
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <vector>
#include "source_base/opt_cg.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "zeqlin_solv.h"

namespace LR
{
    // inline void solve_Z_CG(double* const Z, const double* const R, const int& ld, const int& nstates,
    //     std::function<void(const double* const, double* const)> f_LZ, const bool test_force)
    // Opt_CG's interfaces have no const qualifier for the pointer
    inline void solve_Z_CG(double* const Z, double* R, const int& ld, const int& nstates,
        std::function<void(const double* const, double* const)> f_LZ, const bool test_force)
    {
        ModuleBase::TITLE("Z_vector", "solve_Z_CG");
        ModuleBase::timer::start("Z_vector", "solve_Z_CG");
        const int maxiter = 100;
        double tol = 1e-6;
        double residual = 10.;
        const int size = nstates * ld;
        container::Tensor P = LR_Util::newTensor<double>({ nstates, ld });   // step length
        container::Tensor LP = LR_Util::newTensor<double>({ nstates, ld }); //f_LZ(P)
        ModuleBase::zeros(Z, size);

        ModuleBase::Opt_CG cg;
        cg.allocate(size);
        cg.init_b(R);
        std::cout << "Start solving Z-vector equaiton with CG method ..." << std::endl;
        for (int iter = 0; iter < maxiter; ++iter)
        {
            if (residual < tol)
            {
                break;
            }
            cg.next_direct(LP.data<double>(), 0, P.data<double>());
            std::cout << "iter=" << iter << " residual=" << cg.get_residual() << std::endl;
            // std::cout << "Z=" << std::endl;
            // LR_Util::print_value(P.data<double>(), nstates, ld);
            // std::cout << "LZ before=" << std::endl;
            // LR_Util::print_value(LP.data<double>(), nstates, ld);
            f_LZ(P.data<double>(), LP.data<double>());  // L: act each operators on P
            // std::cout << "LZ=" << std::endl;
            // LR_Util::print_value(LP.data<double>(), nstates, ld);
            int ifPD = 0;   //???
            double step = cg.step_length(LP.data<double>(), P.data<double>(), ifPD);
            // for (int i = 0; i < size; ++i) Z[i] += step * P[i];
            std::transform(Z, Z + size, P.data<double>(), Z, [step](const double& z, const double& p) { return z + step * p; });
            residual = cg.get_residual();
        }
        if (!(residual < tol))
        {
            // Opt_CG's residual precedes its last step; check the returned vector itself.
            f_LZ(Z, LP.data<double>());
            double residual_squared = 0.0;
            for (int i = 0; i < size; ++i)
            {
                const double difference = LP.data<double>()[i] - R[i];
                residual_squared += difference * difference;
            }
            Parallel_Reduce::reduce_all(residual_squared);
            residual = std::sqrt(residual_squared);
            if (!(residual < tol))
            {
                std::ostringstream message;
                message << "Z-vector CG did not converge after " << maxiter
                        << " iterations (residual=" << residual << ", tolerance=" << tol
                        << "). Excited-state forces are not valid.";
                throw std::runtime_error(message.str());
            }
        }
        if (test_force)
        {
            std::cout << "Final Z-vector:" << std::endl;
            LR_Util::print_value(Z, nstates, ld);
        }
        ModuleBase::timer::end("Z_vector", "solve_Z_CG");
    }

    inline void solve_Z_CG(std::complex<double>* const Z, std::complex<double>* R, const int& ld, const int& nstates,
        std::function<void(const std::complex<double>* const, std::complex<double>* const)> f_LZ, const bool test_force)
    {
        throw std::runtime_error("complex Z-vector solver is not implemented yet");
    }

    /// @brief Global length of one state's Z-vector: $\sum_\sigma n_k n_{occ,\sigma} n_{virt,\sigma}$.
    /// `THam` only has to expose `nk`, `nocc` and `nvirt`.
    template<typename THam>
    inline int zvec_global_dim(THam& hm, const int nspin_x)
    {
        int n_global = 0;
        for (int is = 0;is < nspin_x;++is) { n_global += hm.nk * hm.nocc[is] * hm.nvirt[is]; }
        return n_global;
    }

#ifdef __MPI
    /// @brief Gather one state's Z-vector from the `hm.pX` layout (`local`, length `ld`) into the
    /// global vector `full` (length `zvec_global_dim`), replicated on every rank. Collective.
    template<typename T, typename THam>
    inline void zvec_local_to_full(THam& hm, const int nspin_x, const T* const local, T* const full)
    {
        const int n_global = zvec_global_dim(hm, nspin_x);
        // `gather_2d_to_full` sums over the ranks, so the entries a rank does not own must be zero
        std::fill(full, full + n_global, T(0));
        int loffset = 0;
        int goffset = 0;
        for (int is = 0;is < nspin_x;++is)
        {
            const int npairs = hm.nocc[is] * hm.nvirt[is];
            const int lsize = hm.pX[is].get_local_size();
            for (int ik = 0;ik < hm.nk;++ik)
            {
                LR_Util::gather_2d_to_full(hm.pX[is], local + loffset + ik * lsize,
                    full + goffset + ik * npairs, false, hm.nvirt[is], hm.nocc[is]);
            }
            loffset += hm.nk * lsize;
            goffset += hm.nk * npairs;
        }
    }

    /// @brief Inverse of `zvec_local_to_full`: pick this rank's entries of the global vector
    /// `full` into the `hm.pX` layout `local`. No communication.
    template<typename T, typename THam>
    inline void zvec_full_to_local(THam& hm, const int nspin_x, const T* const full, T* const local)
    {
        int loffset = 0;
        int goffset = 0;
        for (int is = 0;is < nspin_x;++is)
        {
            const int npairs = hm.nocc[is] * hm.nvirt[is];
            const int lsize = hm.pX[is].get_local_size();
            for (int ik = 0;ik < hm.nk;++ik)
            {
                LR_Util::scatter_full_to_2d(hm.pX[is], full + goffset + ik * npairs,
                    local + loffset + ik * lsize, false);
            }
            loffset += hm.nk * lsize;
            goffset += hm.nk * npairs;
        }
    }
#endif

    /// @brief Solve the Z-vector equation with a dense LAPACK solve.
    ///
    /// Works for both the closed-shell (`nspin_x == 1`) and the open-shell (`nspin_x == 2`,
    /// X = [up | down]) layouts: everything is expressed through per-spin segment sizes, which
    /// collapse to the single-block case when `nspin_x == 1`.
    /// `THam` only has to expose `matrix()`, `nk`, `nocc`, `nvirt` and `pX`.
    /// @attention Every rank builds the full Hessian and solves the same system redundantly;
    /// see `solve_Z_scalapack` / `solve_Z_elpa` for the distributed solves.
    template<typename T, typename THam>
    inline void solve_Z_lapack(T* const Z, const T* const R, const int& ld, const int& nstates,
        THam& hm, const int nspin_x, const bool test_force)
    {
        ModuleBase::TITLE("Z_vector", "solve_Z_lapack");
        const int n_global = zvec_global_dim(hm, nspin_x);
        int ld_expect = 0;
        for (int is = 0;is < nspin_x;++is) { ld_expect += hm.nk * hm.pX[is].get_local_size(); }
        assert(ld == ld_expect);

        std::vector<T> hessian_full = hm.matrix();   // MO-hessian, A+B
        std::vector<T> Z_full(static_cast<std::size_t>(n_global) * nstates, T(0.0));

        // `hessian_full` is replicated, so the right-hand side must be global too.
        // `R` is distributed over `pX` (length `ld` per state), so gather it first --
        // the mirror image of the scatter of `Z_full` below. Reading `n_global` entries
        // straight out of `R` would be an out-of-bounds read as soon as
        // ld < n_global, and a plain segfault on a rank whose local size is 0.
        std::vector<T> R_full(static_cast<std::size_t>(n_global) * nstates, T(0.0));
#ifdef __MPI
        for (int istate = 0; istate < nstates; ++istate)
        {
            zvec_local_to_full(hm, nspin_x, R + istate * ld, R_full.data() + istate * n_global);
        }
#else
        std::copy(R, R + static_cast<std::size_t>(n_global) * nstates, R_full.begin());
#endif

        // use lapack to solve the linear equation
        ModuleBase::timer::start("Z_vector", "lapack_solver");
        LR_Util::lapack_linear_solver(hessian_full.data(), Z_full.data(), R_full.data(), n_global, nstates);
        ModuleBase::timer::end("Z_vector", "lapack_solver");

        // test: print full Z
        if (test_force)
        {
            std::cout << "The full Z-vector solved by LAPACK:" << std::endl;
            LR_Util::print_value(Z_full.data(), nstates, n_global);
        }

        // copy the local part of Z_full to Z
#ifdef __MPI
        for (int istate = 0; istate < nstates; ++istate)
        {
            zvec_full_to_local(hm, nspin_x, Z_full.data() + istate * n_global, Z + istate * ld);
        }
#else
        std::copy(Z_full.begin(), Z_full.end(), Z);
#endif
        if (test_force)
        {
            std::cout << "The local Z-vector solved by LAPACK:" << std::endl;
            LR_Util::print_value(Z, nstates, ld);
        }
    }

#ifdef __MPI
    /// @brief Build this rank's block of the Z-vector Hessian on the 2D block-cyclic layout `ph`
    /// (n_global x n_global), one column at a time: column g is `hm.hPsi` of the g-th unit vector,
    /// i.e. exactly the operator the CG solve applies. Only O(n_global) is replicated per rank
    /// (one column in flight), instead of the O(n_global^2) of `hm.matrix()`.
    /// Collective: every rank walks every column, since `hPsi` and the gather communicate.
    template<typename T, typename THam>
    std::vector<T> zvec_hessian_2d(THam& hm, const int nspin_x, const int ld, const Parallel_2D& ph)
    {
        ModuleBase::TITLE("Z_vector", "zvec_hessian_2d");
        ModuleBase::timer::start("Z_vector", "zvec_hessian_2d");
        const int n_global = ph.get_global_row_size();
        std::vector<T> h_loc(ph.get_local_size(), T(0));
        std::vector<T> e_full(n_global, T(0));
        std::vector<T> e_loc(ld, T(0));
        std::vector<T> he_loc(ld, T(0));
        std::vector<T> he_full(n_global, T(0));
        for (int gcol = 0; gcol < n_global; ++gcol)
        {
            e_full[gcol] = T(1);
            zvec_full_to_local(hm, nspin_x, e_full.data(), e_loc.data());
            e_full[gcol] = T(0);
            hm.hPsi(e_loc.data(), he_loc.data(), ld, 1);
            zvec_local_to_full(hm, nspin_x, he_loc.data(), he_full.data());
            const int lcol = ph.global2local_col(gcol);
            if (lcol < 0) { continue; }
            T* const h_col = h_loc.data() + static_cast<std::size_t>(lcol) * ph.get_row_size();
            for (int lrow = 0; lrow < ph.get_row_size(); ++lrow) { h_col[lrow] = he_full[ph.local2global_row(lrow)]; }
        }
        ModuleBase::timer::end("Z_vector", "zvec_hessian_2d");
        return h_loc;
    }

    /// @brief Distributed dense Z-vector solve: the Hessian (`zvec_hessian_2d`) and the
    /// right-hand side are laid out 2D block-cyclically on the BLACS grid of `hm.pX[0]`, and
    /// `linear_solver` (`scalapack_linear_solver`, `scalapack_cholesky_linear_solver` or
    /// `elpa_linear_solver`) solves them in place.
    template<typename T, typename THam>
    void solve_Z_2d(T* const Z, const T* const R, const int ld, const int nstates,
        THam& hm, const int nspin_x,
        void (*linear_solver)(T*, T*, const Parallel_2D&, const Parallel_2D&))
    {
        ModuleBase::TITLE("Z_vector", "solve_Z_2d");
        ModuleBase::timer::start("Z_vector", "solve_Z_2d");
        const int n_global = zvec_global_dim(hm, nspin_x);
        const Parallel_2D& px0 = hm.pX[0];
        // 32 suits both ScaLAPACK and ELPA; shrink it for small systems so no rank is left empty
        const int nproc_dim = std::max(px0.get_dim0(), px0.get_dim1());
        const int nb = std::max(1, std::min(32, n_global / nproc_dim));
        Parallel_2D ph;
        LR_Util::setup_2d_division(ph, nb, n_global, n_global, px0.blacs_ctxt);
        Parallel_2D pz;
        LR_Util::setup_2d_division(pz, nb, n_global, nstates, px0.blacs_ctxt);

        std::vector<T> h_loc = zvec_hessian_2d<T>(hm, nspin_x, ld, ph);

        // right-hand side: pX layout -> 2D block-cyclic on pz
        std::vector<T> z_loc(pz.get_local_size(), T(0));
        std::vector<T> r_full(n_global, T(0));
        for (int istate = 0; istate < nstates; ++istate)
        {
            zvec_local_to_full(hm, nspin_x, R + istate * ld, r_full.data());
            const int lcol = pz.global2local_col(istate);
            if (lcol < 0) { continue; }
            T* const z_col = z_loc.data() + static_cast<std::size_t>(lcol) * pz.get_row_size();
            for (int lrow = 0; lrow < pz.get_row_size(); ++lrow) { z_col[lrow] = r_full[pz.local2global_row(lrow)]; }
        }

        linear_solver(h_loc.data(), z_loc.data(), ph, pz);

        // solution: 2D block-cyclic on pz -> pX layout
        std::vector<T> z_full(static_cast<std::size_t>(n_global) * nstates, T(0));
        LR_Util::gather_2d_to_full(pz, z_loc.data(), z_full.data(), false, n_global, nstates);
        for (int istate = 0; istate < nstates; ++istate)
        {
            zvec_full_to_local(hm, nspin_x, z_full.data() + static_cast<std::size_t>(istate) * n_global, Z + istate * ld);
        }
        ModuleBase::timer::end("Z_vector", "solve_Z_2d");
    }
#endif

    /// @brief Distributed dense Z-vector solve with ScaLAPACK LU (p?gesv). Same system as
    /// `solve_Z_lapack`, but neither the Hessian nor the factorization is replicated.
    template<typename T, typename THam>
    inline void solve_Z_scalapack(T* const Z, const T* const R, const int ld, const int nstates,
        THam& hm, const int nspin_x)
    {
#ifdef __MPI
        solve_Z_2d(Z, R, ld, nstates, hm, nspin_x, &scalapack_linear_solver<T>);
#else
        throw std::runtime_error("Z-vector solver 'scalapack' needs an MPI build; use 'lapack' or 'cg'");
#endif
    }

    /// @brief Distributed dense Z-vector solve with a ScaLAPACK Cholesky factorization (p?potrf),
    /// about half the flops of the LU in `solve_Z_scalapack`. Like `solve_Z_elpa` it needs the
    /// orbital Hessian to be positive definite, but a failure is an error on every rank, not a hang.
    template<typename T, typename THam>
    inline void solve_Z_scalapack_chol(T* const Z, const T* const R, const int ld, const int nstates,
        THam& hm, const int nspin_x)
    {
#ifdef __MPI
        solve_Z_2d(Z, R, ld, nstates, hm, nspin_x, &scalapack_cholesky_linear_solver<T>);
#else
        throw std::runtime_error("Z-vector solver 'scalapack_chol' needs an MPI build; use 'lapack' or 'cg'");
#endif
    }

    /// @brief Distributed dense Z-vector solve with an ELPA Cholesky factorization: the orbital
    /// Hessian A+B is symmetric positive definite at a stable ground state, so this needs about
    /// half the flops of the LU in `solve_Z_scalapack`. Fails loudly if the Hessian is not positive
    /// definite (an unstable ground state).
    template<typename T, typename THam>
    inline void solve_Z_elpa(T* const Z, const T* const R, const int ld, const int nstates,
        THam& hm, const int nspin_x)
    {
#ifdef __MPI
        solve_Z_2d(Z, R, ld, nstates, hm, nspin_x, &elpa_linear_solver<T>);
#else
        throw std::runtime_error("Z-vector solver 'elpa' needs an MPI build; use 'lapack' or 'cg'");
#endif
    }

    /// @brief Run the configured solver (and then, for testing, every supported one).
    /// Shared by the closed- and open-shell paths; `THamL` only needs `hPsi` plus what
    /// `solve_Z_lapack` reads.
    template<typename T, typename THamL>
    inline void solve_zeq_with(T* const Z, T* const R, const int ld, const int nstates,
        THamL& ops_L, const int nspin_x, const std::string& zvec_solver, const bool test_force)
    {
        for (int i = 0; i < nstates * ld; ++i) { Z[i] = T(0.0); }   // clear Z
        if (zvec_solver == "cg")
        {
            solve_Z_CG(Z, R, ld, nstates,
                [&ops_L, ld, nstates](const T* const in, T* const out)
                { ops_L.hPsi(in, out, ld, nstates); }, test_force);
        }
        else if (zvec_solver == "lapack") { solve_Z_lapack(Z, R, ld, nstates, ops_L, nspin_x, test_force); }
        else if (zvec_solver == "scalapack") { solve_Z_scalapack(Z, R, ld, nstates, ops_L, nspin_x); }
        else if (zvec_solver == "scalapack_chol") { solve_Z_scalapack_chol(Z, R, ld, nstates, ops_L, nspin_x); }
        else if (zvec_solver == "elpa") { solve_Z_elpa(Z, R, ld, nstates, ops_L, nspin_x); }
        else { throw std::runtime_error("Unsupported Z-vector solver: " + zvec_solver); }
    }

    /// Builds the Z-vector equation's RHS from `ops_R` and solves it with `ops_L`. Extracted to
    /// a template function (rather than a generic lambda, which needs C++14) so the closed- and
    /// open-shell call sites in `Z_vector_equation` -- which pass different `Z_vector_R`/`UR` and
    /// `Z_vector_L`/`UL` types -- can share this body under the repository's C++11 baseline.
    template<typename T, typename TOpsR, typename TOpsL>
    void build_and_solve_zeq(TOpsR& ops_R, TOpsL& ops_L, const int nspin_x,
        const T* const X, container::Tensor& R, T* const Z,
        const int nloc_per_band, const int nstates, const std::string& zvec_solver, const bool test_force)
    {
        ModuleBase::timer::start("Z_vector", "Z_vector_R");
        ops_R.hPsi(X, R.template data<T>(), nloc_per_band, nstates);  // act each operator on X
        ModuleBase::timer::end("Z_vector", "Z_vector_R");
        // std::cout << "The right side of the Z-vector equation:" << std::endl;
        // LR_Util::print_value(R.template data<T>(), nstates, nloc_per_band);
        solve_zeq_with(Z, R.template data<T>(), nloc_per_band, nstates, ops_L, nspin_x, zvec_solver, test_force);
    }

    template<typename T>
    void Z_vector_equation(const GradientInputs<T>& inputs,
        const T* const X,
        T* const Z,
        const int& nstates,
        std::weak_ptr<PotHxcLR> pot,
        const std::string& spin_type,
        const std::string& in_dir,
        const bool openshell,
        const std::string& zvec_solver)
    {
        const int& nspin = inputs.nspin;
        const int& naos = inputs.nbasis;
        const std::vector<int>& nocc = inputs.nocc;
        const std::vector<int>& nvirt = inputs.nvirt;
        const UnitCell& ucell = inputs.ucell;
        const std::vector<double>& orb_cutoff = inputs.orb_cutoff;
        const Grid_Driver& gd = inputs.gd;
        const K_Vectors& kv = inputs.kv;
        const std::vector<Parallel_2D>& px = inputs.px;
        const Parallel_2D& pc = inputs.pc;
        const Parallel_Orbitals& pmat = inputs.pmat;
        const std::string& dft_functional = inputs.dft_functional;
        std::weak_ptr<PotHxcLR> pot_hxc_gs = inputs.pot_hxc_gs;
#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri = inputs.exx_lri;
        const double& exx_alpha = inputs.hybrid_alpha;
#endif
        const std::string& xc_kernel = inputs.xc_kernel;
        const psi::Psi<T>& psi_ks = inputs.psi_ks;
        const ModuleBase::matrix& eig_ks = inputs.eig_ks;
        const std::string& out_dir = inputs.out_dir;
        const std::string& ks_solver = inputs.ks_solver;
        const bool test_force = inputs.test_force;
        ModuleBase::TITLE("Z_vector", "Z_vector");
        const int nk = kv.get_nks() / nspin;
        const int nloc_per_band = openshell
            ? nk * (px[0].get_local_size() + px[1].get_local_size())
            : nk * px[0].get_local_size();
        container::Tensor R = LR_Util::newTensor<T>({ nstates, nloc_per_band });
        R.zero();

        if (openshell)
        {
            Z_vector_UR<T> ops_R(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot, pot_hxc_gs, kv, px, pc, pmat, ks_solver, dft_functional);
            Z_vector_UL<T> ops_L(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot_hxc_gs, kv, px, pc, pmat, dft_functional);
            build_and_solve_zeq(ops_R, ops_L, /*nspin_x=*/2, X, R, Z, nloc_per_band, nstates, zvec_solver, test_force);
        }
        else
        {
            Z_vector_R<T> ops_R(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot, pot_hxc_gs, kv, px, pc, pmat, in_dir, out_dir, dft_functional, spin_type);
            Z_vector_L<T> ops_L(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot_hxc_gs, kv, px, pc, pmat, spin_type, in_dir, out_dir, dft_functional);
            build_and_solve_zeq(ops_R, ops_L, /*nspin_x=*/1, X, R, Z, nloc_per_band, nstates, zvec_solver, test_force);
        }
    }
}

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_ZEQ_SOLVER_HPP
