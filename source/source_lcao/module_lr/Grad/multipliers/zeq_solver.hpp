#pragma once
#include "zeq_solver.h"
#include "source_base/opt_cg.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"

namespace LR
{
    // inline void solve_Z_CG(double* const Z, const double* const R, const int& ld, const int& nstates,
    //     std::function<void(const double* const, double* const)> f_LZ)
    // Opt_CG's interfaces have no const qualifier for the pointer
    inline void solve_Z_CG(double* const Z, double* R, const int& ld, const int& nstates,
        std::function<void(const double* const, double* const)> f_LZ)
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
        int final_iter = 0;
        std::cout << "Start solving Z-vector equaiton with CG method ..." << std::endl;
        for (int iter = 0; iter < maxiter; ++iter)
        {
            if (residual < tol)
            {
                final_iter = iter;
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
        std::cout << "Final Z-vector:" << std::endl;
        LR_Util::print_value(Z, nstates, ld);
        ModuleBase::timer::end("Z_vector", "solve_Z_CG");
    }

    inline void solve_Z_CG(std::complex<double>* const Z, std::complex<double>* R, const int& ld, const int& nstates,
        std::function<void(const std::complex<double>* const, std::complex<double>* const)> f_LZ)
    {
        throw std::runtime_error("complex Z-vector solver is not implemented yet");
    }

    /// @brief Solve the Z-vector equation with a dense LAPACK solve.
    ///
    /// Works for both the closed-shell (`nspin_x == 1`) and the open-shell (`nspin_x == 2`,
    /// X = [up | down]) layouts: everything is expressed through per-spin segment sizes, which
    /// collapse to the single-block case when `nspin_x == 1`.
    /// `THam` only has to expose `matrix()`, `nk`, `nocc`, `nvirt` and `pX`.
    template<typename T, typename THam>
    inline void solve_Z_lapack(T* const Z, const T* const R, const int& ld, const int& nstates,
        const THam& hm, const int nspin_x = 1)
    {
        ModuleBase::TITLE("Z_vector", "solve_Z_lapack");
        std::vector<int> npairs(nspin_x), ldim_is(nspin_x), gdim_is(nspin_x);
        int n_global = 0, ld_expect = 0;
        for (int is = 0;is < nspin_x;++is)
        {
            npairs[is] = hm.nocc[is] * hm.nvirt[is];
            gdim_is[is] = hm.nk * npairs[is];
            ldim_is[is] = hm.nk * hm.pX[is].get_local_size();
            n_global += gdim_is[is];
            ld_expect += ldim_is[is];
        }
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
            int loffset = istate * ld;
            int goffset = istate * n_global;
            for (int is = 0;is < nspin_x;++is)
            {
                for (int ik = 0;ik < hm.nk;++ik)
                {
                    LR_Util::gather_2d_to_full(hm.pX[is],
                        R + loffset + ik * hm.pX[is].get_local_size(),
                        R_full.data() + goffset + ik * npairs[is],
                        false, hm.nvirt[is], hm.nocc[is]);
                }
                loffset += ldim_is[is];
                goffset += gdim_is[is];
            }
        }
#else
        std::copy(R, R + static_cast<std::size_t>(n_global) * nstates, R_full.begin());
#endif

        // use lapack to solve the linear equation
        ModuleBase::timer::start("Z_vector", "lapack_solver");
        LR_Util::lapack_linear_solver(hessian_full.data(), Z_full.data(), R_full.data(), n_global, nstates);
        ModuleBase::timer::end("Z_vector", "lapack_solver");

        // test: print full Z
        std::cout << "The full Z-vector solved by LAPACK:" << std::endl;
        LR_Util::print_value(Z_full.data(), nstates, n_global);

        // copy the local part of Z_full to Z
        for (int istate = 0; istate < nstates; ++istate)
        {
            int loffset = istate * ld;
            int goffset = istate * n_global;
            for (int is = 0;is < nspin_x;++is)
            {
                for (int ik = 0;ik < hm.nk;++ik)
                {
                    LR_Util::scatter_full_to_2d(hm.pX[is],
                        Z_full.data() + goffset + ik * npairs[is],
                        Z + loffset + ik * hm.pX[is].get_local_size(), false);
                }
                loffset += ldim_is[is];
                goffset += gdim_is[is];
            }
        }
        std::cout << "The local Z-vector solved by LAPACK:" << std::endl;
        LR_Util::print_value(Z, nstates, ld);
    }

    /// @brief Run the configured solver (and then, for testing, every supported one).
    /// Shared by the closed- and open-shell paths; `THamL` only needs `hPsi` plus what
    /// `solve_Z_lapack` reads.
    template<typename T, typename THamL>
    inline void solve_zeq_with(T* const Z, T* const R, const int ld, const int nstates,
        const THamL& ops_L, const int nspin_x, const std::string& zvec_solver)
    {
        for (int i = 0; i < nstates * ld; ++i) { Z[i] = T(0.0); }   // clear Z
        if (zvec_solver == "cg")
        {
            solve_Z_CG(Z, R, ld, nstates,
                [&ops_L, ld, nstates](const T* const in, T* const out)
                { ops_L.hPsi(in, out, ld, nstates); });
        }
        else if (zvec_solver == "lapack") { solve_Z_lapack(Z, R, ld, nstates, ops_L, nspin_x); }
        else { throw std::runtime_error("Unsupported Z-vector solver: " + zvec_solver); }
    }

    template<typename T>
    void Z_vector_equation(const T* const X,
        T* const Z,
        const std::string& xc_kernel,
        const int& nstates,
        const int& nspin,
        const int& naos,
        const std::vector<int>& nocc,
        const std::vector<int>& nvirt,
        const UnitCell& ucell,
        const std::vector<double>& orb_cutoff,
        const Grid_Driver& gd,
        const psi::Psi<T>& psi_ks,
        const ModuleBase::matrix& eig_ks,
#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri,
        const double& exx_alpha,
#endif 
        std::weak_ptr<PotHxcLR> pot,
        std::weak_ptr<PotHxcLR> pot_hxc_gs,
        const K_Vectors& kv,
        const std::vector<Parallel_2D>& px,
        const Parallel_2D& pc,
        const Parallel_Orbitals& pmat,
        const std::string& spin_type,
        const bool openshell = false,
        const std::string& zvec_solver = "cg")
    {
        ModuleBase::TITLE("Z_vector", "Z_vector");
        const int nk = kv.get_nks() / nspin;
        const int nloc_per_band = openshell
            ? nk * (px[0].get_local_size() + px[1].get_local_size())
            : nk * px[0].get_local_size();
        container::Tensor R = LR_Util::newTensor<T>({ nstates, nloc_per_band });
        R.zero();

        auto build_and_solve = [&](auto& ops_R, auto& ops_L, const int nspin_x)
            {
                ModuleBase::timer::start("Z_vector", "Z_vector_R");
                ops_R.hPsi(X, R.template data<T>(), nloc_per_band, nstates);  // act each operator on X
                ModuleBase::timer::end("Z_vector", "Z_vector_R");
                std::cout << "The right side of the Z-vector equation:" << std::endl;
                LR_Util::print_value(R.template data<T>(), nstates, nloc_per_band);
                solve_zeq_with(Z, R.template data<T>(), nloc_per_band, nstates, ops_L, nspin_x, zvec_solver);
            };

        if (openshell)
        {
            Z_vector_UR<T> ops_R(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot, pot_hxc_gs, kv, px, pc, pmat);
            Z_vector_UL<T> ops_L(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot_hxc_gs, kv, px, pc, pmat);
            build_and_solve(ops_R, ops_L, /*nspin_x=*/2);
        }
        else
        {
            Z_vector_R<T> ops_R(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot, pot_hxc_gs, kv, px, pc, pmat, spin_type);
            Z_vector_L<T> ops_L(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot_hxc_gs, kv, px, pc, pmat, spin_type);
            build_and_solve(ops_R, ops_L, /*nspin_x=*/1);
        }
    }
}
