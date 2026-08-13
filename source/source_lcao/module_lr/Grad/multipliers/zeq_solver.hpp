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

    template<typename T>
    inline void solve_Z_lapack(T* const Z, const T* const R, const int& ld, const int& nstates,
        const HamiltLR<T>& hm)
    {
        ModuleBase::TITLE("Z_vector", "solve_Z_lapack");
        assert(ld == hm.nk * hm.pX[0].get_local_size());

        std::vector<T> hessian_full = hm.matrix();   // MO-hessian, A+B

        const int n_global = hm.nk * hm.nocc[0] * hm.nvirt[0];
        std::vector<T> Z_full = std::vector<T>(n_global * nstates, T(0.0));

        // use lapack to solve the linear equation
        ModuleBase::timer::start("Z_vector", "lapack_solver");
        LR_Util::lapack_linear_solver(hessian_full.data(), Z_full.data(), R, n_global, nstates);
        ModuleBase::timer::end("Z_vector", "lapack_solver");

        // test: print full Z
        std::cout << "The full Z-vector solved by LAPACK:" << std::endl;
        LR_Util::print_value(Z_full.data(), nstates, n_global);

        // copy the local part of Z_full to Z
        for (int istate = 0; istate < nstates; ++istate)
        {
            const int global_offset = istate * n_global;
            const int offset = istate * ld;
            LR_Util::scatter_full_to_2d(hm.pX[0], Z_full.data() + global_offset, Z + offset, false);
        }
        std::cout << "The local Z-vector solved by LAPACK:" << std::endl;
        LR_Util::print_value(Z, nstates, ld);
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
        const std::string& zvec_solver = "cg")
    {
        ModuleBase::TITLE("Z_vector", "Z_vector");
        const int nk = kv.get_nks() / nspin;
        // 1. the right-hand side of Z-vector equation
        const int nloc_per_band = nk * px[0].get_local_size();
        container::Tensor R = LR_Util::newTensor<T>({ nstates, nloc_per_band });
        R.zero();
        Z_vector_R<T> ops_R(xc_kernel, nspin, naos, nocc, nvirt,
            ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
            exx_lri, exx_alpha,
#endif 
            pot, pot_hxc_gs, kv, px, pc, pmat, spin_type);
        ModuleBase::timer::start("Z_vector", "Z_vector_R");
        ops_R.hPsi(X, R.data<T>(), nloc_per_band, nstates);  // act each operators on X
        ModuleBase::timer::end("Z_vector", "Z_vector_R");
        std::cout << "The right side of the Z-vector equation:" << std::endl;
        LR_Util::print_value(R.data<T>(), nstates, nloc_per_band);

        // 2. the left-hand side of Z-vector equation
        // Z-vector (need a init?)
        Z_vector_L<T> ops_L(xc_kernel, nspin, naos, nocc, nvirt,
            ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
            exx_lri, exx_alpha,
#endif 
            pot_hxc_gs, kv, px, pc, pmat, spin_type);

        // 3. solve Z-vector equation
        auto solve = [&](const std::string& solver)
            {
                // clear Z
                for (int i = 0; i < nstates * nloc_per_band; ++i) { Z[i] = T(0.0); }
                if (solver == "cg")
                {
                    solve_Z_CG(Z, R.data<T>(), nloc_per_band, nstates,
                        std::bind(&HamiltLR<T>::hPsi, &ops_L, std::placeholders::_1, std::placeholders::_2, nloc_per_band, nstates));
                }
                else if (solver == "lapack")
                {
                    solve_Z_lapack(Z, R.data<T>(), nloc_per_band, nstates, ops_L);
                }
                else
                {
                    throw std::runtime_error("Unsupported Z-vector solver: " + solver);
                }
            };

        solve(zvec_solver);

        // test: try supported solvers one by one
        const std::vector<std::string> supported_solvers = { "cg", "lapack" };
        for (const auto& s : supported_solvers)
            solve(s);

        // test: set Z to 0
        // for (int i = 0; i < nstates * nloc_per_band; ++i) { Z[i] = T(0.0); }
    }
}