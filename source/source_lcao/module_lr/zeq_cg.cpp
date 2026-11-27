#include "zeq_solver.h"
#include "source_base/opt_cg.h"
#include <algorithm>
#include <cmath>
#include <sstream>
namespace LR
{
    void solve_Z_CG(double* const Z, double* R, const int& ld, const int& nstates,
        std::function<void(const double* const, double* const)> f_LZ, const bool test_force)
    {
        ModuleBase::TITLE("Z_vector", "solve_Z_CG");
        ModuleBase::timer::start("Z_vector", "solve_Z_CG");
        const int maxiter = 100;
        double tol = 1e-6;
        double residual_squared = 0.0;
        for (int i = 0; i < nstates * ld; ++i)
        {
            residual_squared += R[i] * R[i];
        }
        Parallel_Reduce::reduce_all(residual_squared);
        double residual = std::sqrt(residual_squared);
        const int size = nstates * ld;
        container::Tensor P = LR_Util::newTensor<double>({ nstates, ld });   // step length
        container::Tensor LP = LR_Util::newTensor<double>({ nstates, ld }); //f_LZ(P)
        std::fill_n(Z, size, 0.0);

        ModuleBase::Opt_CG cg;
        cg.allocate(size);
        cg.init_b(R);
        std::cout << "Start solving Z-vector equaiton with CG method ..." << std::endl;
        int iterations = 0;
        for (; iterations < maxiter; ++iterations)
        {
            if (residual < tol)
            {
                break;
            }
            cg.next_direct(LP.data<double>(), 0, P.data<double>());
            residual = cg.get_residual();
            std::cout << "iter=" << iterations << " residual=" << residual << std::endl;
            // std::cout << "Z=" << std::endl;
            // LR_Util::print_value(P.data<double>(), nstates, ld);
            // std::cout << "LZ before=" << std::endl;
            // LR_Util::print_value(LP.data<double>(), nstates, ld);
            f_LZ(P.data<double>(), LP.data<double>());  // L: act each operators on P
            // std::cout << "LZ=" << std::endl;
            // LR_Util::print_value(LP.data<double>(), nstates, ld);
            int ifPD = 0;   //???
            double step = cg.step_length(LP.data<double>(), P.data<double>(), ifPD);
            const double* direction = P.data<double>();
            for (int i = 0; i < size; ++i)
            {
                Z[i] += step * direction[i];
            }
            residual = cg.get_residual();
        }
        // Recursive convergence is insufficient: certify the vector actually returned.
        f_LZ(Z, LP.data<double>());
        residual_squared = 0.0;
        for (int i = 0; i < size; ++i)
        {
            const double difference = LP.data<double>()[i] - R[i];
            residual_squared += difference * difference;
        }
        Parallel_Reduce::reduce_all(residual_squared);
        residual = std::sqrt(residual_squared);
        std::cout << "Z-vector CG true residual=" << residual << " tolerance=" << tol << std::endl;
        if (!(residual < tol))
        {
            std::ostringstream message;
            message << "Z-vector CG did not converge after " << iterations
                    << " iterations (true residual=" << residual << ", tolerance=" << tol
                    << "). Excited-state forces are not valid.";
            ModuleBase::timer::end("Z_vector", "solve_Z_CG");
            throw std::runtime_error(message.str());
        }
        if (test_force)
        {
            std::cout << "Final Z-vector:" << std::endl;
            LR_Util::print_value(Z, nstates, ld);
        }
        ModuleBase::timer::end("Z_vector", "solve_Z_CG");
    }

    void solve_Z_CG(std::complex<double>* const Z, std::complex<double>* R, const int& ld, const int& nstates,
        std::function<void(const std::complex<double>* const, std::complex<double>* const)> f_LZ, const bool test_force)
    {
        throw std::runtime_error("complex Z-vector solver is not implemented yet");
    }

}
