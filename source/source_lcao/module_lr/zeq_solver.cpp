#include "zeq_solver.h"
#include "hamilt_zeq_l.h"
#include "hamilt_zeq_r.h"
namespace LR
{
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
            // Concatenated spin-channel dimensions are prepared before operator construction.
            const int global_dimension = nk * (nocc[0] * nvirt[0] + nocc[1] * nvirt[1]);
            Z_vector_UR<T> ops_R(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot, pot_hxc_gs, kv, px, pc, pmat, ks_solver, dft_functional, nk, nloc_per_band, global_dimension);
            Z_vector_UL<T> ops_L(xc_kernel, nspin, naos, nocc, nvirt,
                ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot_hxc_gs, kv, px, pc, pmat, dft_functional, nk, nloc_per_band, global_dimension);
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
template void Z_vector_equation<double>(const GradientInputs<double>&, const double* const, double* const, const int&, std::weak_ptr<PotHxcLR>, const std::string&, const std::string&, const bool, const std::string&);
template void Z_vector_equation<std::complex<double>>(const GradientInputs<std::complex<double>>&, const std::complex<double>* const, std::complex<double>* const, const int&, std::weak_ptr<PotHxcLR>, const std::string&, const std::string&, const bool, const std::string&);
}
