#ifndef ABACUS_LR_GRADIENT_INPUTS_H
#define ABACUS_LR_GRADIENT_INPUTS_H
#include "lr_force.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
namespace LR
{
// Explicit, borrowed inputs for one gradient evaluation. No solver ownership or caches.
template <typename T>
struct GradientInputs
{
    const UnitCell& ucell;
    const K_Vectors& kv;
    const Grid_Driver& gd;
    const std::vector<double>& orb_cutoff;
    const Parallel_Orbitals& pmat;
    const Parallel_Blkcyclic& pc;
    const std::vector<Parallel_Blkcyclic>& px;
    const psi::Psi<T>& psi_ks;
    const ModuleBase::matrix& eig_ks;
    const std::vector<int>& nocc;
    const std::vector<int>& nvirt;
    int nspin;
    int nk;
    int nbasis;
    int nloc;
    const std::string& xc_kernel;
    const std::string& dft_functional;
    const std::string& ks_solver;
    bool test_force;
    bool excited_relax;
    const std::string& out_dir;
    int my_rank;
    const std::vector<std::string>& spin_types;
    const std::vector<std::shared_ptr<PotHxcLR>>& pot;
    const std::shared_ptr<PotHxcLR>& pot_hxc_gs;
    std::ofstream& ofs;
    double response_tolerance;   ///< Absolute global combined CG residual norm in Ry.
    int response_max_iterations; ///< Maximum CG response updates.
#ifdef __EXX
    const std::shared_ptr<Exx_LRI<T>>& exx_lri;
    double hybrid_alpha;
#endif
};
template <typename T>
std::vector<ModuleBase::matrix> evaluate_closed_shell_force(
    const GradientInputs<T> &inputs,
    LR_Force<T> &force_terms,
    const module_dm::DensityMatrix<T, double> &dm_gs,
    const ct::Tensor &Xz,
    const ct::Tensor &Z,
    const std::vector<double> &omega,
    int label_begin,
    int ispin,
    bool emit_state_forces);
/**
 * @brief Contract ordinary or polarized unrestricted LR force contributions.
 * @param[in] inputs Borrowed geometry, full spin windows, kernels and controls.
 * @param[in,out] force_terms Geometry-local analytic force contractions.
 * @param[in] dm_gs Spin-resolved ground-state density.
 * @param[in] Xz Packed up/down amplitudes; one row per requested quadratic state.
 * @param[in] Z Matching orbital-response rows.
 * @param[in] omega Energies in Ry, one per row; may be common pair energies.
 * @param[in] label_begin Zero-based state label for ordinary force output.
 * @param[in] emit_state_forces True for physical state output; false for NAC polarization.
 * @return Owned atom-major force matrices in Ry/bohr.
 */
template <typename T>
std::vector<ModuleBase::matrix> evaluate_open_shell_force(
    const GradientInputs<T> &inputs,
    LR_Force<T> &force_terms,
    const module_dm::DensityMatrix<T, double> &dm_gs,
    const ct::Tensor &Xz,
    const ct::Tensor &Z,
    const std::vector<double> &omega,
    int label_begin,
    bool emit_state_forces);
}
#endif
