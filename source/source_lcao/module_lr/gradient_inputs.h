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
    const Parallel_2D& pc;
    const std::vector<Parallel_2D>& px;
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
#ifdef __EXX
    const std::shared_ptr<Exx_LRI<T>>& exx_lri;
    double hybrid_alpha;
#endif
};
template <typename T>
std::vector<ModuleBase::matrix> evaluate_closed_shell_force(
    const GradientInputs<T>& inputs, LR_Force<T>& force_terms,
    const module_dm::DensityMatrix<T, double>& dm_gs,
    const ct::Tensor& Xz, const ct::Tensor& Z,
    const std::vector<double>& omega, int label_begin, int ispin);
template <typename T>
std::vector<ModuleBase::matrix> evaluate_open_shell_force(
    const GradientInputs<T>& inputs, LR_Force<T>& force_terms,
    const module_dm::DensityMatrix<T, double>& dm_gs,
    const ct::Tensor& Xz, const ct::Tensor& Z,
    const std::vector<double>& omega, int label_begin);
}
#endif
