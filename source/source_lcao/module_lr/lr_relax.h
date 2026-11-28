#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_RELAX_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_RELAX_H

#include "source_lcao/module_lr/state_track.h"
#include "source_lcao/module_lr/utils/lr_util.h"

#include <ostream>
#include <vector>

class UnitCell;
class TwoCenterIntegrator;
struct Input_para;

namespace LR
{
/// Which excited state (or degenerate multiplet) a geometry relaxation follows.
template <typename T>
class RelaxTarget
{
  public:
    /// Resolve and validate lr_target_state / lr_target_spin; call once the dimensions are final.
    void setup(const Input_para& inp, bool openshell, int nspin, int nstates, std::ostream& log);

    bool active() const;
    int ispin() const;
    int state() const;
    const std::vector<int>& group() const;
    /// Index of the followed state inside `ekb` for the current spin layout.
    int ekb_off(int nstates) const;

    /// Match the old subspace to current energy groups; update `state()` and the references.
    void follow(const Input_para& inp, const UnitCell& ucell, const std::vector<double>& orb_cutoff,
                const TwoCenterIntegrator& overlap, const psi::Psi<T>& psi_ks, const Parallel_Blkcyclic& para_c,
                const Parallel_Orbitals& pmat, const std::vector<Parallel_Blkcyclic>& para_x,
                const std::vector<int>& nocc, const std::vector<int>& nvirt, const T* amplitudes, int nloc,
                int nstates, const std::vector<double>& energies, std::ostream& log);

    /// Fill the multiplet group from the current excitation energies (`ekb` starts at this spin).
    void resolve_group(const Input_para& inp, const double* ekb, int nstates, std::ostream& log);
    /// Excitation energy being minimised: the target's own, or the multiplet average.
    double omega(const double* ekb, int nstates) const;
    /// Store the JT mixed root as the next reference.
    void save_mixed(const T* amplitudes, int nloc, const std::vector<double>& mixing);

  private:
    bool active_ = false;
    bool openshell_ = false;
    int ispin_ = 0;
    int state_ = -1;
    std::vector<int> group_;
    std::vector<T> previous_;
    RootBasis<T> basis_;
};

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_RELAX_H
