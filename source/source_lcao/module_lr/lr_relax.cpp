#include "source_lcao/module_lr/lr_relax.h"

#include "source_base/tool_quit.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_lr/grad_degen.h"
#include "source_lcao/module_lr/lr_amp.h"

#include <algorithm>
#include <complex>
#include <type_traits>

namespace LR
{

template <typename T>
void RelaxTarget<T>::setup(const Input_para& inp, bool openshell, int nspin, int nstates, std::ostream& log)
{
    active_ = (inp.calculation == "relax");
    openshell_ = openshell;
    if (!active_)
    {
        return;
    }

    // The Z-vector equation has no complex solver (see zeq_solver.h), so an
    // excited-state gradient only exists at gamma. Fail here rather than after the SCF.
    if (!std::is_same<T, double>::value)
    {
        ModuleBase::WARNING_QUIT(
            "ESolver_LR", "excited-state relaxation currently requires gamma-only sampling: the complex "
                          "Z-vector solver is not implemented.");
    }
    if (inp.lr_target_state >= nstates)
    {
        ModuleBase::WARNING_QUIT(
            "ESolver_LR", "lr_target_state is beyond the states actually solved (lr_nstates <= 0 expands to "
                          "all particle-hole pairs, which may be fewer than requested).");
    }

    // `openshell` is only settled after the ground-state occupations have been read, which is why
    // this cannot live in the input-file checks
    const std::string& spin = inp.lr_target_spin;
    if (openshell_)
    {
        // An open-shell calculation solves one spin-conserving channel, so there is nothing to
        // choose: whatever lr_target_spin says, this is the state that gets relaxed. Only an
        // explicit `triplet` is worth mentioning -- `singlet` is the default and expresses no
        // intent, and `updown` is already the right name for this channel.
        if (spin == "triplet")
        {
            log << " WARNING: lr_target_spin=triplet is ignored. This is an"
                   " open-shell calculation with a single spin-conserving channel (updown), which is"
                   " what the relaxation will follow."
                << std::endl;
        }
        ispin_ = 0;
    }
    else
    {
        // Closed shell is the opposite case: singlet and triplet are genuinely different states
        // with different gradients, so `updown` here is ambiguous rather than redundant -- it
        // usually means lr_unrestricted was meant to be set.
        if (spin == "updown")
        {
            ModuleBase::WARNING_QUIT(
                "ESolver_LR", "lr_target_spin=updown, but this is a closed-shell calculation, where singlet and "
                              "triplet are separate states with separate gradients. Pick one of them, or set "
                              "lr_unrestricted to run spin-unrestricted.");
        }
        ispin_ = (spin == "triplet") ? 1 : 0;
        if (ispin_ >= nspin)
        {
            ModuleBase::WARNING_QUIT(
                "ESolver_LR", "lr_target_spin=triplet requires nspin=2: the triplet channel is not built here.");
        }
    }
    state_ = inp.lr_target_state; // seed; overlap takes over from step 2
    log << " Excited-state relaxation follows state " << inp.lr_target_state << " of the "
        << (openshell_ ? "updown" : (ispin_ == 1 ? "triplet" : "singlet"))
        << " channel, tracked by cross-geometry orbital and amplitude overlap." << std::endl;
}

template <typename T>
bool RelaxTarget<T>::active() const
{
    return active_;
}

template <typename T>
int RelaxTarget<T>::ispin() const
{
    return ispin_;
}

template <typename T>
int RelaxTarget<T>::state() const
{
    return state_;
}

template <typename T>
const std::vector<int>& RelaxTarget<T>::group() const
{
    return group_;
}

template <typename T>
int RelaxTarget<T>::ekb_off(int nstates) const
{
    return openshell_ ? 0 : ispin_ * nstates;
}

template <typename T>
void RelaxTarget<T>::follow(const Input_para& inp, const UnitCell& ucell, const std::vector<double>& orb_cutoff,
                            const TwoCenterIntegrator& overlap, const psi::Psi<T>& psi_ks,
                            const Parallel_Blkcyclic& para_c, const Parallel_Orbitals& pmat,
                            const std::vector<Parallel_Blkcyclic>& para_x, const std::vector<int>& nocc,
                            const std::vector<int>& nvirt, const T* amplitudes, int nloc, int nstates,
                            const std::vector<double>& energies, std::ostream& log)
{
    if (!active_)
    {
        return;
    }
    // Compare roots in a common electron-hole basis using exact cross-geometry AO overlap.
    const bool reference_first = LR_Util::tolower(inp.lr_degen_mode) == "jt";
    const RootInputs<T> inputs{ucell, orb_cutoff, overlap, psi_ks, para_c, pmat, para_x, nocc, nvirt, openshell_,
                               ispin_, reference_first};
    const bool single_state = LR_Util::tolower(inp.lr_degen_mode) == "state";
    const double group_threshold = single_state ? 0.0 : inp.lr_degen_thr;
    follow_cross_root(inputs, amplitudes, nloc, nstates, inp.lr_target_state, energies, group_threshold, state_,
                      previous_, basis_, log);
}

template <typename T>
void RelaxTarget<T>::resolve_group(const Input_para& inp, const double* ekb, int nstates, std::ostream& log)
{
    // Keyed on `state_`, not on `lr_target_state`: the overlap tracking re-chooses the
    // followed root at every ionic step, and the multiplet has to be the one containing the root
    // actually being followed. They differ as soon as two surfaces have crossed.
    group_.clear();
    if (LR_Util::tolower(inp.lr_degen_mode) == "state")
    {
        return;
    }
    const int off = ekb_off(nstates);
    std::vector<double> omegas(nstates);
    for (int ist = 0; ist < nstates; ++ist)
    {
        omegas[ist] = ekb[off + ist];
    }
    const std::vector<std::vector<int>> groups = group_degenerate_states(omegas, inp.lr_degen_thr);
    for (const std::vector<int>& g : groups)
    {
        if (std::find(g.begin(), g.end(), state_) == g.end())
        {
            continue;
        }
        // a one-member group means the target is not degenerate here, and averaging over it would
        // be the single-state path with extra steps
        const bool average = LR_Util::tolower(inp.lr_degen_mode) == "average";
        if (average && g.size() > 3)
        {
            log << " WARNING: tracked group has " << g.size()
                << " roots; average supports at most three. Following the selected single root." << std::endl;
        }
        else if (!average && g.size() > 3)
        {
            ModuleBase::WARNING_QUIT(
                "ESolver_LR",
                "JT forces support only two- or three-state multiplets; reduce lr_degen_thr or use state mode");
        }
        else if (g.size() > 1)
        {
            group_ = g;
        }
        break;
    }
}

template <typename T>
double RelaxTarget<T>::omega(const double* ekb, int nstates) const
{
    if (group_.empty())
    {
        return ekb[ekb_off(nstates) + state_];
    }
    const int off = ekb_off(nstates);
    double sum = 0.0;
    for (const int ist : group_)
    {
        sum += ekb[off + ist];
    }
    return sum / static_cast<double>(group_.size());
}

template <typename T>
void RelaxTarget<T>::save_mixed(const T* amplitudes, int nloc, const std::vector<double>& mixing)
{
    save_mixed_root(amplitudes, nloc, group_, mixing, previous_);
}

template class RelaxTarget<double>;
template class RelaxTarget<std::complex<double>>;

} // namespace LR
