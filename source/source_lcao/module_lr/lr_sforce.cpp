#include "source_lcao/module_lr/lr_sforce.h"

#include "source_lcao/module_lr/grad_degen.h"

#include <complex>

namespace LR
{

template <typename T>
std::vector<ModuleBase::matrix> state_forces(const GradientInputs<T>& in, const GsForceCtx<T>& gs,
                                             const StateCtx<T>& st, const int ispin, const int istate_only)
{
    const int ist_begin = (istate_only < 0) ? 0 : istate_only;
    const int nst = (istate_only < 0) ? st.nstates : 1;
    // The whole closed-shell gradient runs in the Z window (every virtual band the ground state
    // produced), not the `nvirt` window X was solved in: only there does the Z-vector enforce the
    // Brillouin condition in every occupied-virtual rotation. X is zero-padded into it, so D^X and
    // T come out bit-identical -- and Omega is untouched, so an existing finite-difference
    // reference stays valid.
    const ct::Tensor xz = pad_amplitudes<T>(st.X, st.openshell, st.para_x, st.para_x_z, st.nocc, st.nvirt, st.nk,
                                            st.nloc, st.nloc_z, ispin, ist_begin, nst);
    std::vector<double> omega(nst);
    for (int i = 0; i < nst; ++i)
    {
        omega[i] = st.ekb[ispin * st.nstates + ist_begin + i];
    }
    return force_xz(in, gs, ispin, xz, omega, ist_begin);
}

template <typename T>
std::vector<ModuleBase::matrix> state_forces_os(const GradientInputs<T>& in, const GsForceCtx<T>& gs,
                                                const StateCtx<T>& st, const int istate_only)
{
    // Open shell: there is a single eigenproblem whose vector is the concatenation
    // [up-block | down-block], and every density matrix has two independent channels. The
    // spin-orbital formulas apply verbatim -- unlike the closed-shell singlet/triplet algorithm,
    // X here is normalized over BOTH channels, so it carries no implicit sqrt(2) and none of the
    // collapsed 2/4 factors are needed.
    const int ist_begin = (istate_only < 0) ? 0 : istate_only;
    const int nst = (istate_only < 0) ? st.nstates : 1;
    // Like the closed-shell path, the whole gradient runs in the Z window: X is zero-padded into
    // it (both channels, each re-based), so D^X and T are unchanged while the Z-vector gets every
    // occupied-virtual rotation the AO basis supports.
    const ct::Tensor xz = pad_amplitudes<T>(st.X, st.openshell, st.para_x, st.para_x_z, st.nocc, st.nvirt, st.nk,
                                            st.nloc, st.nloc_z, 0, ist_begin, nst);
    std::vector<double> omega(nst);
    for (int i = 0; i < nst; ++i)
    {
        omega[i] = st.ekb[ist_begin + i];
    }
    return force_xz_os(in, gs, xz, omega, ist_begin);
}

template <typename T>
MplLvc build_mpl_state(const GradientInputs<T>& in, const GsForceCtx<T>& gs, const StateCtx<T>& st, const int ispin,
                       const std::vector<int>& group, const std::vector<ModuleBase::matrix>& diag, std::ofstream& ofs)
{
    // The multiplet's members, each widened into the Z window
    const ct::Tensor xz = pad_group<T>(st.X, st.openshell, st.para_x, st.para_x_z, st.nocc, st.nvirt, st.nk,
                                       st.nloc, st.nloc_z, ispin, group);
    const int ekb_off = st.openshell ? 0 : ispin * st.nstates;
    std::vector<double> omega(group.size());
    for (size_t k = 0; k < group.size(); ++k)
    {
        omega[k] = st.ekb[ekb_off + group[k]];
    }
    // The solver is reached only through this callback; everything else is module code.
    const MplEval eval = [&in, &gs, &st, ispin](const ct::Tensor& x, const std::vector<double>& w, const int label) {
        return st.openshell ? force_xz_os(in, gs, x, w, label) : force_xz(in, gs, ispin, x, w, label);
    };
    return build_mpl<T>(xz, omega, diag, st.nloc_z, ispin, group, in.ucell, eval, ofs);
}

template <typename T>
std::vector<MplLvc> multiplet_grad(const GradientInputs<T>& in, const GsForceCtx<T>& gs, const StateCtx<T>& st,
                                   const int ispin, const double degen_thr,
                                   const std::vector<ModuleBase::matrix>& forces, std::ofstream& ofs)
{
    const int ekb_off = st.openshell ? 0 : ispin * st.nstates;
    std::vector<double> omega(st.nstates);
    for (int ist = 0; ist < st.nstates; ++ist)
    {
        omega[ist] = st.ekb[ekb_off + ist];
    }
    std::vector<MplLvc> out;
    for (const std::vector<int>& group : group_degenerate_states(omega, degen_thr))
    {
        if (group.size() < 2)
        {
            continue;
        }
        std::vector<ModuleBase::matrix> diag;
        for (const int ist : group)
        {
            diag.push_back(forces[ist]);
        }
        out.push_back(build_mpl_state(in, gs, st, ispin, group, diag, ofs));
    }
    return out;
}

template std::vector<ModuleBase::matrix> state_forces<double>(const GradientInputs<double>&,
                                                              const GsForceCtx<double>&, const StateCtx<double>&, int,
                                                              int);
template std::vector<ModuleBase::matrix> state_forces<std::complex<double>>(
    const GradientInputs<std::complex<double>>&, const GsForceCtx<std::complex<double>>&, const StateCtx<std::complex<double>>&,
    int, int);
template std::vector<ModuleBase::matrix> state_forces_os<double>(const GradientInputs<double>&,
                                                                 const GsForceCtx<double>&, const StateCtx<double>&,
                                                                 int);
template std::vector<ModuleBase::matrix> state_forces_os<std::complex<double>>(
    const GradientInputs<std::complex<double>>&, const GsForceCtx<std::complex<double>>&,
    const StateCtx<std::complex<double>>&, int);
template MplLvc build_mpl_state<double>(const GradientInputs<double>&, const GsForceCtx<double>&,
                                        const StateCtx<double>&, int, const std::vector<int>&,
                                        const std::vector<ModuleBase::matrix>&, std::ofstream&);
template MplLvc build_mpl_state<std::complex<double>>(const GradientInputs<std::complex<double>>&,
                                                      const GsForceCtx<std::complex<double>>&,
                                                      const StateCtx<std::complex<double>>&, int,
                                                      const std::vector<int>&, const std::vector<ModuleBase::matrix>&,
                                                      std::ofstream&);
template std::vector<MplLvc> multiplet_grad<double>(const GradientInputs<double>&, const GsForceCtx<double>&,
                                                    const StateCtx<double>&, int, double,
                                                    const std::vector<ModuleBase::matrix>&, std::ofstream&);
template std::vector<MplLvc> multiplet_grad<std::complex<double>>(const GradientInputs<std::complex<double>>&,
                                                                  const GsForceCtx<std::complex<double>>&,
                                                                  const StateCtx<std::complex<double>>&, int, double,
                                                                  const std::vector<ModuleBase::matrix>&, std::ofstream&);

} // namespace LR
