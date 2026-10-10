#include "state_track.h"
#include "state_ovlp.h"
#include "state_group.h"
#include "grad_degen.h"
#include "source_psi/psi.h"
#include "lr_amp.h"
#include "source_base/parallel_2d.h"
#include "source_cell/unitcell.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include <algorithm>
#include <utility>

namespace LR
{
namespace
{
template <typename T>
RootBasis<T> capture_basis(const RootInputs<T>& inputs)
{
    RootBasis<T> basis;
    basis.lattice = inputs.cell.latvec;
    basis.lat0 = inputs.cell.lat0;
    basis.nbasis = inputs.pc.get_global_row_size();
    basis.nbands = inputs.nocc[0] + inputs.nvirt[0];
    basis.nocc = inputs.nocc;
    basis.nvirt = inputs.nvirt;
    basis.layout = {inputs.pc.get_block_size(), inputs.pc.get_dim0(), inputs.pc.get_dim1(),
        inputs.pc.get_coord_row(), inputs.pc.get_coord_col(), inputs.pc.get_row_size(), inputs.pc.get_col_size()};
    for (const auto& px : inputs.px)
    {
        basis.layout.push_back(px.get_block_size());
        basis.layout.push_back(px.get_row_size());
        basis.layout.push_back(px.get_col_size());
    }

    for (int atom = 0; atom < inputs.cell.nat; ++atom)
    {
        const auto position = inputs.cell.get_tau(atom);
        basis.positions.push_back(position);
    }
    // Snapshot only this rank's column-major coefficient block; no AO-by-band gather.
    const int coefficient_count = inputs.pc.get_local_size();
    const T zero = T(0);
    for (int spin = 0; spin < inputs.orbitals.get_nk(); ++spin)
    {
        std::vector<T> coefficients(coefficient_count, zero);
        for (int band = 0; band < inputs.pc.get_col_size(); ++band)
        {
            for (int ao = 0; ao < inputs.pc.get_row_size(); ++ao)
            {
                coefficients[band * inputs.pc.get_row_size() + ao] = inputs.orbitals(spin, band, ao);
            }
        }
        basis.coefficients.push_back(std::move(coefficients));
    }
    return basis;
}


template <typename T>
void project_local_references(const RootInputs<T>& inputs, const RootBasis<T>& old,
    const RootBasis<T>& current, const std::vector<double>& ao, std::vector<T>& previous,
    std::vector<std::vector<T>>& group)
{
    Parallel_2D po;
#ifdef __MPI
    po.set(current.nbands, current.nbands, inputs.pc.get_block_size(), inputs.pc.blacs_ctxt);
#else
    po.set_serial(current.nbands, current.nbands);
#endif
    const std::vector<T> typed_ao(ao.begin(), ao.end());
    std::vector<std::vector<T>*> references{&previous};
    for (auto& root : group) { references.push_back(&root); }
    int offset = 0;
    const std::vector<int> spins = inputs.openshell ? std::vector<int>{0, 1} : std::vector<int>{inputs.channel};
    for (const int spin : spins)
    {
        const Parallel_2D& px = inputs.px[spin];
        const int no = inputs.nocc[spin];
        const int nv = inputs.nvirt[spin];
        const int local_count = px.get_local_size();
        // Build each spin's MO overlap once, then reuse it for the single and group references.
        const auto local_mo = mo_overlap_dist(old.coefficients[spin], typed_ao,
            current.coefficients[spin], inputs.pc, inputs.pmat, po);
        for (auto* reference : references)
        {
            std::vector<T> local_reference(local_count);
            std::copy_n(reference->begin() + offset, local_count, local_reference.begin());
            const auto projected = project_reference_dist(local_reference, local_mo, po, px, no, nv);
            std::copy(projected.begin(), projected.end(), reference->begin() + offset);
        }
        offset += px.get_local_size();
    }
}

template <typename T>
int follow_subspace(const RootInputs<T>& inputs, const T* amplitudes,
    const std::vector<std::vector<int>>& groups, std::vector<T>& previous,
    const std::vector<std::vector<T>>& old_group, std::ostream& log)
{
    const int nlocal = previous.size();
    int nstates = 0;
    for (const auto& group : groups) { nstates += group.size(); }
    const int old_dimension = old_group.size();
    std::vector<T> overlaps(old_dimension * nstates, T(0));
    std::vector<T> reference(nstates, T(0));
    for (int root = 0; root < nstates; ++root)
    {
        for (int i = 0; i < nlocal; ++i)
        {
            reference[root] += LR_Util::get_conj(previous[i]) * amplitudes[root * nlocal + i];
            for (int a = 0; a < old_dimension; ++a)
            {
                overlaps[a * nstates + root] += LR_Util::get_conj(old_group[a][i])
                    * amplitudes[root * nlocal + i];
            }
        }
    }
    const int overlap_count = overlaps.size();
    Parallel_Reduce::reduce_all(overlaps.data(), overlap_count);
    Parallel_Reduce::reduce_all(reference.data(), nstates);
    const std::vector<std::complex<double>> complex_overlap(overlaps.begin(), overlaps.end());
    std::vector<double> reference_size(nstates);
    for (int root = 0; root < nstates; ++root) { reference_size[root] = std::abs(reference[root]); }
    const auto match = match_state_group(complex_overlap, old_dimension, nstates, groups, reference_size, inputs.reference_first);
    for (size_t ig = 0; ig < groups.size(); ++ig)
    {
        log << " EXCITED-STATE RELAX: candidate roots";
        for (const int root : groups[ig]) { log << " " << root; }
        log << "; reference projection weight " << match.reference_weights[ig] << std::endl;
    }
    const auto& selected = groups[match.group];
    int target = selected.front();
    for (const int root : selected)
    {
        if (reference_size[root] > reference_size[target]) { target = root; }
    }
    log << " EXCITED-STATE RELAX: subspace dimensions " << old_dimension << " -> " << selected.size()
        << "; old/new coverage " << match.old_coverage << " / " << match.new_coverage
        << "; reference-first " << inputs.reference_first << "; subspace score " << match.subspace_score
        << "; best/runner-up score " << match.score << " / " << match.runner_up << "; singular values";
    for (const double value : match.singular_values) { log << " " << value; }
    log << "; selected roots";
    for (const int root : selected) { log << " " << root; }
    const auto largest_reference = std::max_element(reference_size.begin(), reference_size.end());
    const int reference_root = std::distance(reference_size.begin(), largest_reference);
    log << "; target " << target << "; reference overlap " << reference_size[target]
        << "; best reference root " << reference_root << " (" << *largest_reference << ")" << std::endl;
    if (!inputs.reference_first && std::find(selected.begin(), selected.end(), reference_root) == selected.end())
    {
        log << " WARNING: the group match excludes the root with largest single/JT reference overlap."
            << std::endl;
    }
    if (groups.size() > 1 && match.score - match.runner_up <= 1e-10)
    {
        log << " WARNING: ambiguous subspace match; the secondary overlap criterion resolves the tie."
            << std::endl;
    }
    // Losing other old-group directions is expected after JT splitting: diagnose the
    // retained actual branch's weight, rather than treating d -> 1 as branch loss.
    const bool low_subspace = match.old_coverage < 0.5 || match.new_coverage < 0.5;
    const bool low_tracking = inputs.reference_first ? match.score < 0.5 : low_subspace;
    if (low_tracking)
    {
        log << " WARNING: low tracking coverage; continuing with the best candidate."
            << " The surface may be discontinuous; consider a smaller step or more roots." << std::endl;
    }
    if (nlocal > 0) { previous.assign(amplitudes + target * nlocal, amplitudes + (target + 1) * nlocal); }
    else { previous.clear(); }
    return target;
}
}

template <typename T>
void follow_cross_root(const RootInputs<T>& inputs, const T* amplitudes,
    const int nlocal, const int nstates, const int initial_state, const std::vector<double>& energies,
    const double group_threshold, int& target,
    std::vector<T>& previous, RootBasis<T>& basis, std::ostream& log)
{
    if (energies.size() != static_cast<size_t>(nstates))
    {
        throw std::invalid_argument("LR root overlap: energy/root counts differ");
    }
    RootBasis<T> current = capture_basis(inputs);
    if (!basis.positions.empty())
    {
        const auto difference = current.lattice - basis.lattice;
        const double cell_change = std::abs(difference.e11) + std::abs(difference.e12) + std::abs(difference.e13)
            + std::abs(difference.e21) + std::abs(difference.e22) + std::abs(difference.e23)
            + std::abs(difference.e31) + std::abs(difference.e32) + std::abs(difference.e33);
        if (cell_change > 1e-12 || std::abs(current.lat0 - basis.lat0) > 1e-12
            || current.nbasis != basis.nbasis || current.nbands != basis.nbands
            || current.nocc != basis.nocc || current.nvirt != basis.nvirt
            || current.layout != basis.layout
            || current.positions.size() != basis.positions.size()
            || current.coefficients.size() != basis.coefficients.size())
        {
            throw std::invalid_argument("LR root overlap: cell, orbital window or MPI layout changed between steps");
        }
        const auto ao = cross_ao_overlap(inputs.cell, basis.positions, inputs.cutoff, inputs.integrator, inputs.pmat);
        project_local_references(inputs, basis, current, ao, previous, basis.group_amplitudes);
        log << " EXCITED-STATE RELAX: root overlap uses old-new AO and MO bases." << std::endl;
    }
    const auto groups = group_degenerate_states(energies, group_threshold);
    if (!basis.group_amplitudes.empty())
    {
        target = follow_subspace(inputs, amplitudes, groups, previous, basis.group_amplitudes, log);
    }
    else { follow_root(amplitudes, nlocal, nstates, initial_state, target, previous, log); }
    for (const auto& group : groups)
    {
        if (std::find(group.begin(), group.end(), target) == group.end()) { continue; }
        for (const int root : group)
        {
            std::vector<T> local(nlocal);
            if (nlocal > 0) { std::copy_n(amplitudes + root * nlocal, nlocal, local.begin()); }
            current.group_amplitudes.push_back(std::move(local));
        }
        break;
    }
    basis = std::move(current);
}

template void follow_cross_root<double>(const RootInputs<double>&, const double*, int, int, int,
    const std::vector<double>&, double, int&, std::vector<double>&, RootBasis<double>&, std::ostream&);
template void follow_cross_root<std::complex<double>>(const RootInputs<std::complex<double>>&,
    const std::complex<double>*, int, int, int, const std::vector<double>&, double, int&, std::vector<std::complex<double>>&,
    RootBasis<std::complex<double>>&, std::ostream&);
}
