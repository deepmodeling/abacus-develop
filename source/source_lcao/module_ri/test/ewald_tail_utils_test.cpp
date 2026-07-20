#include "../ewald_tail_utils.h"

#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

namespace
{
using Cell = std::array<int, 3>;
using EwaldVqDetail::TailDecision;
using EwaldVqDetail::TailKey;
using EwaldVqDetail::TailMode;
using EwaldVqDetail::TailStats;
using EwaldVqDetail::TailTolerances;

void require(const bool condition, const std::string& message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void test_expand_odd_period()
{
    require(EwaldVqDetail::expand_odd_period(Cell{13, 13, 3}, {true, true, false}, 1)
                == Cell{15, 15, 3},
            "period expansion must preserve a nonperiodic direction");
    require(EwaldVqDetail::expand_odd_period(Cell{1, 3, 5}, {true, true, true}, 2)
                == Cell{5, 7, 9},
            "period expansion must add two cells per guard layer");
}

void test_shell_cells()
{
    const auto shell = EwaldVqDetail::make_shell_cells(Cell{3, 3, 1}, Cell{5, 5, 1});
    require(shell.size() == 16, "5x5 shell around 3x3 must contain sixteen cells");
    require(std::set<Cell>(shell.begin(), shell.end()).size() == shell.size(),
            "shell cells must not contain duplicates");
    for (const Cell& cell: shell)
    {
        require(cell[2] == 0, "nonperiodic shell coordinate must remain zero");
        require(!EwaldVqDetail::cell_in_period(cell, Cell{3, 3, 1}),
                "shell must exclude production cells");
        require(EwaldVqDetail::cell_in_period(cell, Cell{5, 5, 1}),
                "shell must stay inside the probe period");
    }
}

void test_hermitian_partner_and_owner()
{
    const TailKey key{0, 1, {2, -1, 0}};
    const TailKey partner = EwaldVqDetail::hermitian_partner(key);
    require(partner == TailKey{1, 0, {-2, 1, 0}},
            "Hermitian partner must swap atoms and negate the cell");
    require(EwaldVqDetail::shell_owner(key, 7) == EwaldVqDetail::shell_owner(partner, 7),
            "Hermitian partners must use the same deterministic owner");
}

void test_tail_decisions()
{
    const TailTolerances tolerances{1e-12, 1e-10, 1e-8};
    TailStats stats;
    stats.reference_norm = 10.0;
    stats.max_norm = 5e-10;
    stats.sum_norm = 2e-7;
    stats.coverage_complete = true;

    require(EwaldVqDetail::decide_tail(stats, tolerances, TailMode::Warn, 0, 3)
                == TailDecision::Warn,
            "warn mode must not enlarge a failed shell");
    require(EwaldVqDetail::decide_tail(stats, tolerances, TailMode::Enlarge, 0, 3)
                == TailDecision::Expand,
            "enlarge mode must expand a failed shell");
    require(EwaldVqDetail::decide_tail(stats, tolerances, TailMode::Strict, 3, 3)
                == TailDecision::Fail,
            "strict mode must fail at the expansion limit");
    require(EwaldVqDetail::decide_tail(stats, tolerances, TailMode::Enlarge, 3, 3)
                == TailDecision::Fail,
            "enlarge mode must fail at the expansion limit");

    stats.max_norm = 5e-11;
    stats.sum_norm = 5e-8;
    stats.hermitian_residual = 2e-9;
    require(!EwaldVqDetail::tail_passes(stats, tolerances),
            "a failed Hermitian invariant must reject an otherwise converged shell");
    require(EwaldVqDetail::decide_tail(stats, tolerances, TailMode::Warn, 0, 3)
                == TailDecision::Fail,
            "warn mode must still fail a Hermitian invariant violation");
    stats.hermitian_residual = 0.0;
    require(EwaldVqDetail::decide_tail(stats, tolerances, TailMode::Enlarge, 0, 3)
                == TailDecision::Accept,
            "maximum, summed, and Hermitian tail criteria must accept");

    stats.coverage_complete = false;
    require(EwaldVqDetail::decide_tail(stats, tolerances, TailMode::Warn, 0, 3)
                == TailDecision::Fail,
            "warn mode must still fail incomplete physical coverage");
    require(EwaldVqDetail::decide_tail(stats, tolerances, TailMode::Strict, 0, 3)
                == TailDecision::Fail,
            "strict mode must reject incomplete physical coverage");
}

void test_adaptive_exact_rmesh_times()
{
    require(EwaldVqDetail::uses_adaptive_realspace_range(TailMode::Enlarge),
            "enlarge mode must use the converged real-space range");
    require(EwaldVqDetail::uses_adaptive_realspace_range(TailMode::Strict),
            "strict mode must use the converged real-space range");
    require(!EwaldVqDetail::uses_adaptive_realspace_range(TailMode::Warn),
            "warn mode must preserve the requested production range");

    require(std::abs(EwaldVqDetail::adaptive_exact_rmesh_times(0.01, {6.0}, {6.0}) - 1.0)
                < 1e-14,
            "adaptive evaluation must keep exact pair tables through the nonoverlap boundary");
    require(std::abs(EwaldVqDetail::adaptive_exact_rmesh_times(0.5, {4.0, 7.0}, {8.0, 6.0})
                     - 5.0 / 3.0)
                < 1e-14,
            "heterogeneous cutoffs must cover every ordered ABFS pair");
    require(std::abs(EwaldVqDetail::adaptive_exact_rmesh_times(2.0, {6.0}, {6.0}) - 2.0)
                < 1e-14,
            "an already larger requested radial mesh must be preserved");
}

void test_missing_block_classification()
{
    using EwaldVqDetail::TailBlockCoverage;
    require(EwaldVqDetail::classify_tail_block(true, true, false) == TailBlockCoverage::Common,
            "two evaluated blocks must be subtracted directly");
    require(EwaldVqDetail::classify_tail_block(true, false, true)
                == TailBlockCoverage::AnalyticCancellation,
            "a nonoverlapping far-field pair may use analytic multipole cancellation");
    require(EwaldVqDetail::classify_tail_block(false, true, true)
                == TailBlockCoverage::AnalyticCancellation,
            "far-field analytic cancellation is symmetric between bare and Gaussian maps");
    require(EwaldVqDetail::classify_tail_block(true, false, false) == TailBlockCoverage::Incomplete,
            "a near-field missing Gaussian block must not be treated as zero");
    require(EwaldVqDetail::classify_tail_block(false, true, false) == TailBlockCoverage::Incomplete,
            "a near-field missing bare block must not be treated as zero");

    const auto missing_gaussian = EwaldVqDetail::plan_tail_block(true, false, 8.0, 10.0, 7.0);
    require(missing_gaussian.coverage == TailBlockCoverage::Common
                && !missing_gaussian.use_bare_multipole
                && missing_gaussian.use_gaussian_multipole,
            "the mixed-support interval must combine exact bare and analytic Gaussian blocks");

    const auto missing_bare = EwaldVqDetail::plan_tail_block(false, true, 8.0, 7.0, 10.0);
    require(missing_bare.coverage == TailBlockCoverage::Common
                && missing_bare.use_bare_multipole
                && !missing_bare.use_gaussian_multipole,
            "the mixed-support interval must combine analytic bare and exact Gaussian blocks");

    const auto both_far = EwaldVqDetail::plan_tail_block(false, false, 12.0, 10.0, 11.0);
    require(both_far.coverage == TailBlockCoverage::AnalyticCancellation
                && !both_far.use_bare_multipole
                && !both_far.use_gaussian_multipole,
            "two absent far-field blocks must cancel without constructing a shared multipole tensor");

    require(EwaldVqDetail::plan_tail_block(true, false, 6.0, 10.0, 7.0).coverage
                == TailBlockCoverage::Incomplete,
            "a missing Gaussian block inside its compact support must remain an error");
}

void test_tail_sample_aggregation()
{
    TailStats stats;
    EwaldVqDetail::accumulate_tail_sample(
        stats,
        EwaldVqDetail::TailSample{TailKey{0, 1, {1, 0, 0}}, 5.0, 4.9, 0.1, 8.0});
    EwaldVqDetail::accumulate_tail_sample(
        stats,
        EwaldVqDetail::TailSample{TailKey{1, 1, {2, 0, 0}}, 3.0, 2.7, 0.3, 12.0});
    require(stats.shell_blocks == 2, "tail aggregation must count every sampled block");
    require(std::abs(stats.sum_norm - 0.4) < 1e-14, "tail aggregation must sum block norms");
    require(std::abs(stats.max_norm - 0.3) < 1e-14, "tail aggregation must retain the maximum norm");
    require(stats.worst_key == TailKey{1, 1, {2, 0, 0}}, "tail aggregation must retain the worst key");
    require(std::abs(stats.worst_bare_norm - 3.0) < 1e-14,
            "tail aggregation must retain the worst bare norm");
    require(std::abs(stats.worst_gaussian_norm - 2.7) < 1e-14,
            "tail aggregation must retain the worst Gaussian norm");
    require(std::abs(stats.worst_distance - 12.0) < 1e-14,
            "tail aggregation must retain the worst-key distance");
}

void test_explicit_key_grouping_and_far_field_boundary()
{
    const std::vector<TailKey> keys{
        TailKey{1, 0, {2, 0, 0}},
        TailKey{0, 1, {-2, 0, 0}},
        TailKey{1, 1, {3, 0, 0}},
    };
    const auto grouped = EwaldVqDetail::group_tail_keys_by_atom(keys);
    require(grouped.size() == 2, "explicit tail keys must be grouped by center atom");
    require(grouped.at(0).size() == 1 && grouped.at(1).size() == 2,
            "grouping must preserve every explicit J,R key exactly once");
    require(EwaldVqDetail::far_multipoles_equivalent(12.0, 10.0, 11.0),
            "nonoverlap of both compact supports permits far-field cancellation");
    require(!EwaldVqDetail::far_multipoles_equivalent(10.5, 10.0, 11.0),
            "overlap of either compact support requires explicit evaluation");
}
} // namespace

int main()
{
    try
    {
        test_expand_odd_period();
        test_shell_cells();
        test_hermitian_partner_and_owner();
        test_tail_decisions();
        test_adaptive_exact_rmesh_times();
        test_missing_block_classification();
        test_tail_sample_aggregation();
        test_explicit_key_grouping_and_far_field_boundary();
    }
    catch (const std::exception& error)
    {
        std::cerr << "Ewald tail utility test failed: " << error.what() << '\n';
        return 1;
    }
    std::cout << "Ewald tail utility tests passed\n";
    return 0;
}
