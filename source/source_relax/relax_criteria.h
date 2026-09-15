#ifndef RELAX_CRITERIA_H
#define RELAX_CRITERIA_H

#include <string>

/**
 * @brief INPUT-derived settings that the relaxation algorithms need.
 *
 * These values used to be read straight out of the global PARAM inside the
 * algorithms. That made the algorithms impossible to unit test without
 * mutating global state, which in turn is why their tests had to switch off
 * access control with `#define private public`.
 *
 * They are now filled once by the relaxation driver and passed down
 * explicitly. Leaf functions still take only the individual values they use;
 * this struct exists to keep the plumbing signatures readable.
 */
struct Relax_Criteria
{
    // The defaults below deliberately mirror the corresponding Input_para
    // defaults, so that a caller (or a test) which leaves a field alone gets
    // exactly the behaviour it got when these values were read from PARAM.
    double force_thr = -1;        ///< Force convergence threshold, Ry/Bohr
    double force_thr_ev = -1;     ///< The same threshold in eV/Angstrom, reconciled by ReadInput
    double stress_thr = 0.5;      ///< Stress convergence threshold, kbar
    bool fixed_ibrav = false;     ///< Keep the Bravais lattice type fixed while relaxing the cell
    std::string out_level = "ie"; ///< Output verbosity; "ie" prints per-step energy to stdout
    int test_relax_method = 0;    ///< Debug verbosity for the relaxation algorithms
};

#endif
