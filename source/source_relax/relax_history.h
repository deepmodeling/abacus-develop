#ifndef RELAX_HISTORY_H
#define RELAX_HISTORY_H

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

/**
 * @file relax_history.h
 * @brief Shared formatting of per-step convergence history for the running log.
 *
 * Both relaxation paths (Relax in relax_sync.cpp and IonCellOptimizer in
 * relax_nsync.cpp) record the largest force / stress of each step and print a
 * summary when relaxation finishes. This free function keeps that summary
 * layout identical across the two paths.
 */

/**
 * @brief Format a per-step history for the running log.
 *
 * Layout rules:
 * - Up to max_full entries are printed in full, per_line values per line.
 * - Beyond max_full, only the first and last keep entries are shown, with an
 *   omission marker in between, so very long runs stay readable.
 */
inline std::string format_relax_history(const std::vector<double>& hist,
                                        const int per_line = 5,
                                        const int max_full = 100,
                                        const int keep = 10)
{
    std::ostringstream out;
    out << std::scientific << std::setprecision(3);

    const int n = static_cast<int>(hist.size());
    const bool truncate = (n > max_full);
    const int head = truncate ? keep : n;

    int printed = 0;
    auto print_one = [&out, &printed, per_line](const double value) {
        if (printed % per_line == 0)
        {
            out << "\n  ";
        }
        else
        {
            out << " ";
        }
        out << value;
        ++printed;
    };

    for (int i = 0; i < head; ++i)
    {
        print_one(hist[i]);
    }
    if (truncate)
    {
        const int omitted = n - 2 * keep;
        out << "\n   ... (omitted " << omitted << " step(s)) ...";
        for (int i = n - keep; i < n; ++i)
        {
            print_one(hist[i]);
        }
    }
    out << "\n";
    return out.str();
}

#endif // RELAX_HISTORY_H
