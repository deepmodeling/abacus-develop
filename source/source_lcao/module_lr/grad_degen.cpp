#include "grad_degen.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

namespace LR
{
    std::vector<std::vector<int>> group_degenerate_states(const std::vector<double>& omega,
        const double thr)
    {
        const int nst = static_cast<int>(omega.size());
        std::vector<std::vector<int>> groups;
        if (nst == 0) { return groups; }

        // Casida output happens to come out ascending, but nothing here should depend on that.
        std::vector<int> order(nst);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(),
            [&omega](const int a, const int b) { return omega[a] < omega[b]; });

        if (thr <= 0.0)
        {
            for (int i = 0; i < nst; ++i) { groups.push_back(std::vector<int>(1, order[i])); }
            return groups;
        }

        double anchor = omega[order[0]];
        groups.push_back(std::vector<int>(1, order[0]));
        for (int i = 1; i < nst; ++i)
        {
            const int ist = order[i];
            // measured against the group's first (lowest) member, so the spread inside a group is
            // bounded by `thr` however many members it collects
            if (omega[ist] - anchor < thr)
            {
                groups.back().push_back(ist);
            }
            else
            {
                groups.push_back(std::vector<int>(1, ist));
                anchor = omega[ist];
            }
        }
        for (std::vector<int>& g : groups) { std::sort(g.begin(), g.end()); }
        return groups;
    }

    std::vector<std::pair<int, int>> degenerate_pairs(const int d)
    {
        std::vector<std::pair<int, int>> pairs;
        if (d < 2) { return pairs; }
        pairs.reserve(static_cast<size_t>(d) * (d - 1) / 2);
        for (int k = 0; k < d; ++k)
        {
            for (int l = k + 1; l < d; ++l) { pairs.push_back(std::make_pair(k, l)); }
        }
        return pairs;
    }

}
