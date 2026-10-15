#ifndef GHOST_EXCHANGE_PLAN_H
#define GHOST_EXCHANGE_PLAN_H

#include <vector>

struct GhostExchangePlan
{
    std::vector<int> sendnum;
    std::vector<int> recvnum;
    std::vector<int> firstrecv;
    std::vector<std::vector<int> > send_atom_indices;
    std::vector<int> sendproc;
    std::vector<int> recvproc;
    // For a single rank, map every local atom (owned and ghost) to its owned
    // atom index. Distributed DeePMD obtains this information from exchange
    // metadata instead.
    std::vector<int> mapping;
    int nprocs;

    GhostExchangePlan() : nprocs(1) {}
};

#endif // GHOST_EXCHANGE_PLAN_H
