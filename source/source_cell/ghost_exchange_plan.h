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
    int nprocs;

    GhostExchangePlan() : nprocs(1) {}
};

#endif // GHOST_EXCHANGE_PLAN_H
