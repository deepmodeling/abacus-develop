#ifndef GHOST_EXCHANGE_PLAN_H
#define GHOST_EXCHANGE_PLAN_H

#include <vector>

struct GhostExchangePlanView
{
    int nswap;
    int* sendnum;
    int* recvnum;
    int* firstrecv;
    int* sendproc;
    int* recvproc;
    int* mapping;
    std::vector<int*> sendlist_storage;
};

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

    GhostExchangePlanView make_view();

    bool is_valid(const int natom) const
    {
        const int nswap = static_cast<int>(sendnum.size());
        const bool communication_sizes_match = static_cast<int>(recvnum.size()) == nswap &&
                                               static_cast<int>(firstrecv.size()) == nswap &&
                                               static_cast<int>(send_atom_indices.size()) == nswap &&
                                               static_cast<int>(sendproc.size()) == nswap &&
                                               static_cast<int>(recvproc.size()) == nswap;
        const bool mapping_matches_atoms = nprocs != 1 || static_cast<int>(mapping.size()) == natom;
        return communication_sizes_match && mapping_matches_atoms;
    }
};

#endif // GHOST_EXCHANGE_PLAN_H
