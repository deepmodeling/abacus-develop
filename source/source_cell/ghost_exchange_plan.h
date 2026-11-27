#ifndef GHOST_EXCHANGE_PLAN_H
#define GHOST_EXCHANGE_PLAN_H

#include <cstddef>
#include <limits>
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

    bool is_valid(const std::size_t natom) const
    {
        const std::size_t nswap = sendnum.size();
        const bool communication_sizes_match = recvnum.size() == nswap && firstrecv.size() == nswap
                                                && send_atom_indices.size() == nswap
                                                && sendproc.size() == nswap && recvproc.size() == nswap;
        const bool nswap_fits_int = nswap <= static_cast<std::size_t>(std::numeric_limits<int>::max());
        const bool mapping_matches_atoms = nprocs != 1 || mapping.size() == natom;
        return communication_sizes_match && nswap_fits_int && mapping_matches_atoms;
    }
};

#endif // GHOST_EXCHANGE_PLAN_H
