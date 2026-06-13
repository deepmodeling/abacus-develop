#ifndef EXX_WAVE_REDISTRIBUTOR_H
#define EXX_WAVE_REDISTRIBUTOR_H

#include "source_basis/module_pw/pw_basis_k.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace hamilt
{

template <typename T>
class ExxWaveRedistributorCpu
{
  public:
    void setup(const ModulePW::PW_Basis_K* wfc_basis, const ModulePW::PW_Basis_K* exx_basis)
    {
        wfc_ = wfc_basis;
        exx_ = exx_basis;
        plans_.clear();
        plans_.resize(exx_->nks);
        if (exx_->poolnproc <= 1)
        {
            return;
        }
#ifdef __MPI
        for (int ik = 0; ik < exx_->nks; ++ik)
        {
            plans_[ik] = build_plan(ik);
        }
#endif
    }

    void wfc_to_exx(int ik, const T* wfc_local, T* exx_local) const
    {
        std::fill(exx_local, exx_local + exx_->npwk[ik], T(0));
        if (exx_->poolnproc <= 1)
        {
            return;
        }
#ifdef __MPI
        apply_plan(plans_[ik].forward, wfc_local, exx_local, false, T(1));
#endif
    }

    void exx_to_wfc_add(int ik, const T* exx_local, T* wfc_local, T factor) const
    {
        if (exx_->poolnproc <= 1)
        {
            return;
        }
#ifdef __MPI
        apply_plan(plans_[ik].reverse, exx_local, wfc_local, true, factor);
#endif
    }

  private:
    struct Owner
    {
        Owner() = default;
        Owner(int rank_in, int local_in) : rank(rank_in), local(local_in) {}
        int rank = -1;
        int local = -1;
    };

    struct Entry
    {
        Entry() = default;
        Entry(std::tuple<int, int, int> g_in, Owner owner_in) : g(std::move(g_in)), owner(owner_in) {}
        std::tuple<int, int, int> g;
        Owner owner;
    };

    struct TransferPlan
    {
        std::vector<int> send_counts;
        std::vector<int> send_displs;
        std::vector<int> recv_counts;
        std::vector<int> recv_displs;
        std::vector<int> send_src;
        std::vector<int> recv_dst;
    };

    struct KPlan
    {
        TransferPlan forward;
        TransferPlan reverse;
    };

#ifdef __MPI
    static MPI_Datatype mpi_type()
    {
        if (std::is_same<T, std::complex<double>>::value)
        {
            return MPI_DOUBLE_COMPLEX;
        }
        return MPI_COMPLEX;
    }
#endif

    static void prefix_counts(const std::vector<int>& counts, std::vector<int>& displs)
    {
        displs.assign(counts.size(), 0);
        for (std::size_t i = 1; i < counts.size(); ++i)
        {
            displs[i] = displs[i - 1] + counts[i - 1];
        }
    }

    std::vector<Entry> gather_entries(const ModulePW::PW_Basis_K* basis, int ik) const
    {
        std::vector<Entry> entries;
        const int npwk = basis->npwk[ik];
#ifdef __MPI
        std::vector<int> local(static_cast<std::size_t>(npwk) * 4);
        for (int ig = 0; ig < npwk; ++ig)
        {
            const auto g = basis->getgdirect(ik, ig);
            local[4 * ig] = static_cast<int>(std::lround(g.x));
            local[4 * ig + 1] = static_cast<int>(std::lround(g.y));
            local[4 * ig + 2] = static_cast<int>(std::lround(g.z));
            local[4 * ig + 3] = ig;
        }

        std::vector<int> recv_counts(basis->poolnproc, 0);
        const int local_count = static_cast<int>(local.size());
        MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, basis->pool_world);
        std::vector<int> displs;
        prefix_counts(recv_counts, displs);
        const int total = displs.back() + recv_counts.back();
        std::vector<int> gathered(total);
        MPI_Allgatherv(local.data(),
                       local_count,
                       MPI_INT,
                       gathered.data(),
                       recv_counts.data(),
                       displs.data(),
                       MPI_INT,
                       basis->pool_world);

        entries.reserve(static_cast<std::size_t>(total) / 4);
        for (int rank = 0; rank < basis->poolnproc; ++rank)
        {
            const int begin = displs[rank];
            const int end = begin + recv_counts[rank];
            for (int i = begin; i < end; i += 4)
            {
                entries.push_back({std::make_tuple(gathered[i], gathered[i + 1], gathered[i + 2]),
                                   Owner{rank, gathered[i + 3]}});
            }
        }
#else
        entries.reserve(npwk);
        for (int ig = 0; ig < npwk; ++ig)
        {
            const auto g = basis->getgdirect(ik, ig);
            entries.push_back({std::make_tuple(static_cast<int>(std::lround(g.x)),
                                               static_cast<int>(std::lround(g.y)),
                                               static_cast<int>(std::lround(g.z))),
                               Owner{0, ig}});
        }
#endif
        return entries;
    }

    TransferPlan build_transfer(const std::map<std::tuple<int, int, int>, Owner>& source_owner,
                                const std::vector<Entry>& target_entries) const
    {
        TransferPlan plan;
        const int nproc = exx_->poolnproc;
        const int rank = exx_->poolrank;
        std::vector<std::vector<int>> send_by_rank(nproc);
        std::vector<std::vector<int>> recv_by_rank(nproc);

        for (const Entry& target: target_entries)
        {
            const auto it = source_owner.find(target.g);
            if (it == source_owner.end())
            {
                continue;
            }
            const Owner src = it->second;
            const Owner dst = target.owner;
            if (src.rank == rank)
            {
                send_by_rank[dst.rank].push_back(src.local);
            }
            if (dst.rank == rank)
            {
                recv_by_rank[src.rank].push_back(dst.local);
            }
        }

        plan.send_counts.assign(nproc, 0);
        plan.recv_counts.assign(nproc, 0);
        for (int ip = 0; ip < nproc; ++ip)
        {
            plan.send_counts[ip] = static_cast<int>(send_by_rank[ip].size());
            plan.recv_counts[ip] = static_cast<int>(recv_by_rank[ip].size());
        }
        prefix_counts(plan.send_counts, plan.send_displs);
        prefix_counts(plan.recv_counts, plan.recv_displs);

        for (int ip = 0; ip < nproc; ++ip)
        {
            plan.send_src.insert(plan.send_src.end(), send_by_rank[ip].begin(), send_by_rank[ip].end());
            plan.recv_dst.insert(plan.recv_dst.end(), recv_by_rank[ip].begin(), recv_by_rank[ip].end());
        }
        return plan;
    }

    KPlan build_plan(int ik) const
    {
        const auto wfc_entries = gather_entries(wfc_, ik);
        const auto exx_entries = gather_entries(exx_, ik);
        std::map<std::tuple<int, int, int>, Owner> wfc_owner;
        std::map<std::tuple<int, int, int>, Owner> exx_owner;
        for (const Entry& entry: wfc_entries)
        {
            wfc_owner.emplace(entry.g, entry.owner);
        }
        for (const Entry& entry: exx_entries)
        {
            exx_owner.emplace(entry.g, entry.owner);
        }

        KPlan plan;
        plan.forward = build_transfer(wfc_owner, exx_entries);
        plan.reverse = build_transfer(exx_owner, wfc_entries);
        return plan;
    }

    void apply_plan(const TransferPlan& plan, const T* source, T* target, bool add, T factor) const
    {
#ifdef __MPI
        std::vector<T> send_buffer(plan.send_src.size());
        for (std::size_t i = 0; i < plan.send_src.size(); ++i)
        {
            send_buffer[i] = source[plan.send_src[i]];
        }
        std::vector<T> recv_buffer(plan.recv_dst.size());
        MPI_Alltoallv(send_buffer.data(),
                      plan.send_counts.data(),
                      plan.send_displs.data(),
                      mpi_type(),
                      recv_buffer.data(),
                      plan.recv_counts.data(),
                      plan.recv_displs.data(),
                      mpi_type(),
                      exx_->pool_world);
        if (add)
        {
            for (std::size_t i = 0; i < recv_buffer.size(); ++i)
            {
                target[plan.recv_dst[i]] += factor * recv_buffer[i];
            }
        }
        else
        {
            for (std::size_t i = 0; i < recv_buffer.size(); ++i)
            {
                target[plan.recv_dst[i]] = recv_buffer[i];
            }
        }
#endif
    }

    const ModulePW::PW_Basis_K* wfc_ = nullptr;
    const ModulePW::PW_Basis_K* exx_ = nullptr;
    std::vector<KPlan> plans_;
};

} // namespace hamilt

#endif // EXX_WAVE_REDISTRIBUTOR_H
