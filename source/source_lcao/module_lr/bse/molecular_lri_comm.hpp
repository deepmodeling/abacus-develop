//=======================
// AUTHOR : Ziqing Guan
// DATE :   2026-03-22
//=======================
#include "molecular_lri.h"
#include <algorithm>
#include <cstddef> // offsetof
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MPI
#include <mpi.h>
#endif
namespace BSE
{

struct BlockHead
{
    int gr0;   // global row start
    int gc0;   // global col start
    int nr;    // rows in block
    int nc;    // cols in block
};

#ifdef __MPI
inline MPI_Datatype mpi_type_blockhead()
{
    static MPI_Datatype dt = MPI_DATATYPE_NULL;
    static bool committed = false;
    if (!committed)
    {
        int blen[4] = {1, 1, 1, 1};
        MPI_Aint disp[4];
        MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};

        disp[0] = static_cast<MPI_Aint>(offsetof(BlockHead, gr0));
        disp[1] = static_cast<MPI_Aint>(offsetof(BlockHead, gc0));
        disp[2] = static_cast<MPI_Aint>(offsetof(BlockHead, nr));
        disp[3] = static_cast<MPI_Aint>(offsetof(BlockHead, nc));

        MPI_Type_create_struct(4, blen, disp, types, &dt);
        MPI_Type_commit(&dt);
        committed = true;
    }
    return dt;
}
#endif

/// @brief W[k_ai][k_bj] to 2d local matrix WA[aik1, bjk2]
template <typename T>
void MolecularLRI<T>::transform_k_2dlocal(std::vector<T>& m_2d,
    const std::map<int, std::map<int, RI::Tensor<T>>>& m_lri,
    const Parallel_2D& pm_2d, const double factor)
{
    ModuleBase::TITLE("MolecularLRI", "transform_k_2dlocal");
    ModuleBase::timer::start("MolecularLRI", "transform_k_2dlocal");
    const int npair = this->nocc * this->nvirt;
    const double fac_base = factor * 2.0; // factor 2 for Ha → Ry; k-point weight multiplies below
    const int nb = pm_2d.get_block_size();
#ifdef __MPI
    MPI_Datatype mpitype_blockhead = mpi_type_blockhead();
    // 0. outer loop of k1: communicate per 64 k1
    constexpr int comm_nk1 = 64;
    std::vector<int> send_head_counts(GlobalV::NPROC, 0), recv_head_counts(GlobalV::NPROC, 0);
    std::vector<int> send_buffer_counts(GlobalV::NPROC, 0), recv_buffer_counts(GlobalV::NPROC, 0);
    std::vector<std::size_t> send_buffer_counts_c(GlobalV::NPROC, 0);
    std::vector<int> shdispls(GlobalV::NPROC, 0), rhdispls(GlobalV::NPROC, 0); //displacements of block heads
    std::vector<int> sbdispls(GlobalV::NPROC, 0), rbdispls(GlobalV::NPROC, 0); //displacements of buffer
    std::vector<int> cursor_head(GlobalV::NPROC, 0), cursor_buffer(GlobalV::NPROC, 0);
    std::vector<BlockHead> send_heads, recv_heads;
    std::vector<T> send_buffers, recv_buffers;
    int max_send_head_total = 0, max_recv_head_total = 0;
    int max_send_buffer_total = 0, max_recv_buffer_total = 0;
    for (int k1_start = 0; k1_start < this->nk; k1_start += comm_nk1)
    {
        std::fill(send_head_counts.begin(), send_head_counts.end(), 0);
        std::fill(send_buffer_counts.begin(), send_buffer_counts.end(), 0);
        std::fill(send_buffer_counts_c.begin(), send_buffer_counts_c.end(), 0);
        std::fill(recv_head_counts.begin(), recv_head_counts.end(), 0);
        std::fill(recv_buffer_counts.begin(), recv_buffer_counts.end(), 0);
        const int k1_end = std::min(k1_start+comm_nk1, this->nk);
        // 1. calculate and coummunicate block counts, then calculate block displacements
        for (int kai = k1_start; kai < k1_end; ++kai)
        {
            if ( !this->is_local_k1[kai] ) continue;
            const int row_base = kai * npair;
            for (const int kbj : this->LR_lri.k2_indices)
            {
                const int col_base = kbj * npair;
                for (int j = 0; j < npair; )
                {
                    const int global_col = col_base + j;
                    const int j_next = std::min(global_col/nb*nb + nb - col_base, npair);
                    for (int i = 0; i < npair; )
                    {
                        const int global_row = row_base + i;
                        const int i_next = std::min(global_row/nb*nb + nb - row_base, npair);
                        if (!pm_2d.in_this_processor(global_row, global_col))
                        {
                            const int owner = pm_2d.owner_processor(global_row, global_col);
                            ++send_head_counts[owner];
                            send_buffer_counts_c[owner] += std::size_t(i_next - i) * std::size_t(j_next - j);
                        }
                        i = i_next;
                    }
                    j = j_next;
                }
            }
        }
        for (int p = 0; p < GlobalV::NPROC; ++p)
        {
            if (send_buffer_counts_c[p] > std::numeric_limits<int>::max())
            {
                throw std::overflow_error("in transform_k_2dlocal: overflow converting to int!");
            }
            send_buffer_counts[p] = static_cast<int>(send_buffer_counts_c[p]);
        }
        assert(send_head_counts.at(GlobalV::MY_RANK) == 0);
        assert(send_buffer_counts.at(GlobalV::MY_RANK) == 0);
        MPI_Alltoall(send_head_counts.data(), 1, MPI_INT, recv_head_counts.data(), 1, MPI_INT, pm_2d.comm());
        MPI_Alltoall(send_buffer_counts.data(), 1, MPI_INT,
                     recv_buffer_counts.data(), 1, MPI_INT, pm_2d.comm());
        
        int send_head_total = 0, recv_head_total = 0, send_buffer_total = 0, recv_buffer_total = 0;
        for (int p = 0; p < GlobalV::NPROC; ++p)
        {
            shdispls[p] = send_head_total; send_head_total += send_head_counts[p];
            rhdispls[p] = recv_head_total; recv_head_total += recv_head_counts[p];
            sbdispls[p] = send_buffer_total; send_buffer_total += send_buffer_counts[p];
            rbdispls[p] = recv_buffer_total; recv_buffer_total += recv_buffer_counts[p];
        }
        if (send_head_total > max_send_head_total)
            { max_send_head_total = send_head_total; send_heads.reserve(max_send_head_total); }
        if (recv_head_total > max_recv_head_total)
            { max_recv_head_total = recv_head_total; recv_heads.reserve(max_recv_head_total); }
        if (send_buffer_total > max_send_buffer_total)
            { max_send_buffer_total = send_buffer_total; send_buffers.reserve(max_send_buffer_total); }
        if (recv_buffer_total > max_recv_buffer_total)
            { max_recv_buffer_total = recv_buffer_total; recv_buffers.reserve(max_recv_buffer_total); }

        // 2. prepare block heads and buffers for send and recv
        send_heads.resize(send_head_total);
        recv_heads.resize(recv_head_total);
        send_buffers.resize(send_buffer_total);
        recv_buffers.resize(recv_buffer_total);
        cursor_head = shdispls;
        cursor_buffer = sbdispls;
        for (int kai = k1_start; kai < k1_end; ++kai)
        {
            if (!this->is_local_k1[kai] ) continue;
            const int row_base = kai * npair;
            const double fac = fac_base * this->kv.wk[kai]; // k-point weight, normalized as sum = 1
            for (const int kbj : this->LR_lri.k2_indices)
            {
                const int col_base = kbj * npair;
                const RI::Tensor<T>& m_kai_kbj = m_lri.at(kai).at(kbj);
                for (int j = 0; j < npair; )
                {
                    const int global_col = col_base + j;
                    const int j_next = std::min(global_col/nb*nb + nb - col_base, npair);
                    for (int i = 0; i < npair; )
                    {
                        const int global_row = row_base + i;
                        const int i_next = std::min(global_row/nb*nb + nb - row_base, npair);
                        const int owner = pm_2d.owner_processor(global_row, global_col);
                        if (owner == GlobalV::MY_RANK)
                        {
                            const int lr0 = pm_2d.global2local_row(global_row);
                            const int lc0 = pm_2d.global2local_col(global_col);
                            const int lld = pm_2d.get_row_size();
                            for(int jj = j; jj < j_next; ++jj)
                            {
                                for(int ii = i; ii < i_next; ++ii)
                                {
                                    const int lr = lr0 + (ii - i);
                                    const int lc = lc0 + (jj - j);
                                    const std::size_t idx_2d = lr + std::size_t(lc) * std::size_t(lld);
                                    m_2d[idx_2d] += (*m_kai_kbj.data)[ii + jj * npair] * fac;
                                }
                            }
                        }
                        else
                        {
                            BlockHead& head = send_heads[cursor_head[owner]++];
                            head.gr0 = global_row;
                            head.gc0 = global_col;
                            head.nr = i_next - i;
                            head.nc = j_next - j;
                            for(int jj = j; jj < j_next; ++jj)
                            {
                                for (int ii = i; ii < i_next; ++ii)
                                {
                                    send_buffers[cursor_buffer[owner]++] = (*m_kai_kbj.data)[ii + jj * npair] * fac;
                                }
                            }
                        }
                        i = i_next;
                    }
                    j = j_next;
                }
            }
        }
        for (int p = 0; p < GlobalV::NPROC; ++p)
        {
            assert(cursor_head[p] == shdispls[p] + send_head_counts[p]);
            assert(cursor_buffer[p] == sbdispls[p] + send_buffer_counts[p]);
        }
        // 3. communicate block heads and buffers
        MPI_Alltoallv(send_heads.data(), send_head_counts.data(), shdispls.data(), mpitype_blockhead,
                    recv_heads.data(), recv_head_counts.data(), rhdispls.data(), mpitype_blockhead,
                    pm_2d.comm());
        MPI_Alltoallv(send_buffers.data(), send_buffer_counts.data(), sbdispls.data(), LR_Util::MPIType<T>::value(),
                    recv_buffers.data(), recv_buffer_counts.data(), rbdispls.data(), LR_Util::MPIType<T>::value(),
                    pm_2d.comm());
        // 4. unpack recv head and buffer
        int buf_cursor = 0;
        for (int iblock = 0; iblock < recv_head_total; ++iblock)
        {
            const BlockHead& rh = recv_heads[iblock];
            const int nr = rh.nr;
            const int nc = rh.nc;
            const int lr = pm_2d.global2local_row(rh.gr0);
            const int lc = pm_2d.global2local_col(rh.gc0);
            const int lld = pm_2d.get_row_size();
            for (int j = 0; j < nc; ++j)
            {
                for (int i = 0; i < nr; ++i)
                {
                    const int idx_buffer = i + j * nr;
                    const std::size_t idx_2d = (lr + i) + std::size_t(lc + j) * std::size_t(lld);
                    m_2d[idx_2d] += recv_buffers[idx_buffer + buf_cursor];
                }
            }
            buf_cursor += nr * nc;
        }
        assert(buf_cursor == recv_buffer_total);
    }
#else
    // gather all Wk
    auto gather_matrix = [&](std::vector<T>& target,
        const std::valarray<T>& value,
        int k1_step,
        int k2_step,
        double factor) -> void
    {
        for (int j = 0; j < npair; ++j)
        {
            for (int i = 0; i < npair; ++i)
            {
                const std::size_t idx_target = (k1_step + i) + std::size_t(k2_step + j) * std::size_t(this->ndim);
                const int idx_value = i + j * npair;
                target[idx_target] += value[idx_value] * factor;
            }
        }
    };

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) collapse(2)
    #endif
    for (const int kai : this->k1_indices)
    {
        const int k1_step = kai * npair;
        for (const int kbj : this->k2_indices)
        {
            const int k2_step = kbj * npair;
            const RI::Tensor<T>& m_kai_kbj = m_lri.at(kai).at(kbj);
            gather_matrix(m_2d, *m_kai_kbj.data, k1_step, k2_step, fac_base * this->kv.wk[kai]);
        }
    }
#endif
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "transform_k_2dlocal");
    ModuleBase::timer::end("MolecularLRI", "transform_k_2dlocal");
}
}