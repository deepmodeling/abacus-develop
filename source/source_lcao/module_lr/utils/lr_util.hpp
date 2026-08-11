#pragma once
#include "lr_util.h"
#include "source_cell/unitcell.h"
#include "source_base/constants.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include <algorithm>
#include <cstddef>
#include <limits>
namespace LR_Util
{
    /// =================PHYSICS====================

    template <typename TCell>
    int cal_nelec(const TCell& ucell) {
        int nelec = 0;
        for (int it = 0; it < ucell.ntype; ++it) {
            nelec += ucell.atoms[it].ncpp.zv * ucell.atoms[it].na;
        }
        return nelec;
    }

    /// =================ALGORITHM====================
    ///for lack of make_unique in c++11 
    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args &&... args)
    {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    //====== newers and deleters========
    /// @brief  new 2d pointer
    /// @tparam T
    /// @param size1
    /// @param size2
    template <typename T>
    void _allocate_2order_nested_ptr(T**& p2, size_t size1, size_t size2)
    {
        p2 = new T * [size1];
        for (size_t i = 0; i < size1; ++i)
        {
            p2[i] = new T[size2];
        }
    };

    /// @brief  delete 2d pointer 
    /// @tparam T 
    /// @param p2 
    /// @param size 
    template <typename T>
    void _deallocate_2order_nested_ptr(T** p2, size_t size)
    {
        if (p2 != nullptr)
        {
            for (size_t i = 0; i < size; ++i)
            {
                if (p2[i] != nullptr) { delete[] p2[i]; }
            }
            delete[] p2;
        }
    };

    inline double get_conj(const double& x)
    {
        return x;
    }
    inline std::complex<double> get_conj(const std::complex<double>& x)
    {
        return std::conj(x);
    }
    template <typename T>
    void matsym(const T* in, const int n, T* out)
    {
        for (int i = 0; i < n; ++i) {
            out[i * n + i] = 0.5 * in[i * n + i] + 0.5 * get_conj(in[i * n + i]);
        }
        for (int i = 0;i < n;++i) {
            for (int j = i + 1;j < n;++j)
            {
                out[i * n + j] = 0.5 * (in[i * n + j] + get_conj(in[j * n + i]));
                out[j * n + i] = get_conj(out[i * n + j]);
            }
        }
    }
    template <typename T>
    void matsym(T* inout, const int n)
    {
        for (int i = 0; i < n; ++i) {
            inout[i * n + i] = 0.5 * (inout[i * n + i] + get_conj(inout[i * n + i]));
        }
        for (int i = 0;i < n;++i) {
            for (int j = i + 1;j < n;++j)
            {
                inout[i * n + j] = 0.5 * (inout[i * n + j] + get_conj(inout[j * n + i]));
                inout[j * n + i] = get_conj(inout[i * n + j]);
            }
        }
    }
    template<typename T>
    bool is_hermitian(const T* mat, const Parallel_2D& pmat, const double threshold){
        bool is_herm = true;
        int n = pmat.get_global_row_size();
        assert(n == pmat.get_global_col_size());
        double work_dummy = 0.0;
        const char norm_type = 'F';
        double norm1(0.0), norm2(0.0);
#ifdef __MPI
        const int one = 1;
        T alpha = 1.0;
        T beta = 0.0;
        if (pmat.get_local_size() > std::numeric_limits<int>::max())
        {
            std::cerr<< "Warning: pmat in RANK " + std::to_string(GlobalV::MY_RANK) + " has "
             + std::to_string(pmat.get_local_size()) + " elements, which may overflow during `tranc`, please use more mpi." << std::endl;
        }
        std::vector<T> herm_mat(pmat.get_local_size());
        ScalapackConnector::tranc(n, n,
                alpha, const_cast<T*>(mat), one, one, pmat.desc,
                beta, herm_mat.data(), one, one, pmat.desc);
        ModuleBase::TITLE("LR_Util", "is_hermitian: tranc");
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "is_hermitian: tranc");
        T diff;
        for (std::size_t i = 0;i < pmat.get_local_size();++i) {
            diff = mat[i] - herm_mat[i];
            if (std::abs(diff) > threshold) { is_herm = false; }
            norm1 += static_cast<double>(std::norm(diff));
            norm2 += static_cast<double>(std::norm(mat[i] + herm_mat[i]));
        }
        MPI_Allreduce(MPI_IN_PLACE, &norm1, 1, MPI_DOUBLE, MPI_SUM, pmat.comm());
        MPI_Allreduce(MPI_IN_PLACE, &norm2, 1, MPI_DOUBLE, MPI_SUM, pmat.comm());
        norm1 = std::sqrt(norm1);
        norm2 = std::sqrt(norm2);
        // We don't use `lange` here because minus_mat and sum_mat are memory-consuming.
        // double norm1 = ScalapackConnector::lange(norm_type, n, n, minus_mat.data(), one, one, pmat.desc, &work_dummy);
        // double norm2 = ScalapackConnector::lange(norm_type, n, n, sum_mat.data(), one, one, pmat.desc, &work_dummy);        
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "Frobenius norm");
        int local_flag = is_herm ? 1 : 0;
        int global_flag = 0;
        MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_LAND, pmat.comm());
        is_herm = (global_flag != 0);
#else
        assert(pmat.is_serial);
        std::vector<T> minus_mat(pmat.get_local_size());
        std::vector<T> sum_mat(pmat.get_local_size());
        for (std::size_t i = 0;i < n;++i) {
            for (std::size_t j = i;j < n;++j) {
                minus_mat[i * n + j] = mat[i * n + j] - get_conj(mat[j * n + i]);
                minus_mat[j * n + i] = -get_conj(minus_mat[i * n + j]);
                if (std::abs(minus_mat[i * n + j]) > threshold) { is_herm = false; }
                sum_mat[i * n + j] = mat[i * n + j] + get_conj(mat[j * n + i]);
                sum_mat[j * n + i] = get_conj(sum_mat[i * n + j]);
            }
        }
        norm1 = LapackConnector::lange(norm_type, n, n, minus_mat.data(), n, &work_dummy);
        norm2 = LapackConnector::lange(norm_type, n, n, sum_mat.data(), n, &work_dummy);
#endif
        if (GlobalV::MY_RANK == 0) {
            std::cout << "|  Hermitian check: ||A - A^H||_F = " << norm1 << ", ||A + A^H||_F = " << norm2 << std::endl;
            std::cout << "|   ||A - A^H||_F / ||A + A^H||_F = " << norm1 / norm2 << std::endl;
        }
        return is_herm;
    }

    template<typename T>
    bool is_symmetric(const T* mat, const Parallel_2D& pmat, const double threshold){
        bool is_sym = true;
        int n = pmat.get_global_row_size();
        assert(n == pmat.get_global_col_size());
        double work_dummy = 0.0;
        const char norm_type = 'F';
        double norm1(0.0), norm2(0.0);
#ifdef __MPI
        const int one = 1;
        T alpha = 1.0;
        T beta = 0.0;
        if (pmat.get_local_size() > std::numeric_limits<int>::max())
        {
            std::cerr<< "Warning: pmat in RANK " + std::to_string(GlobalV::MY_RANK) + " has "
             + std::to_string(pmat.get_local_size()) + " elements, which may overflow during `tranu`, please use more mpi." << std::endl;
        }
        std::vector<T> trans_mat(pmat.get_local_size());
        ScalapackConnector::tranu(n, n,
                alpha, const_cast<T*>(mat), one, one, pmat.desc,
                beta, trans_mat.data(), one, one, pmat.desc);
        ModuleBase::TITLE("LR_Util", "is_symmetric: tranu");
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "is_hermitian: tranu");
        T diff;
        for (std::size_t i = 0;i < pmat.get_local_size();++i) {
            diff = mat[i] - trans_mat[i];
            if (std::abs(diff) > threshold) { is_sym = false; }
            norm1 += static_cast<double>(std::norm(diff));
            norm2 += static_cast<double>(std::norm(mat[i] + trans_mat[i]));
        }
        MPI_Allreduce(MPI_IN_PLACE, &norm1, 1, MPI_DOUBLE, MPI_SUM, pmat.comm());
        MPI_Allreduce(MPI_IN_PLACE, &norm2, 1, MPI_DOUBLE, MPI_SUM, pmat.comm());
        norm1 = std::sqrt(norm1);
        norm2 = std::sqrt(norm2);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "Frobenius norm");
        int local_flag = is_sym ? 1 : 0;
        int global_flag = 0;
        MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_LAND, pmat.comm());
        is_sym = (global_flag != 0);
#else
        assert(pmat.is_serial);
        std::vector<T> minus_mat(pmat.get_local_size());
        std::vector<T> sum_mat(pmat.get_local_size());
        for (std::size_t i = 0;i < n;++i) {
            for (std::size_t j = i;j < n;++j) {
                minus_mat[i * n + j] = mat[i * n + j] - mat[j * n + i];
                minus_mat[j * n + i] = -minus_mat[i * n + j];
                if (std::abs(minus_mat[i * n + j]) > threshold) {is_sym = false; }
                sum_mat[i * n + j] = mat[i * n + j] + mat[j * n + i];
                sum_mat[j * n + i] = sum_mat[i * n + j];
            }
        }
        norm1 = LapackConnector::lange(norm_type, n, n, minus_mat.data(), n, &work_dummy);
        norm2 = LapackConnector::lange(norm_type, n, n, sum_mat.data(), n, &work_dummy);
#endif
        if (GlobalV::MY_RANK == 0) {
            std::cout << "|  Symmetric check: ||B - B^T||_F = " << norm1 << ", ||B + B^T||_F = " << norm2 << std::endl;
            std::cout << "|   ||B - B^T||_F / ||B + B^T||_F = " << norm1 / norm2 << std::endl;
        }
        return is_sym;
    }

    /// get the Psi wrapper of the selected spin from the Psi object
    template<typename T>
    psi::Psi<T> get_psi_spin(const psi::Psi<T>& psi_in, const int& is, const int& nk)
    {
        return psi::Psi<T>(&psi_in(is * nk, 0, 0), 
                           nk, 
                           psi_in.get_nbands(),
                           psi_in.get_nbasis(),
                           true);
    }

    /// psi(nk=1, nbands=nb, nk * nbasis) -> psi(nb, nk, nbasis) without memory copy
    template<typename T, typename Device>
    psi::Psi<T, Device> k1_to_bfirst_wrapper(const psi::Psi<T, Device>& psi_kfirst, int nk_in, int nbasis_in)
    {
        assert(psi_kfirst.get_nk() == 1);
        assert(nk_in * nbasis_in == psi_kfirst.get_nbasis());

        int ib_now = psi_kfirst.get_current_b();
        psi_kfirst.fix_b(0);    // for get_pointer() to get the head pointer
        psi::Psi<T, Device> psi_bfirst(psi_kfirst.get_pointer(), 
                                       nk_in, 
                                       psi_kfirst.get_nbands(), 
                                       nbasis_in, 
                                       nbasis_in, 
                                       false);
        psi_kfirst.fix_b(ib_now);
        return psi_bfirst;
    }

    ///  psi(nb, nk, nbasis) -> psi(nk=1, nbands=nb, nk * nbasis)  without memory copy
    template<typename T, typename Device>
    psi::Psi<T, Device> bfirst_to_k1_wrapper(const psi::Psi<T, Device>& psi_bfirst)
    {
        int ib_now = psi_bfirst.get_current_b();
        int ik_now = psi_bfirst.get_current_k();

        psi_bfirst.fix_kb(0, 0);    // for get_pointer() to get the head pointer
        psi::Psi<T, Device> psi_kfirst(psi_bfirst.get_pointer(), 
                                       1, 
                                       psi_bfirst.get_nbands(), 
                                       psi_bfirst.get_nk() * psi_bfirst.get_nbasis(), 
                                       psi_bfirst.get_nk() * psi_bfirst.get_nbasis(),
                                       true);
        psi_bfirst.fix_kb(ik_now, ib_now);
        return psi_kfirst;
    }
//=================2D-block Parallel===============

    /// @brief assign global X to 2d-matrix, its col is band(excition state), and row is { spin, k-point, occ, virt }
    /// @attention pX is 2d-blocked as {occ, virt}, this assignment is used to calculate transition density matrix c_b X_{bj} c_j
    /// @todo this function is a merge version of HamiltULR::global2local and HamiltLR::global2local, they should be replaced
    template <typename T>
    void global2local_X(T* local_X, const T* global_X, const int nband, const int nk, 
        const std::vector<int>& nocc, const std::vector<int>& nvirt, const std::vector<Parallel_2D>& pX,
        const bool openshell)
    {
        const int nspin_X = openshell ? 2 : 1;
        const std::vector<int> npairs = { nocc[0] * nvirt[0], nocc[1] * nvirt[1] };
        const int gdim = openshell ? nk * (npairs[0] + npairs[1] ) : nk * npairs[0];
        const int ldim = openshell ? nk * (pX[0].get_local_size() + pX[1].get_local_size()) : nk * pX[0].get_local_size();
        
        for (int ib = 0;ib < nband;++ib)
        {
            const int loffset_b = ib * ldim;
            const int goffset_b = ib * gdim;
            for (int is = 0;is < nspin_X;++is)
            {
                const int loffset_bs = loffset_b + is * nk * pX[0].get_local_size();
                const int goffset_bs = goffset_b + is * nk * npairs[0];                    
                for (int ik = 0;ik < nk;++ik)
                {
                    const int loffset = loffset_bs + ik * pX[is].get_local_size();
                    const int goffset = goffset_bs + ik * npairs[is];
                    for (int lo = 0;lo < pX[is].get_col_size();++lo)
                    {
                        const int go = pX[is].local2global_col(lo);
                        for (int lv = 0;lv < pX[is].get_row_size();++lv)
                        {
                            const int gv = pX[is].local2global_row(lv);
                            local_X[loffset + lo * pX[is].get_row_size() + lv] = global_X[goffset + go * nvirt[is] + gv];
                        }
                    }
                }
            }
        }
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "global2local_X");
    }
#ifdef __MPI
    struct BlockHead
    {
        int ivirt;    // {nvirt, nocc} global row
        int iocc;    // {nvirt, nocc} global col
        int ik;
        int src_rank;
        void set(int iv_, int io_, int ik_, int src_rank_)
        {
            ivirt = iv_;
            iocc = io_;
            ik = ik_;
            src_rank = src_rank_;
        }
    };
    inline MPI_Datatype mpi_type_blockhead()
    {
        static MPI_Datatype dt = MPI_DATATYPE_NULL;
        static bool committed = false;
        if (!committed)
        {
            int blen[4] = {1, 1, 1, 1};
            MPI_Aint disp[4];
            MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    
            disp[0] = static_cast<MPI_Aint>(offsetof(BlockHead, ivirt));
            disp[1] = static_cast<MPI_Aint>(offsetof(BlockHead, iocc));
            disp[2] = static_cast<MPI_Aint>(offsetof(BlockHead, ik));
            disp[3] = static_cast<MPI_Aint>(offsetof(BlockHead, src_rank));
    
            MPI_Type_create_struct(4, blen, disp, types, &dt);
            MPI_Type_commit(&dt);
            committed = true;
        }
        return dt;
    }

    /// @brief assign X in pA to X in pX
    /// @note if use Cpxgemr2d, the communication time will be much more than MPIAlltoall by hand
    template <typename T>
    void pA2pX(T* X_pX, const T* X_pA, const int nband, const int nk, 
        const std::vector<int>& nocc, const std::vector<int>& nvirt,
        const std::vector<Parallel_2D>& pX, const Parallel_2D& pA,
        const int row_offset, const int col_offset, const bool openshell)
    {
        ModuleBase::TITLE("LR_Util", "pA2pX");
        ModuleBase::timer::start("LR_Util", "pA2pX");
        const int nspin_X = openshell ? 2 : 1;
        const std::vector<int> npairs = { nocc[0] * nvirt[0], nocc[1] * nvirt[1] };
        const std::size_t gdim = openshell ? nk * (npairs[0] + npairs[1] ) : nk * npairs[0];
        const std::size_t ldim = openshell ? nk * (pX[0].get_local_size() + pX[1].get_local_size()) : nk * pX[0].get_local_size();
        assert(pA.get_global_row_size() == gdim || pA.get_global_row_size() == 2 * gdim);
        assert(pA.get_global_col_size() == gdim || pA.get_global_col_size() == 2 * gdim);
        for (int is = 0;is < nspin_X;++is)
        {
            assert(pX[is].get_global_row_size() == nvirt[is]);
            assert(pX[is].get_global_col_size() == nocc[is]);
        }
        MPI_Datatype mpitype_blockhead = mpi_type_blockhead();

        // 0. outer loop of row, communicate per 64 kpoints
        constexpr int comm_nk = 64;
        std::vector<int> send_head_counts(GlobalV::NPROC, 0), recv_head_counts(GlobalV::NPROC, 0);
        std::vector<int> send_buffer_counts(GlobalV::NPROC, 0), recv_buffer_counts(GlobalV::NPROC, 0);
        std::vector<int> shdispls(GlobalV::NPROC, 0), rhdispls(GlobalV::NPROC, 0); // displacements of block heads
        std::vector<int> sbdispls(GlobalV::NPROC, 0), rbdispls(GlobalV::NPROC, 0); // displacements of buffer
        std::vector<int> cursor_head(GlobalV::NPROC, 0), cursor_buffer(GlobalV::NPROC, 0);
        std::vector<BlockHead> send_heads, recv_heads;
        std::vector<T> send_buffers, recv_buffers;
        int max_send_head_total = 0, max_recv_head_total = 0;
        int max_send_buffer_total = 0, max_recv_buffer_total = 0;

        // 1. pre-calculate local bands and global band list on each processor
        const int gb_start = col_offset;
        const int gb_end = col_offset + nband;
        std::vector<int> lnbands(GlobalV::NPROC, 0);
        std::vector<std::vector<int>> gblist(GlobalV::NPROC);
        const int nbpA = pA.get_block_size();
        for (int gb = gb_start; gb < gb_end; ++gb)
        {
            int proc_col = (gb / nbpA) % pA.dim1;
            for (int pr = 0; pr < pA.dim0; ++pr)
            {
                int ownerA = proc_col + pr * pA.dim1;
                ++lnbands[ownerA];
                gblist[ownerA].push_back(gb);
            }
        }
        const int lnband = lnbands[GlobalV::MY_RANK];
        if (lnband == 0)
        {
            throw std::runtime_error("rank" + std::to_string(GlobalV::MY_RANK) + "lnband = 0, please check the pA and pX distribution!");
        }
        const int lb_start = pA.global2local_col(gblist[GlobalV::MY_RANK].front());
        const int lb_end = lb_start + lnband;        

        for (int is = 0; is < nspin_X; ++is)
        {
            const std::size_t gspin_base = is * nk * npairs[0];
            const std::size_t lspin_base = is * nk * pX[0].get_local_size();
            const std::size_t nrow = comm_nk * npairs[is];
            const std::size_t grow_start = row_offset + gspin_base;
            const std::size_t grow_end = grow_start + nk * npairs[is];
            for (int irow_start = grow_start; irow_start < grow_end; irow_start += nrow)
            {
                std::fill(send_head_counts.begin(), send_head_counts.end(), 0);
                std::fill(send_buffer_counts.begin(), send_buffer_counts.end(), 0);
                
                // 2. calculate and coummunicate block counts, then calculate block displacements
                int irow_end = std::min(irow_start + nrow, grow_end);
                for (int irow = irow_start; irow < irow_end; ++irow)
                {
                    const int lrpA = pA.global2local_row(irow);
                    if (lrpA == -1) continue;
                    const int ir = irow - row_offset - gspin_base;
                    const int ik = ir / npairs[is];
                    const int io = (ir-ik*npairs[is]) / nvirt[is];
                    const int iv = ir % nvirt[is];
                    int ownerX = pX[is].owner_processor(iv, io);
                    if (ownerX != GlobalV::MY_RANK)
                    {
                        ++send_head_counts[ownerX];
                    }
                }
                for (int p = 0; p < GlobalV::NPROC; ++p)
                {
                    std::size_t nbuffer = size_t(lnband) * send_head_counts[p];
                    if (nbuffer > std::numeric_limits<int>::max())
                    {
                        throw std::overflow_error("in pA2pX: overflow converting to int!");
                    }
                    send_buffer_counts[p] = static_cast<int>(nbuffer);
                }

                assert(send_head_counts.at(GlobalV::MY_RANK) == 0);
                assert(send_buffer_counts.at(GlobalV::MY_RANK) == 0);
                MPI_Alltoall(send_head_counts.data(), 1, MPI_INT, recv_head_counts.data(), 1, MPI_INT, pA.comm());
                MPI_Alltoall(send_buffer_counts.data(), 1, MPI_INT,
                             recv_buffer_counts.data(), 1, MPI_INT, pA.comm());

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
                
                // 3. prepare block heads and buffers
                send_heads.resize(send_head_total);
                recv_heads.resize(recv_head_total);
                send_buffers.resize(send_buffer_total);
                recv_buffers.resize(recv_buffer_total);
                cursor_head = shdispls;
                cursor_buffer = sbdispls;

                for (int irow = irow_start; irow < irow_end; ++irow)
                {
                    const int lr = pA.global2local_row(irow);
                    if (lr == -1) continue;
                    const std::size_t ir = irow - row_offset - gspin_base;
                    const std::size_t ik = ir / npairs[is];
                    const int io = (ir-ik*npairs[is]) / nvirt[is];
                    const int iv = ir % nvirt[is];
                    const int ownerX = pX[is].owner_processor(iv, io);
                    if (ownerX == GlobalV::MY_RANK)
                    {
                        for (std::size_t lb = lb_start; lb < lb_end; ++lb)
                        {
                            const int gb = pA.local2global_col(lb);
                            const std::size_t lb_base = std::size_t(gb - col_offset) * std::size_t(ldim);
                            const std::size_t lrowX = pX[is].global2local_row(iv);
                            const std::size_t lcolX = pX[is].global2local_col(io);
                            const std::size_t lX_indx = lrowX + lcolX * pX[is].get_row_size() + lb_base 
                                                        + lspin_base + ik * pX[is].get_local_size();
                            X_pX[lX_indx] = X_pA[static_cast<std::size_t>(lr) + lb * pA.get_row_size()];
                        }
                    }
                    else
                    {
                        BlockHead& head = send_heads[cursor_head[ownerX]++];
                        head.set(iv, io, ik, GlobalV::MY_RANK);
                        for (std::size_t lb = lb_start; lb < lb_end; ++lb)
                        {
                            send_buffers[cursor_buffer[ownerX]++] = X_pA[static_cast<std::size_t>(lr) + lb * pA.get_row_size()];
                        }
                    }
                }
                for (int p = 0; p < GlobalV::NPROC; ++p)
                {
                    assert(cursor_head[p] == shdispls[p] + send_head_counts[p]);
                    assert(cursor_buffer[p] == sbdispls[p] + send_buffer_counts[p]);
                }
                // 4. communicate block heads and buffers
                MPI_Alltoallv(send_heads.data(), send_head_counts.data(), shdispls.data(), mpitype_blockhead,
                              recv_heads.data(), recv_head_counts.data(), rhdispls.data(), mpitype_blockhead,
                              pA.comm());
                MPI_Alltoallv(send_buffers.data(), send_buffer_counts.data(), sbdispls.data(), LR_Util::MPIType<T>::value(),
                              recv_buffers.data(), recv_buffer_counts.data(), rbdispls.data(), LR_Util::MPIType<T>::value(),
                              pA.comm());

                // 5. unpack received buffers to pX
                for (int src = 0; src < GlobalV::NPROC; ++src)
                {
                    const int nh = recv_head_counts[src];
                    if (nh == 0) continue;
                    int buf_cursor = rbdispls[src];
                    const int ib_start = rhdispls[src];
                    const int ib_end = ib_start + nh;
                    for (int iblock = ib_start; iblock < ib_end; ++iblock)
                    {
                        const BlockHead& rh = recv_heads[iblock];
                        assert(rh.src_rank == src);
                        const std::size_t ik = rh.ik;
                        const std::size_t lr = pX[is].global2local_row(rh.ivirt);
                        const std::size_t lc = pX[is].global2local_col(rh.iocc);
                        for (int ib = 0; ib < lnbands[src]; ++ib)
                        {
                            const std::size_t lb_base = std::size_t(gblist[src][ib] - col_offset) * std::size_t(ldim);
                            const std::size_t lX_indx = lr + lc * pX[is].get_row_size() + lb_base 
                                                        + lspin_base + ik * pX[is].get_local_size();
                            X_pX[lX_indx] = recv_buffers[ib + buf_cursor];
                        }
                        buf_cursor += lnbands[src];
                    }
                    assert(buf_cursor == rbdispls[src] + recv_buffer_counts[src]);
                }
            } // end of irow_start
        } // end of is
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "pA2pX");
        ModuleBase::timer::end("LR_Util", "pA2pX");
    }

    template <typename T>
    void gather_2d_to_full(const Parallel_2D& pv, const T* submat, T* fullmat,
        const bool row_major, const std::size_t global_nrow, const std::size_t global_ncol)
    {
        assert(pv.get_global_row_size() == global_nrow);
        assert(pv.get_global_col_size() == global_ncol);
        // copy
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int j = 0;j < pv.get_col_size();++j) {
            for (int i = 0;i < pv.get_row_size();++i) {
                if (row_major) {
                    fullmat[pv.local2global_row(i) * global_ncol + pv.local2global_col(j)] = submat[i * pv.get_col_size() + j];
                } else {
                    fullmat[pv.local2global_col(j) * global_nrow + pv.local2global_row(i)] = submat[j * pv.get_row_size() + i];
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, fullmat, global_nrow * global_ncol, LR_Util::MPIType<T>::value(), MPI_SUM, pv.comm());
    };
#endif

}
