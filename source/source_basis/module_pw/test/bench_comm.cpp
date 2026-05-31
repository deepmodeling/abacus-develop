/**
 * @file bench_comm.cpp
 * @brief MPI microbenchmark: blocking (MPI_Alltoallv) vs nonblocking
 *        (MPI_Isend/Irecv + MPI_Waitsome) gather/scatter.
 *
 * This directly mirrors ABACUS's feat/unblock gatherp_scatters/gathers_scatterp.
 *
 * Compile from abacus-develop root:
 *   /usr/bin/mpicxx.openmpi -std=c++14 -O3 -fopenmp -DOMPI_SKIP_MPICXX \
 *     -o bench_comm source/source_basis/module_pw/test/bench_comm.cpp
 *
 * Run:
 *   mpirun -np N ./bench_comm [n_iter] [nx] [ny] [nz]
 */
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

// ---------------------------------------------------------------------------
struct CommParams {
    int poolnproc, poolrank;
    std::vector<int> numz, nst_per, numg, numr, startg, startr, startz;
    int send_total, recv_total, nst, nz, nplane, fftnxy, nstot, nrxx;
    MPI_Comm comm;
};

CommParams generate_params(int nx, int ny, int nz, MPI_Comm comm) {
    CommParams p;
    p.comm = comm;
    MPI_Comm_size(comm, &p.poolnproc);
    MPI_Comm_rank(comm, &p.poolrank);

    p.numz.resize(p.poolnproc);
    int base_z = nz / p.poolnproc;
    int rem_z  = nz % p.poolnproc;
    for (int ip = 0; ip < p.poolnproc; ++ip)
        p.numz[ip] = base_z + (ip < rem_z ? 1 : 0);

    int fftnx = nx / 2 + 1;
    p.fftnxy = fftnx * ny;
    p.nplane = p.numz[p.poolrank];
    p.nrxx = p.nplane * p.fftnxy;
    p.nz = nz;
    p.nstot = p.fftnxy;
    p.nst_per.resize(p.poolnproc);
    int base_s = p.nstot / p.poolnproc;
    int rem_s  = p.nstot % p.poolnproc;
    for (int ip = 0; ip < p.poolnproc; ++ip)
        p.nst_per[ip] = base_s + (ip < rem_s ? 1 : 0);
    p.nst = p.nst_per[p.poolrank];

    p.numg.resize(p.poolnproc);
    p.numr.resize(p.poolnproc);
    p.startg.resize(p.poolnproc);
    p.startr.resize(p.poolnproc);
    p.startz.resize(p.poolnproc);
    for (int ip = 0; ip < p.poolnproc; ++ip) {
        p.numg[ip] = p.nst_per[p.poolrank] * p.numz[ip];
        p.numr[ip] = p.nst_per[ip] * p.numz[p.poolrank];
    }
    p.startg[0] = 0; p.startr[0] = 0; p.startz[0] = 0;
    for (int ip = 1; ip < p.poolnproc; ++ip) {
        p.startg[ip] = p.startg[ip - 1] + p.numg[ip - 1];
        p.startr[ip] = p.startr[ip - 1] + p.numr[ip - 1];
        p.startz[ip] = p.startz[ip - 1] + p.numz[ip - 1];
    }
    p.send_total = p.startr.back() + p.numr.back();
    p.recv_total = p.startg.back() + p.numg.back();
    return p;
}

// =========================================================================
//  Blocking (original ABACUS develop-branch pattern)
// =========================================================================

template <typename T>
static void blocking_gatherp(const CommParams& p,
                             std::complex<T>* in,
                             std::complex<T>* out,
                             std::complex<T>* workbuf)
{
    if (p.poolnproc == 1) {
        for (int is = 0; is < p.nst; ++is) {
            std::complex<T>* outp = &out[is * p.nz];
            std::complex<T>* inp  = &in[is * p.nz];
            for (int iz = 0; iz < p.nz; ++iz) outp[iz] = inp[iz];
        }
        return;
    }

    std::complex<T>* sendbuf = workbuf;

    // Pack: (nplane, fftnxy) -> (nplane, nstot) in sendbuf
    if (p.nplane > 0) {
#pragma omp parallel for
        for (int istot = 0; istot < p.nstot; ++istot) {
            int ixy = istot;
            std::complex<T>* outp = &sendbuf[istot * p.nplane];
            std::complex<T>* inp  = &in[ixy * p.nplane];
            for (int iz = 0; iz < p.nplane; ++iz) outp[iz] = inp[iz];
        }
    }

    MPI_Datatype mpi_type = (sizeof(T) == sizeof(double))
        ? MPI_DOUBLE_COMPLEX : MPI_COMPLEX;

    // Blocking Alltoallv: send from sendbuf, receive into out (reuse out as recv)
    MPI_Alltoallv(sendbuf, p.numr.data(), p.startr.data(), mpi_type,
                  out,     p.numg.data(), p.startg.data(), mpi_type,
                  p.comm);

    // Unpack: out now has data in (numz[ip], nst) layout, reorganize to (nz, nst)
#pragma omp parallel for collapse(2)
    for (int ip = 0; ip < p.poolnproc; ++ip) {
        for (int is = 0; is < p.nst; ++is) {
            int nzip = p.numz[ip];
            // out temporarily holds received data; unpack in-place is tricky.
            // The ABACUS original uses 'in' (sticks) as recv and unpacks to 'out'.
            // Since we use a workbuf, we receive into 'out' directly in the correct
            // layout, so no second unpack needed — but we DO need to rearrange
            // from (numz[ip],nst) to (nz,nst).
            std::complex<T>* outp = &out[p.startz[ip] + is * p.nz];
            std::complex<T>* inp  = &out[p.startg[ip] + is * nzip];
            for (int izip = 0; izip < nzip; ++izip) outp[izip] = inp[izip];
        }
    }
}

template <typename T>
static void blocking_gathers(const CommParams& p,
                             std::complex<T>* in,
                             std::complex<T>* out,
                             std::complex<T>* workbuf)
{
    if (p.poolnproc == 1) {
        for (int i = 0; i < p.nrxx; ++i) out[i] = std::complex<T>(0, 0);
        for (int is = 0; is < p.nst; ++is) {
            std::complex<T>* outp = &out[is * p.nz];
            std::complex<T>* inp  = &in[is * p.nz];
            for (int iz = 0; iz < p.nz; ++iz) outp[iz] = inp[iz];
        }
        return;
    }

    std::complex<T>* sendbuf = workbuf;

    // Pack: (nz, nst) -> (numz[ip], nst) in sendbuf
#pragma omp parallel for collapse(2)
    for (int ip = 0; ip < p.poolnproc; ++ip) {
        for (int is = 0; is < p.nst; ++is) {
            int nzip = p.numz[ip];
            std::complex<T>* outp = &sendbuf[p.startg[ip] + is * nzip];
            std::complex<T>* inp  = &in[p.startz[ip] + is * p.nz];
            for (int izip = 0; izip < nzip; ++izip) outp[izip] = inp[izip];
        }
    }

    MPI_Datatype mpi_type = (sizeof(T) == sizeof(double))
        ? MPI_DOUBLE_COMPLEX : MPI_COMPLEX;

    // Blocking Alltoallv: send from sendbuf, receive into out
    MPI_Alltoallv(sendbuf, p.numg.data(), p.startg.data(), mpi_type,
                  out,     p.numr.data(), p.startr.data(), mpi_type,
                  p.comm);

    // Unpack: out has received data in (numr layout), repack to (nplane, fftnxy)
    // But first we need to preserve the received data before clearing out
    // Actually in the ABACUS original, the recv was into 'in', not 'out'.
    // We'll use a different approach: unpack from the received 'out' data.
    // The received data in 'out' is at offsets startr[ip], in plane-major layout.
    // We need to copy it to the correct (ixy, nplane) positions.

    // Since we received directly into out, and the zeroing + unpack would lose data,
    // let's do a careful unpack: out[startr[ip] + is*nplane] -> out[ixy * nplane + iz]
    // But this is an in-place operation on overlapping regions, so we need a temp buffer.

    // Simpler: use a temp buffer for the received data
    std::vector<std::complex<T>> temp(p.send_total);
    for (int i = 0; i < p.send_total; ++i) temp[i] = out[i];

    // Zero output
#pragma omp parallel for schedule(static)
    for (int i = 0; i < p.nrxx; ++i) out[i] = std::complex<T>(0, 0);

    // Unpack
    std::vector<int> istot_offsets(p.poolnproc, 0);
    for (int ip = 1; ip < p.poolnproc; ++ip)
        istot_offsets[ip] = istot_offsets[ip - 1] + p.nst_per[ip - 1];

#pragma omp parallel for
    for (int ip = 0; ip < p.poolnproc; ++ip) {
        int peer_nst = p.nst_per[ip];
        if (peer_nst == 0 || p.nplane == 0) continue;
        int istot0 = istot_offsets[ip];
        for (int is = 0; is < peer_nst; ++is) {
            int istot = istot0 + is;
            int ixy = istot;
            std::complex<T>* outp = &out[ixy * p.nplane];
            std::complex<T>* inp  = &temp[p.startr[ip] + is * p.nplane];
            for (int iz = 0; iz < p.nplane; ++iz) outp[iz] = inp[iz];
        }
    }
}

// =========================================================================
//  Nonblocking (ABACUS feat/unblock branch pattern)
// =========================================================================

template <typename T>
static void nonblocking_gatherp(const CommParams& p,
                                std::complex<T>* in,
                                std::complex<T>* out,
                                std::complex<T>* workbuf)
{
    if (p.poolnproc == 1) {
        for (int is = 0; is < p.nst; ++is) {
            std::complex<T>* outp = &out[is * p.nz];
            std::complex<T>* inp  = &in[is * p.nz];
            for (int iz = 0; iz < p.nz; ++iz) outp[iz] = inp[iz];
        }
        return;
    }

    std::complex<T>* sendbuf = workbuf;
    std::complex<T>* recvbuf = workbuf + p.send_total;

    // Pack
    if (p.nplane > 0) {
#pragma omp parallel for
        for (int istot = 0; istot < p.nstot; ++istot) {
            int ixy = istot;
            std::complex<T>* outp = &sendbuf[istot * p.nplane];
            std::complex<T>* inp  = &in[ixy * p.nplane];
            for (int iz = 0; iz < p.nplane; ++iz) outp[iz] = inp[iz];
        }
    }

    MPI_Datatype mpi_type = (sizeof(T) == sizeof(double))
        ? MPI_DOUBLE_COMPLEX : MPI_COMPLEX;

    std::vector<MPI_Request> recv_req(p.poolnproc, MPI_REQUEST_NULL);
    std::vector<MPI_Request> send_req(p.poolnproc, MPI_REQUEST_NULL);
    int active_recvs = 0, active_sends = 0;

    for (int ip = 0; ip < p.poolnproc; ++ip) {
        if (ip == p.poolrank || p.numg[ip] == 0) continue;
        MPI_Irecv(&recvbuf[p.startg[ip]], p.numg[ip], mpi_type,
                  ip, 0, p.comm, &recv_req[ip]);
        ++active_recvs;
    }
    for (int ip = 0; ip < p.poolnproc; ++ip) {
        if (ip == p.poolrank || p.numr[ip] == 0) continue;
        MPI_Isend(&sendbuf[p.startr[ip]], p.numr[ip], mpi_type,
                  ip, 0, p.comm, &send_req[ip]);
        ++active_sends;
    }

    // Self-copy
    for (int i = 0; i < p.numg[p.poolrank]; ++i)
        recvbuf[p.startg[p.poolrank] + i] = sendbuf[p.startr[p.poolrank] + i];

    auto unpack = [&](int ip) {
        int nzip = p.numz[ip];
#pragma omp parallel for
        for (int is = 0; is < p.nst; ++is) {
            std::complex<T>* outp = &out[is * p.nz + p.startz[ip]];
            std::complex<T>* inp  = &recvbuf[p.startg[ip] + is * nzip];
            for (int izip = 0; izip < nzip; ++izip) outp[izip] = inp[izip];
        }
    };
    unpack(p.poolrank);

    std::vector<MPI_Status> recv_status(p.poolnproc);
    std::vector<int> recv_indices(p.poolnproc, MPI_UNDEFINED);
    while (active_recvs > 0) {
        int outcount = 0;
        MPI_Waitsome(p.poolnproc, recv_req.data(), &outcount,
                     recv_indices.data(), recv_status.data());
        if (outcount == MPI_UNDEFINED) break;
        for (int idx = 0; idx < outcount; ++idx) unpack(recv_indices[idx]);
        active_recvs -= outcount;
    }
    if (active_sends > 0)
        MPI_Waitall(p.poolnproc, send_req.data(), MPI_STATUSES_IGNORE);
}

template <typename T>
static void nonblocking_gathers(const CommParams& p,
                                std::complex<T>* in,
                                std::complex<T>* out,
                                std::complex<T>* workbuf)
{
    if (p.poolnproc == 1) {
        for (int i = 0; i < p.nrxx; ++i) out[i] = std::complex<T>(0, 0);
        for (int is = 0; is < p.nst; ++is) {
            std::complex<T>* outp = &out[is * p.nz];
            std::complex<T>* inp  = &in[is * p.nz];
            for (int iz = 0; iz < p.nz; ++iz) outp[iz] = inp[iz];
        }
        return;
    }

    int send_count = p.recv_total;
    std::complex<T>* sendbuf = workbuf;
    std::complex<T>* recvbuf = workbuf + send_count;

    // Pack
#pragma omp parallel for collapse(2)
    for (int ip = 0; ip < p.poolnproc; ++ip) {
        for (int is = 0; is < p.nst; ++is) {
            int nzip = p.numz[ip];
            std::complex<T>* outp = &sendbuf[p.startg[ip] + is * nzip];
            std::complex<T>* inp  = &in[p.startz[ip] + is * p.nz];
            for (int izip = 0; izip < nzip; ++izip) outp[izip] = inp[izip];
        }
    }

    MPI_Datatype mpi_type = (sizeof(T) == sizeof(double))
        ? MPI_DOUBLE_COMPLEX : MPI_COMPLEX;

    std::vector<MPI_Request> recv_req(p.poolnproc, MPI_REQUEST_NULL);
    std::vector<MPI_Request> send_req(p.poolnproc, MPI_REQUEST_NULL);
    int active_recvs = 0, active_sends = 0;

    for (int ip = 0; ip < p.poolnproc; ++ip) {
        if (ip == p.poolrank || p.numr[ip] == 0) continue;
        MPI_Irecv(&recvbuf[p.startr[ip]], p.numr[ip], mpi_type,
                  ip, 0, p.comm, &recv_req[ip]);
        ++active_recvs;
    }
    for (int ip = 0; ip < p.poolnproc; ++ip) {
        if (ip == p.poolrank || p.numg[ip] == 0) continue;
        MPI_Isend(&sendbuf[p.startg[ip]], p.numg[ip], mpi_type,
                  ip, 0, p.comm, &send_req[ip]);
        ++active_sends;
    }

    // Zero output
#pragma omp parallel for schedule(static)
    for (int i = 0; i < p.nrxx; ++i) out[i] = std::complex<T>(0, 0);

    // Self-copy
    for (int i = 0; i < p.numr[p.poolrank]; ++i)
        recvbuf[p.startr[p.poolrank] + i] = sendbuf[p.startg[p.poolrank] + i];

    std::vector<int> istot_offsets(p.poolnproc, 0);
    for (int ip = 1; ip < p.poolnproc; ++ip)
        istot_offsets[ip] = istot_offsets[ip - 1] + p.nst_per[ip - 1];

    auto unpack = [&](int ip) {
        int peer_nst = p.nst_per[ip];
        if (peer_nst == 0 || p.nplane == 0) return;
        int istot0 = istot_offsets[ip];
#pragma omp parallel for
        for (int is = 0; is < peer_nst; ++is) {
            int istot = istot0 + is;
            int ixy = istot;
            std::complex<T>* outp = &out[ixy * p.nplane];
            std::complex<T>* inp  = &recvbuf[p.startr[ip] + is * p.nplane];
            for (int iz = 0; iz < p.nplane; ++iz) outp[iz] = inp[iz];
        }
    };
    unpack(p.poolrank);

    std::vector<MPI_Status> recv_status(p.poolnproc);
    std::vector<int> recv_indices(p.poolnproc, MPI_UNDEFINED);
    while (active_recvs > 0) {
        int outcount = 0;
        MPI_Waitsome(p.poolnproc, recv_req.data(), &outcount,
                     recv_indices.data(), recv_status.data());
        if (outcount == MPI_UNDEFINED) break;
        for (int idx = 0; idx < outcount; ++idx) unpack(recv_indices[idx]);
        active_recvs -= outcount;
    }
    if (active_sends > 0)
        MPI_Waitall(p.poolnproc, send_req.data(), MPI_STATUSES_IGNORE);
}

// =========================================================================
//  Benchmark driver
// =========================================================================

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n_iter = 1000;
    int nx = 32, ny = 32, nz = 32;
    if (argc > 1) n_iter = atoi(argv[1]);
    if (argc > 2) nx = atoi(argv[2]);
    if (argc > 3) ny = atoi(argv[3]);
    if (argc > 4) nz = atoi(argv[4]);

    CommParams p = generate_params(nx, ny, nz, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Grid: %dx%dx%d  Iterations: %d  MPI ranks: %d  OMP threads: %d\n",
               nx, ny, nz, n_iter, size, omp_get_max_threads());
        printf("Buffer: send_total=%d  recv_total=%d  nst=%d  nz=%d  nplane=%d\n",
               p.send_total, p.recv_total, p.nst, p.nz, p.nplane);
    }

    // Allocate buffers
    int64_t work_size = (p.send_total + p.recv_total) * 2;
    if (work_size < 2) work_size = 2;
    int64_t plane_size = (p.nrxx > 0) ? p.nrxx : 1;
    int64_t sticks_size = (p.nst * p.nz > 0) ? p.nst * p.nz : 1;

    std::vector<std::complex<double>> workbuf(work_size);
    std::vector<std::complex<double>> plane_in(plane_size);
    std::vector<std::complex<double>> sticks(sticks_size);
    std::vector<std::complex<double>> plane_out(plane_size);

    // Fill input data
    for (int64_t i = 0; i < plane_size; ++i)
        plane_in[i] = std::complex<double>(sin(i * 0.01), cos(i * 0.01));

    int warmup = std::max(10, n_iter / 10);

    // Warmup nonblocking
    for (int i = 0; i < warmup; ++i) {
        nonblocking_gatherp<double>(p, plane_in.data(), sticks.data(), workbuf.data());
        nonblocking_gathers<double>(p, sticks.data(), plane_out.data(), workbuf.data());
    }

    // Warmup blocking
    for (int i = 0; i < warmup; ++i) {
        blocking_gatherp<double>(p, plane_in.data(), sticks.data(), workbuf.data());
        blocking_gathers<double>(p, sticks.data(), plane_out.data(), workbuf.data());
    }

    // ---- TIMED RUNS ----
    MPI_Barrier(MPI_COMM_WORLD);

    // Blocking gatherp
    double t0 = MPI_Wtime();
    for (int i = 0; i < n_iter; ++i)
        blocking_gatherp<double>(p, plane_in.data(), sticks.data(), workbuf.data());
    double t_blk_fwd = MPI_Wtime() - t0;

    // Nonblocking gatherp
    t0 = MPI_Wtime();
    for (int i = 0; i < n_iter; ++i)
        nonblocking_gatherp<double>(p, plane_in.data(), sticks.data(), workbuf.data());
    double t_nb_fwd = MPI_Wtime() - t0;

    // Blocking gathers
    t0 = MPI_Wtime();
    for (int i = 0; i < n_iter; ++i)
        blocking_gathers<double>(p, sticks.data(), plane_out.data(), workbuf.data());
    double t_blk_rev = MPI_Wtime() - t0;

    // Nonblocking gathers
    t0 = MPI_Wtime();
    for (int i = 0; i < n_iter; ++i)
        nonblocking_gathers<double>(p, sticks.data(), plane_out.data(), workbuf.data());
    double t_nb_rev = MPI_Wtime() - t0;

    // Gather timing stats
    double buf[8], rbuf[8];
    buf[0] = t_blk_fwd; buf[1] = t_nb_fwd;
    buf[2] = t_blk_rev; buf[3] = t_nb_rev;
    buf[4] = t_blk_fwd; buf[5] = t_nb_fwd;  // for max
    buf[6] = t_blk_rev; buf[7] = t_nb_rev;

    MPI_Reduce(buf, rbuf, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(buf + 4, rbuf + 4, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double fwd_blk = rbuf[0] / size / n_iter;
        double fwd_nb  = rbuf[1] / size / n_iter;
        double rev_blk = rbuf[2] / size / n_iter;
        double rev_nb  = rbuf[3] / size / n_iter;
        double tot_blk = fwd_blk + rev_blk;
        double tot_nb  = fwd_nb + rev_nb;

        printf("\n--- Average per-iteration times (us) ---\n");
        printf("%-22s %10s %10s %8s\n", "Operation", "Blocking", "Nonblock", "Speedup");
        printf("gatherp_scatters      %10.2f %10.2f %7.2fx\n",
               fwd_blk * 1e6, fwd_nb * 1e6, fwd_blk / fwd_nb);
        printf("gathers_scatterp      %10.2f %10.2f %7.2fx\n",
               rev_blk * 1e6, rev_nb * 1e6, rev_blk / rev_nb);
        printf("Total roundtrip       %10.2f %10.2f %7.2fx\n",
               tot_blk * 1e6, tot_nb * 1e6, tot_blk / tot_nb);

        double fwd_blk_max = rbuf[4] / n_iter;
        double fwd_nb_max  = rbuf[5] / n_iter;
        double rev_blk_max = rbuf[6] / n_iter;
        double rev_nb_max  = rbuf[7] / n_iter;
        printf("\n--- Max (slowest rank) per-iteration (us) ---\n");
        printf("gatherp  blocking=%10.2f  nonblocking=%10.2f\n",
               fwd_blk_max * 1e6, fwd_nb_max * 1e6);
        printf("gathers  blocking=%10.2f  nonblocking=%10.2f\n",
               rev_blk_max * 1e6, rev_nb_max * 1e6);

        double tot_speedup = tot_blk / tot_nb;
        printf("\n=== Overall roundtrip speedup: %.2fx ===\n", tot_speedup);
        if (tot_speedup > 1.02)      printf("Nonblocking is BETTER\n");
        else if (tot_speedup < 0.98) printf("Blocking is BETTER\n");
        else                         printf("Performance is comparable\n");
    }

    MPI_Finalize();
    return 0;
}
