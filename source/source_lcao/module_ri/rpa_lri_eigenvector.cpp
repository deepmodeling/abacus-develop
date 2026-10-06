#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_eigen_vector_legacy(const Parallel_Orbitals& parav, const psi::Psi<T>& psi)
{
    const int nks_tot = this->runtime.input.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int npsin_tmp = this->runtime.input.nspin == 2 ? 2 : 1;
    const int nbands = parav.get_wfc_global_nbands();
    const int nbasis = parav.get_wfc_global_nbasis();
    const std::complex<double> zero(0.0, 0.0);
    const std::size_t values_per_iw = static_cast<std::size_t>(nbands) * npsin_tmp;
    // Preserve the legacy text writer when reader-v1 output is disabled.
#ifdef __MPI
    const MPI_Comm mpi_comm = parav.comm();
    const int mpi_rank = parav.get_coord_row() * parav.get_dim1() + parav.get_coord_col();
    const int mpi_size = parav.get_dim0() * parav.get_dim1();
    const auto check_mpi = [mpi_comm](const int local_error, const std::string& context) {
        const int local_failed = local_error == MPI_SUCCESS ? 0 : 1;
        int any_failed = 0;
        if (MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, mpi_comm) != MPI_SUCCESS || any_failed != 0)
        {
            throw std::runtime_error(context);
        }
    };
    std::vector<int> output_basis_counts(mpi_size, nbasis / mpi_size);
    for (int ip = 0; ip < nbasis % mpi_size; ++ip)
    {
        ++output_basis_counts[ip];
    }
    std::vector<int> output_basis_offsets(mpi_size + 1, 0);
    std::vector<int> output_basis_owner(nbasis);
    for (int ip = 0; ip < mpi_size; ++ip)
    {
        output_basis_offsets[ip + 1] = output_basis_offsets[ip] + output_basis_counts[ip];
        std::fill(output_basis_owner.begin() + output_basis_offsets[ip],
                  output_basis_owner.begin() + output_basis_offsets[ip + 1],
                  ip);
    }
    const int local_nw = output_basis_counts[mpi_rank];
#else
    const int mpi_rank = 0;
    const int local_nw = nbasis;
#endif
    const std::size_t local_size = static_cast<std::size_t>(local_nw) * values_per_iw;
    std::vector<std::complex<double>> local_wfc(local_size);
#ifdef __MPI
    struct WfcPackIndex
    {
        int spin;
        int band;
        int basis;
    };

    const int local_band_count = parav.ncol_bands;
    const unsigned long long send_size_wide
        = static_cast<unsigned long long>(local_band_count) * psi.get_nbasis() * npsin_tmp;
    const unsigned long long recv_size_wide = static_cast<unsigned long long>(local_nw) * values_per_iw;
    const unsigned long long max_alltoallv_count = std::numeric_limits<int>::max();
    check_mpi(send_size_wide <= max_alltoallv_count && recv_size_wide <= max_alltoallv_count ? MPI_SUCCESS
                                                                                             : MPI_ERR_COUNT,
              "RPA eigenvector buffer exceeds MPI_Alltoallv count range.");

    std::vector<int> send_counts(mpi_size, 0);
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (parav.global2local_col(ib) < 0)
        {
            continue;
        }
        for (int ir = 0; ir < psi.get_nbasis(); ++ir)
        {
            send_counts[output_basis_owner[parav.local2global_row(ir)]] += npsin_tmp;
        }
    }
    std::vector<int> send_displacements(mpi_size, 0);
    std::vector<int> recv_counts(mpi_size);
    std::vector<int> recv_displacements(mpi_size, 0);
    for (int ip = 1; ip < mpi_size; ++ip)
    {
        send_displacements[ip] = send_displacements[ip - 1] + send_counts[ip - 1];
    }
    check_mpi(MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, mpi_comm),
              "Failed to exchange RPA eigenvector redistribution counts.");
    for (int ip = 1; ip < mpi_size; ++ip)
    {
        recv_displacements[ip] = recv_displacements[ip - 1] + recv_counts[ip - 1];
    }
    const int send_size = static_cast<int>(send_size_wide);
    const int recv_size = static_cast<int>(recv_size_wide);
    check_mpi(recv_displacements.back() + recv_counts.back() == recv_size ? MPI_SUCCESS : MPI_ERR_COUNT,
              "RPA eigenvector Alltoallv redistribution map is inconsistent.");

    std::vector<WfcPackIndex> pack_indices(send_size);
    std::vector<unsigned long long> send_targets(send_size);
    std::vector<unsigned long long> recv_targets(recv_size);
    std::vector<int> send_positions = send_displacements;
    for (int ib = 0; ib < nbands; ++ib)
    {
        const int ib_local = parav.global2local_col(ib);
        if (ib_local < 0)
        {
            continue;
        }
        for (int ir = 0; ir < psi.get_nbasis(); ++ir)
        {
            const int iw = parav.local2global_row(ir);
            const int destination = output_basis_owner[iw];
            for (int is = 0; is < npsin_tmp; ++is)
            {
                const int position = send_positions[destination]++;
                pack_indices[position] = {is, ib_local, ir};
                send_targets[position]
                    = static_cast<unsigned long long>(iw - output_basis_offsets[destination]) * values_per_iw
                      + ib * npsin_tmp + is;
            }
        }
    }
    std::complex<double> dummy(0.0, 0.0);
    unsigned long long dummy_target = 0;
    check_mpi(MPI_Alltoallv(send_size > 0 ? send_targets.data() : &dummy_target,
                            send_counts.data(),
                            send_displacements.data(),
                            MPI_UNSIGNED_LONG_LONG,
                            recv_size > 0 ? recv_targets.data() : &dummy_target,
                            recv_counts.data(),
                            recv_displacements.data(),
                            MPI_UNSIGNED_LONG_LONG,
                            mpi_comm),
              "Failed to exchange RPA eigenvector destination indices.");
    std::vector<std::complex<double>> send_wfc(send_size);
    std::vector<std::complex<double>> recv_wfc(recv_size);
#endif
    std::string output_buffer;
    constexpr std::size_t output_line_bytes = 61;
    output_buffer.reserve(local_size * output_line_bytes + 32);

    const std::string filename = this->runtime.input.rpa_outdir + "KS_eigenvector.txt";
#ifdef __MPI
    MPI_File file = MPI_FILE_NULL;
    unsigned long long total_file_bytes = 0;
    const unsigned long long max_mpi_offset = static_cast<unsigned long long>(std::numeric_limits<MPI_Offset>::max());
    const char dummy_buffer = '\0';
    check_mpi(MPI_File_open(mpi_comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file),
              "Failed to open " + filename + ".");
    check_mpi(MPI_File_set_size(file, 0), "Failed to truncate " + filename + ".");
#else
    std::ofstream ofs(filename.c_str(), std::ios::out);
#endif

    for (int ik = 0; ik < nks_tot; ik++)
    {
#ifdef __MPI
        // 4.1 Pack and redistribute all bands and spins for this k point.
        for (int index = 0; index < send_size; ++index)
        {
            const WfcPackIndex& source = pack_indices[index];
            send_wfc[index] = psi(ik + nks_tot * source.spin, source.band, source.basis);
        }
        check_mpi(MPI_Alltoallv(send_size > 0 ? send_wfc.data() : &dummy,
                                send_counts.data(),
                                send_displacements.data(),
                                MPI_DOUBLE_COMPLEX,
                                recv_size > 0 ? recv_wfc.data() : &dummy,
                                recv_counts.data(),
                                recv_displacements.data(),
                                MPI_DOUBLE_COMPLEX,
                                mpi_comm),
                  "Failed to redistribute RPA eigenvectors with MPI_Alltoallv.");
        for (int index = 0; index < recv_size; ++index)
        {
            local_wfc[static_cast<std::size_t>(recv_targets[index])] = recv_wfc[index];
        }
#else
        for (int is = 0; is < npsin_tmp; is++)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int iw = 0; iw < psi.get_nbasis(); ++iw)
                {
                    local_wfc[iw * values_per_iw + ib * npsin_tmp + is] = psi(ik + nks_tot * is, ib, iw);
                }
            }
        }
#endif
        // 4.2 Format this rank's contiguous basis interval in basis-band-spin order.
        output_buffer.clear();
        if (mpi_rank == 0)
        {
            output_buffer += std::to_string(ik + 1);
            output_buffer += '\n';
        }
        char line_buffer[128];
        for (int iw = 0; iw < local_nw; ++iw)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int is = 0; is < npsin_tmp; is++)
                {
                    const std::complex<double>& value
                        = local_wfc[static_cast<std::size_t>(iw) * values_per_iw + ib * npsin_tmp + is];
                    const int line_length = std::snprintf(line_buffer,
                                                          sizeof(line_buffer),
                                                          "%30.15f%30.15f\n",
                                                          value.real(),
                                                          value.imag());
                    if (line_length != static_cast<int>(output_line_bytes))
                    {
                        throw std::runtime_error("Failed to format an RPA eigenvector value.");
                    }
                    output_buffer.append(line_buffer, static_cast<std::size_t>(line_length));
                }
            }
        }

#ifdef __MPI
        // 4.3 Fixed-width value lines make each rank's file offset deterministic;
        // only the k-point header length varies.
        const unsigned long long local_bytes = static_cast<unsigned long long>(output_buffer.size());
        const unsigned long long header_bytes = std::to_string(ik + 1).size() + 1;
        const unsigned long long bytes_per_basis = values_per_iw * output_line_bytes;
        if (bytes_per_basis > (max_mpi_offset - header_bytes) / std::max(nbasis, 1)
            || header_bytes + static_cast<unsigned long long>(nbasis) * bytes_per_basis
                   > max_mpi_offset - total_file_bytes)
        {
            throw std::runtime_error("RPA eigenvector file layout exceeds MPI_Offset range.");
        }
        const unsigned long long kpoint_bytes
            = header_bytes + static_cast<unsigned long long>(nbasis) * bytes_per_basis;
        const unsigned long long rank_offset_bytes
            = mpi_rank == 0
                  ? 0
                  : header_bytes + static_cast<unsigned long long>(output_basis_offsets[mpi_rank]) * bytes_per_basis;
        const unsigned long long file_offset_bytes = total_file_bytes + rank_offset_bytes;
        // Rank 0 has the largest basis slice and additionally owns the k-point header.
        const unsigned long long max_local_bytes
            = header_bytes + static_cast<unsigned long long>(output_basis_counts[0]) * bytes_per_basis;
        const unsigned long long max_write_count = static_cast<unsigned long long>(std::numeric_limits<int>::max());
        // MPI_File_write_at_all is collective, so every rank executes the same
        // number of chunks; ranks without data in a chunk pass a zero count.
        for (unsigned long long buffer_offset = 0; buffer_offset < max_local_bytes; buffer_offset += max_write_count)
        {
            const unsigned long long remaining = local_bytes > buffer_offset ? local_bytes - buffer_offset : 0;
            const int write_count = static_cast<int>(std::min(remaining, max_write_count));
            const char* write_buffer
                = write_count > 0 ? output_buffer.data() + static_cast<std::size_t>(buffer_offset) : &dummy_buffer;
            const unsigned long long write_offset
                = write_count > 0 ? file_offset_bytes + buffer_offset : file_offset_bytes;
            check_mpi(MPI_File_write_at_all(file,
                                            static_cast<MPI_Offset>(write_offset),
                                            write_buffer,
                                            write_count,
                                            MPI_CHAR,
                                            MPI_STATUS_IGNORE),
                      "Failed to write " + filename + ".");
        }
        total_file_bytes += kpoint_bytes;
#else
        ofs << output_buffer;
#endif
    } // ik
#ifdef __MPI
    check_mpi(MPI_File_sync(file), "Failed to sync " + filename + ".");
    check_mpi(MPI_File_close(&file), "Failed to close " + filename + ".");
#endif
}

template void RPA_LRI<double, double>::out_eigen_vector_legacy(const Parallel_Orbitals& parav,
                                                               const psi::Psi<double>& psi);
template void RPA_LRI<std::complex<double>, double>::out_eigen_vector_legacy(const Parallel_Orbitals& parav,
                                                                             const psi::Psi<std::complex<double>>& psi);
