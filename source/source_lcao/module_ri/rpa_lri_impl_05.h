            RpaLriDetail::write_ks_eigenvector_v1_mpi(io_comm,
                                                      parav,
                                                      psi,
                                                      nks_tot,
                                                      npsin_tmp,
                                                      this->runtime.input.nspin,
                                                      this->runtime.input.nbands,
                                                      this->runtime.nlocal,
                                                      this->runtime.input.rpa_outdir + "KS_wfc_0.dat");
        }
        catch (...)
        {
            ModuleBase::timer::end("RPA_LRI", "out_eigen_vector_v1_mpi_io");
            throw;
        }
        ModuleBase::timer::end("RPA_LRI", "out_eigen_vector_v1_mpi_io");
#else
        struct KSEigenRecord
        {
            std::int32_t ik = 0;
            std::int64_t payload_offset = 0;
            std::vector<std::complex<double>> payload;
        };

        std::vector<KSEigenRecord> records;
        records.reserve(static_cast<std::size_t>(nks_tot));

        for (int ik = 0; ik < nks_tot; ik++)
        {
            std::vector<ModuleBase::ComplexMatrix> is_wfc_ib_iw(npsin_tmp);
            for (int is = 0; is < npsin_tmp; is++)
            {
                is_wfc_ib_iw[is].create(this->runtime.input.nbands, this->runtime.nlocal);
                for (int ib_global = 0; ib_global < this->runtime.input.nbands; ++ib_global)
                {
                    std::vector<std::complex<double>> wfc_iks(this->runtime.nlocal, zero);

                    const int ib_local = parav.global2local_col(ib_global);

                    if (ib_local >= 0)
                    {
                        for (int ir = 0; ir < psi.get_nbasis(); ir++)
                        {
                            wfc_iks[parav.local2global_row(ir)] = psi(ik + nks_tot * is, ib_local, ir);
                        }
                    }

                    for (int iw = 0; iw < this->runtime.nlocal; iw++)
                    {
                        is_wfc_ib_iw[is](ib_global, iw) = wfc_iks[iw];
                    }
                }
            }

            if (this->runtime.rank == 0)
            {
                KSEigenRecord record;
                record.ik = static_cast<std::int32_t>(ik + 1);
                record.payload.reserve(static_cast<std::size_t>(npsin_tmp)
                                       * static_cast<std::size_t>(this->runtime.input.nbands)
                                       * static_cast<std::size_t>(this->runtime.nlocal));
                if (this->runtime.input.nspin == 4)
                {
                    if (this->runtime.nlocal % 2 != 0)
                    {
                        throw std::runtime_error("SOC KS eigenvector output expects an even basis size.");
                    }
                    const int nlocal_ao = this->runtime.nlocal / 2;
                    for (int isoc = 0; isoc < 2; ++isoc)
                    {
                        for (int ib = 0; ib < this->runtime.input.nbands; ++ib)
                        {
                            for (int iw = 0; iw < nlocal_ao; ++iw)
                            {
                                record.payload.push_back(is_wfc_ib_iw[0](ib, iw * 2 + isoc));
                            }
                        }
                    }
                }
                else
                {
                    for (int is = 0; is < npsin_tmp; ++is)
                    {
                        for (int ib = 0; ib < this->runtime.input.nbands; ++ib)
                        {
                            for (int iw = 0; iw < this->runtime.nlocal; ++iw)
                            {
                                record.payload.push_back(is_wfc_ib_iw[is](ib, iw));
                            }
                        }
                    }
                }
                records.push_back(std::move(record));
            }
        }

        if (this->runtime.rank == 0)
        {
            const std::string out_name = this->runtime.input.rpa_outdir + "KS_wfc_0.dat";
            const std::int64_t record_bytes = static_cast<std::int64_t>(sizeof(std::int32_t))
                + static_cast<std::int64_t>(sizeof(std::int64_t));
            std::int64_t offset = 6 * static_cast<std::int64_t>(sizeof(std::int32_t))
                + static_cast<std::int64_t>(records.size()) * record_bytes;
            for (auto& record: records)
            {
                record.payload_offset = offset;
                offset += static_cast<std::int64_t>(record.payload.size() * sizeof(std::complex<double>));
            }

            std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            if (!ofs.good())
            {
                throw std::runtime_error("Failed to open " + out_name);
            }

            const std::int32_t marker = RpaLriDetail::LIBRPA_KS_EIGENVECTOR_V1_MARKER;
            const std::int32_t kind = RpaLriDetail::LIBRPA_KS_EIGENVECTOR_V1_KIND_COMPLEX_DOUBLE;
            const std::int32_t nkpoints_local = RpaLriDetail::checked_i32_from_size(records.size(),
                                                                                    "KS eigenvector k-point count");
            const std::int32_t nspins = RpaLriDetail::checked_i32_from_int(npsin_tmp, "KS eigenvector spin count");
            const std::int32_t nstates = RpaLriDetail::checked_i32_from_int(this->runtime.input.nbands,
                                                                            "KS eigenvector state count");
            const std::int32_t nbasis_wfc = RpaLriDetail::checked_i32_from_int(this->runtime.nlocal,
                                                                               "KS eigenvector basis count");
            RpaLriDetail::write_scalar(ofs, marker, out_name);
            RpaLriDetail::write_scalar(ofs, kind, out_name);
            RpaLriDetail::write_scalar(ofs, nkpoints_local, out_name);
            RpaLriDetail::write_scalar(ofs, nspins, out_name);
            RpaLriDetail::write_scalar(ofs, nstates, out_name);
            RpaLriDetail::write_scalar(ofs, nbasis_wfc, out_name);
            for (const auto& record: records)
            {
                RpaLriDetail::write_scalar(ofs, record.ik, out_name);
                RpaLriDetail::write_scalar(ofs, record.payload_offset, out_name);
            }
            for (const auto& record: records)
            {
                RpaLriDetail::checked_write(ofs,
                                            record.payload.data(),
                                            record.payload.size() * sizeof(std::complex<double>),
                                            out_name);
            }
            ofs.close();
        }
#endif
        return;
    }

    // Preserve the legacy text writer when reader-v1 output is disabled.
#ifdef __MPI
    const MPI_Comm mpi_comm = parav.comm();
    const int mpi_rank = parav.get_coord_row() * parav.get_dim1() + parav.get_coord_col();
    const int mpi_size = parav.get_dim0() * parav.get_dim1();
    const auto check_mpi = [mpi_comm](const int local_error, const std::string& context) {
        const int local_failed = local_error == MPI_SUCCESS ? 0 : 1;
        int any_failed = 0;
        if (MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, mpi_comm) != MPI_SUCCESS
            || any_failed != 0)
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
    const unsigned long long recv_size_wide
        = static_cast<unsigned long long>(local_nw) * values_per_iw;
    const unsigned long long max_alltoallv_count = std::numeric_limits<int>::max();
    check_mpi(send_size_wide <= max_alltoallv_count && recv_size_wide <= max_alltoallv_count
                  ? MPI_SUCCESS
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
    check_mpi(recv_displacements.back() + recv_counts.back() == recv_size
                  ? MPI_SUCCESS
                  : MPI_ERR_COUNT,
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
                    = static_cast<unsigned long long>(iw - output_basis_offsets[destination])
                          * values_per_iw
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
    const unsigned long long max_mpi_offset
        = static_cast<unsigned long long>(std::numeric_limits<MPI_Offset>::max());
    const char dummy_buffer = '\0';
    check_mpi(MPI_File_open(mpi_comm,
                            filename.c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY,
                            MPI_INFO_NULL,
                            &file),
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
                    local_wfc[iw * values_per_iw + ib * npsin_tmp + is]
                        = psi(ik + nks_tot * is, ib, iw);
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
        const unsigned long long local_bytes
            = static_cast<unsigned long long>(output_buffer.size());
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
                  : header_bytes
                        + static_cast<unsigned long long>(output_basis_offsets[mpi_rank])
                              * bytes_per_basis;
        const unsigned long long file_offset_bytes = total_file_bytes + rank_offset_bytes;
        // Rank 0 has the largest basis slice and additionally owns the k-point header.
        const unsigned long long max_local_bytes
            = header_bytes
              + static_cast<unsigned long long>(output_basis_counts[0]) * bytes_per_basis;
        const unsigned long long max_write_count
            = static_cast<unsigned long long>(std::numeric_limits<int>::max());
        // MPI_File_write_at_all is collective, so every rank executes the same
        // number of chunks; ranks without data in a chunk pass a zero count.
        for (unsigned long long buffer_offset = 0;
             buffer_offset < max_local_bytes;
             buffer_offset += max_write_count)
        {
            const unsigned long long remaining
                = local_bytes > buffer_offset ? local_bytes - buffer_offset : 0;
            const int write_count
                = static_cast<int>(std::min(remaining, max_write_count));
            const char* write_buffer
                = write_count > 0
                      ? output_buffer.data() + static_cast<std::size_t>(buffer_offset)
                      : &dummy_buffer;
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

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_struc(const UnitCell& ucell)
{
    if (this->runtime.rank != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_struc");
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    const ModuleBase::Matrix3 lat = ucell.latvec * ucell.lat0; // in unit of Bohr
    const ModuleBase::Matrix3 G_RPA = ucell.G * (ModuleBase::TWO_PI / ucell.lat0); // in unit of 1/Bohr
    std::ofstream ofs;
    ofs.open(this->runtime.input.rpa_outdir + "stru_out.txt", std::ios::out);
    ofs << std::fixed << std::setprecision(9);
    ofs << lat.e11 << std::setw(15) << lat.e12 << std::setw(15) << lat.e13 << std::endl;
    ofs << lat.e21 << std::setw(15) << lat.e22 << std::setw(15) << lat.e23 << std::endl;
    ofs << lat.e31 << std::setw(15) << lat.e32 << std::setw(15) << lat.e33 << std::endl;

    ofs << G_RPA.e11 << std::setw(15) << G_RPA.e12 << std::setw(15) << G_RPA.e13 << std::endl;
    ofs << G_RPA.e21 << std::setw(15) << G_RPA.e22 << std::setw(15) << G_RPA.e23 << std::endl;
    ofs << G_RPA.e31 << std::setw(15) << G_RPA.e32 << std::setw(15) << G_RPA.e33 << std::endl;

    ofs << ucell.nat << std::endl;
    for (int it = 0; it < ucell.ntype; it++)
