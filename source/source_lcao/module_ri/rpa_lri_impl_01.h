    const unsigned long long payload_bytes
        = checked_mul_u64(static_cast<unsigned long long>(nks_tot),
                          kpoint_bytes,
                          "KS eigenvector v1 payload byte size");
    metadata.total_bytes
        = checked_add_u64(payload_begin, payload_bytes, "KS eigenvector v1 file size");

    const std::int32_t marker = LIBRPA_KS_EIGENVECTOR_V1_MARKER;
    const std::int32_t kind = LIBRPA_KS_EIGENVECTOR_V1_KIND_COMPLEX_DOUBLE;
    const std::int32_t nkpoints_i32 = checked_i32_from_int(nks_tot, "KS eigenvector k-point count");
    const std::int32_t nspins_i32 = checked_i32_from_int(nspins, "KS eigenvector spin count");
    const std::int32_t nstates_i32 = checked_i32_from_int(nbands, "KS eigenvector state count");
    const std::int32_t nbasis_i32 = checked_i32_from_int(nbasis_wfc, "KS eigenvector basis count");
    append_binary_scalar(metadata.bytes, marker);
    append_binary_scalar(metadata.bytes, kind);
    append_binary_scalar(metadata.bytes, nkpoints_i32);
    append_binary_scalar(metadata.bytes, nspins_i32);
    append_binary_scalar(metadata.bytes, nstates_i32);
    append_binary_scalar(metadata.bytes, nbasis_i32);

    metadata.payload_offsets.reserve(static_cast<std::size_t>(nks_tot));
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        const unsigned long long kpoint_offset
            = checked_mul_u64(static_cast<unsigned long long>(ik),
                              kpoint_bytes,
                              "KS eigenvector v1 k-point offset");
        const unsigned long long payload_offset
            = checked_add_u64(payload_begin, kpoint_offset, "KS eigenvector v1 payload offset");
        metadata.payload_offsets.push_back(payload_offset);
        const std::int32_t ik_file = checked_i32_from_int(ik + 1, "KS eigenvector k-point index");
        const std::int64_t payload_offset_i64
            = checked_i64_from_u64(payload_offset, "KS eigenvector v1 payload offset");
        append_binary_scalar(metadata.bytes, ik_file);
        append_binary_scalar(metadata.bytes, payload_offset_i64);
    }
    if (metadata.bytes.size() != static_cast<std::size_t>(payload_begin))
    {
        throw std::runtime_error("KS eigenvector v1 metadata size is inconsistent.");
    }
    return metadata;
}

struct KSEigenvectorPackCursor
{
    int band_local = 0;
    int row_local = 0;
};

template <typename T>
inline int pack_ks_eigenvector_chunk(const Parallel_Orbitals& parav,
                                     const psi::Psi<T>& psi,
                                     const int psi_k,
                                     const bool is_soc,
                                     const int spinor_component,
                                     KSEigenvectorPackCursor& cursor,
                                     std::vector<std::complex<double>>& buffer,
                                     const int requested_count)
{
    int packed = 0;
    while (cursor.band_local < parav.ncol_bands && packed < requested_count)
    {
        while (cursor.row_local < psi.get_nbasis() && packed < requested_count)
        {
            const int ir_local = cursor.row_local++;
            const int iw_global = parav.local2global_row(ir_local);
            if (is_soc && iw_global % 2 != spinor_component)
            {
                continue;
            }
            buffer[static_cast<std::size_t>(packed)]
                = std::complex<double>(psi(psi_k, cursor.band_local, ir_local));
            ++packed;
        }
        if (cursor.row_local == psi.get_nbasis())
        {
            ++cursor.band_local;
            cursor.row_local = 0;
        }
    }
    return packed;
}

template <typename T>
inline void write_ks_eigenvector_v1_mpi(const MPI_Comm mpi_comm,
                                        const Parallel_Orbitals& parav,
                                        const psi::Psi<T>& psi,
                                        const int nks_tot,
                                        const int nspins,
                                        const int nspin_abacus,
                                        const int nbands,
                                        const int nbasis_wfc,
                                        const std::string& final_name)
{
    if (mpi_comm == MPI_COMM_NULL)
    {
        throw std::runtime_error("KS eigenvector MPI-IO writer received MPI_COMM_NULL.");
    }
    int mpi_rank = 0;
    collective_mpi_check(mpi_comm,
                         MPI_Comm_rank(mpi_comm, &mpi_rank),
                         "Failed to query the KS eigenvector MPI rank");

    const bool is_soc = nspin_abacus == 4;
    const int ncomponents = is_soc ? 2 : nspins;
    bool dimensions_valid
        = nks_tot >= 0 && nbands >= 0 && nbasis_wfc >= 0 && nspins > 0
          && (nspin_abacus == 1 || nspin_abacus == 2 || nspin_abacus == 4)
          && nspins == (nspin_abacus == 2 ? 2 : 1)
          && (!is_soc || nbasis_wfc % 2 == 0)
          && psi.get_nbasis() == parav.get_row_size()
          && psi.get_nbands() == parav.ncol_bands
          && parav.get_wfc_global_nbasis() == nbasis_wfc
          && parav.get_wfc_global_nbands() == nbands;
    if (dimensions_valid)
    {
        const unsigned long long required_psi_kpoints
            = checked_mul_u64(static_cast<unsigned long long>(nks_tot),
                              static_cast<unsigned long long>(nspins),
                              "KS eigenvector MPI source k-point count");
        dimensions_valid
            = required_psi_kpoints <= static_cast<unsigned long long>(psi.get_nk());
    }
    collective_require(mpi_comm, dimensions_valid, "Invalid KS eigenvector MPI-IO dimensions");

    const int spatial_basis = is_soc ? nbasis_wfc / 2 : nbasis_wfc;
    std::vector<KSEigenvectorMpiLayout> layouts;
    layouts.reserve(static_cast<std::size_t>(is_soc ? 2 : 1));
    for (int component = 0; component < (is_soc ? 2 : 1); ++component)
    {
        layouts.push_back(make_ks_eigenvector_mpi_layout(
            mpi_comm, parav, nbands, nbasis_wfc, is_soc, component));
    }

    const unsigned long long chunk_elements = 4ULL * 1024ULL * 1024ULL;
    unsigned long long local_buffer_elements = 0;
    for (const auto& layout: layouts)
    {
        local_buffer_elements = std::max(local_buffer_elements,
                                         std::min(layout.local_count, chunk_elements));
    }
    bool buffer_allocated = true;
    std::vector<std::complex<double>> buffer;
    try
    {
        buffer.resize(static_cast<std::size_t>(local_buffer_elements));
    }
    catch (const std::exception&)
    {
        buffer_allocated = false;
    }
    collective_require(mpi_comm,
                       buffer_allocated,
                       "Failed to allocate the bounded KS eigenvector MPI pack buffer");
    std::complex<double> dummy(0.0, 0.0);

    const KSEigenvectorV1Metadata metadata
        = make_ks_eigenvector_v1_metadata(
            nks_tot, nspins, ncomponents, nbands, nbasis_wfc, spatial_basis);
    collective_require(
        mpi_comm,
        metadata.bytes.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "KS eigenvector MPI metadata exceeds the MPI count range");

    const std::string temporary_name = final_name + ".tmp";
    MPI_File file = MPI_FILE_NULL;
    bool file_opened = false;
    try
    {
        int mpi_error = MPI_File_open(mpi_comm,
                                      const_cast<char*>(temporary_name.c_str()),
                                      MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                      MPI_INFO_NULL,
                                      &file);
        collective_mpi_check(mpi_comm, mpi_error, "Failed to open " + temporary_name);
        file_opened = true;
        collective_mpi_check(mpi_comm,
                             MPI_File_set_errhandler(file, MPI_ERRORS_RETURN),
                             "Failed to set the KS eigenvector MPI file error handler");
        collective_mpi_check(
            mpi_comm,
            MPI_File_set_size(
                file,
                checked_mpi_offset_from_u64(metadata.total_bytes, "KS eigenvector MPI file size")),
            "Failed to set the KS eigenvector MPI file size");

        const int metadata_count
            = mpi_rank == 0 ? static_cast<int>(metadata.bytes.size()) : 0;
        MPI_Status metadata_status;
        mpi_error = MPI_File_write_at_all(file,
                                          0,
                                          metadata_count == 0
                                              ? static_cast<void*>(&dummy)
                                              : static_cast<void*>(
                                                  const_cast<char*>(metadata.bytes.data())),
                                          metadata_count,
                                          MPI_BYTE,
                                          &metadata_status);
        collective_mpi_check(mpi_comm, mpi_error, "Failed to write KS eigenvector MPI metadata");

        for (int ik = 0; ik < nks_tot; ++ik)
        {
            for (int component = 0; component < ncomponents; ++component)
            {
                KSEigenvectorMpiLayout& layout = layouts[is_soc ? component : 0];
                const unsigned long long component_offset
                    = checked_mul_u64(static_cast<unsigned long long>(component),
                                      metadata.component_bytes,
                                      "KS eigenvector MPI component offset");
                const unsigned long long displacement
                    = checked_add_u64(metadata.payload_offsets[static_cast<std::size_t>(ik)],
                                      component_offset,
                                      "KS eigenvector MPI view displacement");
                collective_mpi_check(
                    mpi_comm,
                    MPI_File_set_view(
                        file,
                        checked_mpi_offset_from_u64(displacement,
                                                    "KS eigenvector MPI view displacement"),
                        MPI_C_DOUBLE_COMPLEX,
                        layout.filetype,
                        const_cast<char*>("native"),
                        MPI_INFO_NULL),
                    "Failed to set the KS eigenvector MPI file view");

                KSEigenvectorPackCursor cursor;
                const int psi_k = is_soc ? ik : ik + nks_tot * component;
                unsigned long long chunk_begin = 0;
                while (chunk_begin < layout.max_local_count)
                {
                    const unsigned long long local_remaining
                        = chunk_begin < layout.local_count ? layout.local_count - chunk_begin : 0;
                    const int local_chunk_count = static_cast<int>(
                        std::min(local_remaining, chunk_elements));
                    const int packed = pack_ks_eigenvector_chunk(parav,
                                                                 psi,
                                                                 psi_k,
                                                                 is_soc,
                                                                 component,
                                                                 cursor,
                                                                 buffer,
                                                                 local_chunk_count);
                    collective_require(mpi_comm,
                                       packed == local_chunk_count,
                                       "Failed to pack the expected KS eigenvector MPI chunk");

                    MPI_Status payload_status;
                    mpi_error = MPI_File_write_at_all(
                        file,
                        checked_mpi_offset_from_u64(chunk_begin,
                                                    "KS eigenvector MPI chunk offset"),
                        local_chunk_count == 0 ? static_cast<void*>(&dummy)
                                               : static_cast<void*>(buffer.data()),
                        local_chunk_count,
                        MPI_C_DOUBLE_COMPLEX,
                        &payload_status);
                    collective_mpi_check(mpi_comm,
                                         mpi_error,
                                         "Failed to write a KS eigenvector MPI payload chunk");
                    if (layout.max_local_count - chunk_begin <= chunk_elements)
                    {
                        break;
                    }
                    chunk_begin = checked_add_u64(
                        chunk_begin, chunk_elements, "KS eigenvector MPI chunk offset");
                }
            }
        }

        mpi_error = MPI_File_sync(file);
        collective_mpi_check(mpi_comm, mpi_error, "Failed to sync the KS eigenvector MPI file");
        mpi_error = MPI_File_close(&file);
        file_opened = false;
        collective_mpi_check(mpi_comm, mpi_error, "Failed to close the KS eigenvector MPI file");

        int free_error = MPI_SUCCESS;
        for (auto& layout: layouts)
        {
            const int layout_error = free_ks_eigenvector_mpi_layout(layout);
            if (free_error == MPI_SUCCESS && layout_error != MPI_SUCCESS)
            {
                free_error = layout_error;
            }
        }
        collective_mpi_check(mpi_comm,
                             free_error,
                             "Failed to free a KS eigenvector MPI file type");

        const int rename_error
            = mpi_rank == 0 && std::rename(temporary_name.c_str(), final_name.c_str()) != 0
                  ? MPI_ERR_IO
                  : MPI_SUCCESS;
        collective_mpi_check(mpi_comm,
                             rename_error,
                             "Failed to publish the completed KS eigenvector MPI file");
    }
    catch (...)
    {
        if (file_opened)
        {
            MPI_File_close(&file);
        }
        for (auto& layout: layouts)
        {
            free_ks_eigenvector_mpi_layout(layout);
        }
        MPI_Barrier(mpi_comm);
        if (mpi_rank == 0)
        {
            std::remove(temporary_name.c_str());
        }
        MPI_Barrier(mpi_comm);
        throw;
    }

    if (mpi_rank == 0)
    {
        std::cout << "KS eigenvector writer: binary v1 MPI-IO, bounded pack buffer, file "
                  << final_name << std::endl;
    }
}
#endif

inline int sum_int_vector(const std::vector<int>& values)
{
    int sum = 0;
    for (const int value: values)
    {
        if (value > std::numeric_limits<int>::max() - sum)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
        }
        sum += value;
    }
    return sum;
}

template <typename Value>
inline double real_as_double(const Value& value)
{
    return static_cast<double>(value);
}

inline double real_as_double(const std::complex<double>& value)
{
    return value.real();
}

inline std::vector<std::vector<int>>
collect_abfs_l_nchi(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    std::vector<std::vector<int>> abfs_l_nchi;
    abfs_l_nchi.reserve(abfs.size());
    for (const auto& abfs_type : abfs)
    {
        std::vector<int> shell_counts;
        shell_counts.reserve(abfs_type.size());
        for (const auto& abfs_l : abfs_type)
        {
            shell_counts.push_back(static_cast<int>(abfs_l.size()));
        }
        abfs_l_nchi.push_back(std::move(shell_counts));
    }
    return abfs_l_nchi;
}

inline std::vector<std::vector<int>>
collect_wfc_l_nchi(const UnitCell& ucell)
{
    std::vector<std::vector<int>> wfc_l_nchi;
    wfc_l_nchi.reserve(static_cast<std::size_t>(ucell.ntype));
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const auto& atom = ucell.atoms[itype];
        if (atom.nwl < 0 || atom.l_nchi.size() < static_cast<std::size_t>(atom.nwl + 1))
        {
            throw std::runtime_error("LibRPA v1 basis output found inconsistent AO shell counts.");
        }
        std::vector<int> shell_counts;
        shell_counts.reserve(static_cast<std::size_t>(atom.nwl + 1));
        for (int l = 0; l <= atom.nwl; ++l)
        {
            if (atom.l_nchi[static_cast<std::size_t>(l)] < 0)
            {
                throw std::runtime_error("LibRPA v1 basis output found negative AO shell count.");
            }
            shell_counts.push_back(atom.l_nchi[static_cast<std::size_t>(l)]);
        }
        wfc_l_nchi.push_back(std::move(shell_counts));
    }
    return wfc_l_nchi;
}

inline int basis_size_from_shell_counts(const std::vector<int>& shell_counts)
{
    int basis_size = 0;
    for (std::size_t l = 0; l < shell_counts.size(); ++l)
    {
        const int count = shell_counts[l];
        if (count < 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found negative shell count.");
        }
        const int shell_size = static_cast<int>(2 * l + 1);
        if (count > 0 && shell_size > std::numeric_limits<int>::max() / count)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 shell sizes.");
        }
        const int increment = count * shell_size;
        if (increment > std::numeric_limits<int>::max() - basis_size)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 shell sizes.");
        }
        basis_size += increment;
    }
    return basis_size;
}

inline void write_librpa_split_basis_file(const UnitCell& ucell,
                                          const std::vector<int>& type_sizes,
                                          const std::vector<std::vector<int>>& type_shell_counts,
                                          const std::string& filename)
{
    if (type_sizes.size() != static_cast<std::size_t>(ucell.ntype)
        || type_shell_counts.size() != static_cast<std::size_t>(ucell.ntype))
    {
        throw std::runtime_error("LibRPA v1 basis output found inconsistent atom-type count.");
    }

    int total_basis = 0;
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const int type_size = type_sizes[static_cast<std::size_t>(itype)];
        if (type_size <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type basis size.");
        }
        const int shell_size = basis_size_from_shell_counts(type_shell_counts[static_cast<std::size_t>(itype)]);
        if (type_size != shell_size)
        {
            throw std::runtime_error("LibRPA v1 basis output found inconsistent shell layout size.");
        }
        if (ucell.atoms[itype].na < 0 || type_size > std::numeric_limits<int>::max() / ucell.atoms[itype].na)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
        }
        const int increment = type_size * ucell.atoms[itype].na;
        if (increment > std::numeric_limits<int>::max() - total_basis)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
