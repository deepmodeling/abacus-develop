#include "single_r_io.h"

#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <vector>

namespace
{
template <typename T>
struct SparseEntry
{
    int row;
    int col;
    T value;
};

template <typename T>
struct CsrBlock
{
    std::vector<T> values;
    std::vector<int> column_indices;
    std::vector<long long> row_ptr;
};

void write_data(std::ofstream& ofs, const double& data, const int precision)
{
    ofs << " " << std::fixed << std::scientific << std::setprecision(precision) << data;
}

void write_data(std::ofstream& ofs, const std::complex<double>& data, const int precision)
{
    ofs << " (" << std::fixed << std::scientific << std::setprecision(precision)
        << data.real() << "," << data.imag() << ")";
}

template <typename T>
std::vector<SparseEntry<T>> pack_local_entries(const ModuleIO::SparseRBlock<T>& matrix,
                                               const Parallel_Orbitals& pv,
                                               const ModuleIO::SparseWriteOptions& options,
                                               const int nlocal)
{
    std::vector<SparseEntry<T>> entries;
    for (typename ModuleIO::SparseRBlock<T>::const_iterator row = matrix.begin(); row != matrix.end(); ++row)
    {
        if (row->first >= static_cast<size_t>(nlocal))
        {
            ModuleBase::WARNING_QUIT("ModuleIO::output_single_R", "Sparse row index out of range.");
        }
        if (options.reduce && pv.global2local_row(row->first) < 0)
        {
            continue;
        }
        for (typename std::map<size_t, T>::const_iterator value = row->second.begin();
             value != row->second.end();
             ++value)
        {
            if (value->first >= static_cast<size_t>(nlocal))
            {
                ModuleBase::WARNING_QUIT("ModuleIO::output_single_R", "Sparse column index out of range.");
            }
            if (options.reduce || std::abs(value->second) > options.threshold)
            {
                SparseEntry<T> entry;
                entry.row = static_cast<int>(row->first);
                entry.col = static_cast<int>(value->first);
                entry.value = value->second;
                entries.push_back(entry);
            }
        }
    }
    return entries;
}

#ifdef __MPI
template <typename T>
std::vector<SparseEntry<T>> gather_entries(const std::vector<SparseEntry<T>>& local_entries,
                                           const MPI_Comm comm,
                                           const int root)
{
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (local_entries.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
    {
        ModuleBase::WARNING_QUIT("ModuleIO::output_single_R", "Too many local sparse entries for MPI_Gatherv.");
    }
    const int local_count = static_cast<int>(local_entries.size());
    std::vector<int> counts(rank == root ? size : 0);
    MPI_Gather(&local_count, 1, MPI_INT, rank == root ? counts.data() : nullptr, 1, MPI_INT, root, comm);

    std::vector<int> displacements(rank == root ? size : 0);
    int total_count = 0;
    if (rank == root)
    {
        for (int irank = 0; irank < size; ++irank)
        {
            displacements[irank] = total_count;
            if (counts[irank] > std::numeric_limits<int>::max() - total_count)
            {
                ModuleBase::WARNING_QUIT("ModuleIO::output_single_R", "Too many gathered sparse entries.");
            }
            total_count += counts[irank];
        }
    }

    std::vector<int> local_rows(local_count);
    std::vector<int> local_cols(local_count);
    std::vector<T> local_values(local_count);
    for (int index = 0; index < local_count; ++index)
    {
        local_rows[index] = local_entries[index].row;
        local_cols[index] = local_entries[index].col;
        local_values[index] = local_entries[index].value;
    }

    std::vector<int> rows(rank == root ? total_count : 0);
    std::vector<int> cols(rank == root ? total_count : 0);
    std::vector<T> values(rank == root ? total_count : 0);
    MPI_Gatherv(local_rows.data(), local_count, MPI_INT,
                rank == root ? rows.data() : nullptr,
                rank == root ? counts.data() : nullptr,
                rank == root ? displacements.data() : nullptr,
                MPI_INT, root, comm);
    MPI_Gatherv(local_cols.data(), local_count, MPI_INT,
                rank == root ? cols.data() : nullptr,
                rank == root ? counts.data() : nullptr,
                rank == root ? displacements.data() : nullptr,
                MPI_INT, root, comm);
    MPI_Gatherv(local_values.data(), local_count, Parallel_Reduce::MPI_Type<T>::value,
                rank == root ? values.data() : nullptr,
                rank == root ? counts.data() : nullptr,
                rank == root ? displacements.data() : nullptr,
                Parallel_Reduce::MPI_Type<T>::value, root, comm);

    std::vector<SparseEntry<T>> entries(rank == root ? total_count : 0);
    if (rank == root)
    {
        for (int index = 0; index < total_count; ++index)
        {
            entries[index].row = rows[index];
            entries[index].col = cols[index];
            entries[index].value = values[index];
        }
    }
    return entries;
}
#endif

template <typename T>
CsrBlock<T> make_csr(std::vector<SparseEntry<T>> entries, const int nlocal, const double threshold)
{
    std::stable_sort(entries.begin(), entries.end(), [](const SparseEntry<T>& lhs, const SparseEntry<T>& rhs) {
        return lhs.row < rhs.row || (lhs.row == rhs.row && lhs.col < rhs.col);
    });

    CsrBlock<T> csr;
    csr.values.reserve(entries.size());
    csr.column_indices.reserve(entries.size());
    csr.row_ptr.assign(nlocal + 1, 0);
    for (size_t begin = 0; begin < entries.size();)
    {
        size_t end = begin + 1;
        T value = entries[begin].value;
        while (end < entries.size()
               && entries[end].row == entries[begin].row
               && entries[end].col == entries[begin].col)
        {
            value += entries[end].value;
            ++end;
        }
        if (std::abs(value) > threshold)
        {
            csr.values.push_back(value);
            csr.column_indices.push_back(entries[begin].col);
            ++csr.row_ptr[entries[begin].row + 1];
        }
        begin = end;
    }
    for (int row = 1; row <= nlocal; ++row)
    {
        csr.row_ptr[row] += csr.row_ptr[row - 1];
    }
    return csr;
}

template <typename T>
void write_csr(std::ofstream& ofs, const CsrBlock<T>& csr, const ModuleIO::SparseWriteOptions& options)
{
    if (options.binary)
    {
        for (typename std::vector<T>::const_iterator value = csr.values.begin(); value != csr.values.end(); ++value)
        {
            ofs.write(reinterpret_cast<const char*>(&(*value)), sizeof(T));
        }
        for (std::vector<int>::const_iterator col = csr.column_indices.begin(); col != csr.column_indices.end(); ++col)
        {
            ofs.write(reinterpret_cast<const char*>(&(*col)), sizeof(int));
        }
        for (std::vector<long long>::const_iterator ptr = csr.row_ptr.begin(); ptr != csr.row_ptr.end(); ++ptr)
        {
            ofs.write(reinterpret_cast<const char*>(&(*ptr)), sizeof(long long));
        }
    }
    else
    {
        for (typename std::vector<T>::const_iterator value = csr.values.begin(); value != csr.values.end(); ++value)
        {
            write_data(ofs, *value, options.precision);
        }
        ofs << std::endl;
        for (std::vector<int>::const_iterator col = csr.column_indices.begin(); col != csr.column_indices.end(); ++col)
        {
            ofs << " " << *col;
        }
        ofs << std::endl;
        for (std::vector<long long>::const_iterator ptr = csr.row_ptr.begin(); ptr != csr.row_ptr.end(); ++ptr)
        {
            ofs << " " << *ptr;
        }
        ofs << std::endl;
    }
}
} // namespace

template <typename T>
void ModuleIO::output_single_R(std::ofstream& ofs,
                               const SparseRBlock<T>& matrix,
                               const Parallel_Orbitals& pv,
                               const SparseWriteOptions& options)
{
    const int nlocal = pv.get_global_row_size();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::output_single_R",
                                 "Parallel_Orbitals global row size must be positive.");
    }

    std::vector<SparseEntry<T>> entries = pack_local_entries(matrix, pv, options, nlocal);
    bool write_on_this_rank = true;
#ifdef __MPI
    if (options.reduce)
    {
        const MPI_Comm comm = pv.comm();
        if (comm != MPI_COMM_NULL)
        {
            int rank = 0;
            MPI_Comm_rank(comm, &rank);
            entries = gather_entries(entries, comm, 0);
            write_on_this_rank = rank == 0;
        }
    }
#endif
    if (write_on_this_rank)
    {
        write_csr(ofs, make_csr(entries, nlocal, options.threshold), options);
    }
}

template void ModuleIO::output_single_R<double>(std::ofstream&,
                                                 const SparseRBlock<double>&,
                                                 const Parallel_Orbitals&,
                                                 const SparseWriteOptions&);

template void ModuleIO::output_single_R<std::complex<double>>(std::ofstream&,
                                                               const SparseRBlock<std::complex<double>>&,
                                                               const Parallel_Orbitals&,
                                                               const SparseWriteOptions&);
