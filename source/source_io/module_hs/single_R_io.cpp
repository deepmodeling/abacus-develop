#include "single_R_io.h"
#include "source_base/parallel_reduce.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include <vector>
#include <algorithm>
#include <complex>
#include <type_traits>
#include <iomanip>

inline void write_data(std::ofstream& ofs, const double& data)
{
    ofs << " " << data;
}
inline void write_data(std::ofstream& ofs, const std::complex<double>& data)
{
    ofs << " (" << data.real() << "," << data.imag() << ")";
}

template<typename T>
void ModuleIO::output_single_R(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, T>>& XR,
    const double& sparse_threshold,
    const bool& binary,
    const Parallel_Orbitals& pv,
    const bool& reduce)
{
    // Step 1: Collect non-zero elements local to this processor from XR
    std::vector<int> local_rows;
    std::vector<int> local_cols;
    std::vector<T> local_vals;

    for (auto const& row_pair : XR)
    {
        auto const& row = row_pair.first;
        for (auto const& col_pair : row_pair.second)
        {
            auto const& col = col_pair.first;
            auto const& val = col_pair.second;
            if (std::abs(val) > sparse_threshold)
            {
                local_rows.push_back(row);
                local_cols.push_back(col);
                local_vals.push_back(val);
            }
        }
    }

    // Step 2: Handle gathering
#ifdef __MPI
    if (reduce && GlobalV::NPROC > 1)
    {
        int local_size = local_rows.size();
        std::vector<int> recvcounts;
        std::vector<int> displs;
        int total_size = 0;

        if (GlobalV::DRANK == 0)
        {
            recvcounts.resize(GlobalV::NPROC);
        }

        MPI_Gather(&local_size, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<int> global_rows;
        std::vector<int> global_cols;
        std::vector<T> global_vals;

        if (GlobalV::DRANK == 0)
        {
            displs.resize(GlobalV::NPROC);
            displs[0] = 0;
            for (int i = 0; i < GlobalV::NPROC; ++i)
            {
                total_size += recvcounts[i];
                if (i > 0)
                {
                    displs[i] = displs[i - 1] + recvcounts[i - 1];
                }
            }
            global_rows.resize(total_size);
            global_cols.resize(total_size);
            global_vals.resize(total_size);
        }

        MPI_Datatype mpi_type;
        if (std::is_same<T, double>::value)
        {
            mpi_type = MPI_DOUBLE;
        }
        else
        {
            mpi_type = MPI_DOUBLE_COMPLEX;
        }

        MPI_Gatherv(local_rows.data(), local_size, MPI_INT,
                    global_rows.data(), recvcounts.data(), displs.data(), MPI_INT,
                    0, MPI_COMM_WORLD);
        MPI_Gatherv(local_cols.data(), local_size, MPI_INT,
                    global_cols.data(), recvcounts.data(), displs.data(), MPI_INT,
                    0, MPI_COMM_WORLD);
        MPI_Gatherv(local_vals.data(), local_size, mpi_type,
                    global_vals.data(), recvcounts.data(), displs.data(), mpi_type,
                    0, MPI_COMM_WORLD);

        if (GlobalV::DRANK == 0)
        {
            std::vector<std::vector<std::pair<int, T>>> gathered_rows(PARAM.globalv.nlocal);
            for (int i = 0; i < total_size; ++i)
            {
                gathered_rows[global_rows[i]].push_back({global_cols[i], global_vals[i]});
            }

            for (int r = 0; r < PARAM.globalv.nlocal; ++r)
            {
                std::sort(gathered_rows[r].begin(), gathered_rows[r].end(),
                          [](const std::pair<int, T>& a, const std::pair<int, T>& b) { return a.first < b.first; });
            }

            std::vector<int> col_indices;
            std::vector<T> values;
            std::vector<long long> indptr(PARAM.globalv.nlocal + 1);
            indptr[0] = 0;

            for (int r = 0; r < PARAM.globalv.nlocal; ++r)
            {
                for (auto const& col_pair : gathered_rows[r])
                {
                    col_indices.push_back(col_pair.first);
                    values.push_back(col_pair.second);
                }
                indptr[r + 1] = col_indices.size();
            }

            // Write to file
            if (binary)
            {
                if (!values.empty())
                {
                    ofs.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(T));
                    ofs.write(reinterpret_cast<const char*>(col_indices.data()), col_indices.size() * sizeof(int));
                }
                ofs.write(reinterpret_cast<const char*>(indptr.data()), indptr.size() * sizeof(long long));
            }
            else
            {
                ofs << std::fixed << std::scientific << std::setprecision(8);
                for (auto const& val : values)
                {
                    write_data(ofs, val);
                }
                ofs << std::endl;

                for (auto const& col : col_indices)
                {
                    ofs << " " << col;
                }
                ofs << std::endl;

                for (auto const& i : indptr)
                {
                    ofs << " " << i;
                }
                ofs << std::endl;
            }
        }
        return;
    }
#endif

    // Serial mode or no reduce
    if (!reduce || GlobalV::DRANK == 0)
    {
        // Reconstruct local rows and construct CSR directly
        std::vector<std::vector<std::pair<int, T>>> gathered_rows(PARAM.globalv.nlocal);
        for (size_t i = 0; i < local_rows.size(); ++i)
        {
            gathered_rows[local_rows[i]].push_back({local_cols[i], local_vals[i]});
        }

        for (int r = 0; r < PARAM.globalv.nlocal; ++r)
        {
            std::sort(gathered_rows[r].begin(), gathered_rows[r].end(),
                      [](const std::pair<int, T>& a, const std::pair<int, T>& b) { return a.first < b.first; });
        }

        std::vector<int> col_indices;
        std::vector<T> values;
        std::vector<long long> indptr(PARAM.globalv.nlocal + 1);
        indptr[0] = 0;

        for (int r = 0; r < PARAM.globalv.nlocal; ++r)
        {
            for (auto const& col_pair : gathered_rows[r])
            {
                col_indices.push_back(col_pair.first);
                values.push_back(col_pair.second);
            }
            indptr[r + 1] = col_indices.size();
        }

        // Write to file
        if (binary)
        {
            if (!values.empty())
            {
                ofs.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(T));
                ofs.write(reinterpret_cast<const char*>(col_indices.data()), col_indices.size() * sizeof(int));
            }
            ofs.write(reinterpret_cast<const char*>(indptr.data()), indptr.size() * sizeof(long long));
        }
        else
        {
            ofs << std::fixed << std::scientific << std::setprecision(8);
            for (auto const& val : values)
            {
                write_data(ofs, val);
            }
            ofs << std::endl;

            for (auto const& col : col_indices)
            {
                ofs << " " << col;
            }
            ofs << std::endl;

            for (auto const& i : indptr)
            {
                ofs << " " << i;
            }
            ofs << std::endl;
        }
    }
}

template void ModuleIO::output_single_R<double>(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, double>>& XR,
    const double& sparse_threshold,
    const bool& binary,
    const Parallel_Orbitals& pv,
    const bool& reduce);

template void ModuleIO::output_single_R<std::complex<double>>(std::ofstream& ofs,
    const std::map<size_t, std::map<size_t, std::complex<double>>>& XR,
    const double& sparse_threshold,
    const bool& binary,
    const Parallel_Orbitals& pv,
    const bool& reduce);
