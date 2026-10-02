#include "lat_r_csr.h"
#include "source_base/parallel_reduce.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"

#include <complex>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

inline void write_data(std::ofstream& ofs, const double& data, const int precision)
{
    ofs << " " << std::fixed << std::scientific << std::setprecision(precision) << data;
}
inline void write_data(std::ofstream& ofs, const std::complex<double>& data, const int precision)
{
    ofs << " (" << std::fixed << std::scientific << std::setprecision(precision)
        << data.real() << "," << data.imag() << ")";
}

template<typename T>
void ModuleIO::save_lat_r(std::ofstream& ofs,
    const SparseRBlock<T>& XR,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options)
{
    const int nlocal = pv.get_global_row_size();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::save_lat_r",
                                 "Parallel_Orbitals global row size must be positive.");
    }

    std::vector<long long> indptr;
    indptr.reserve(nlocal + 1);
    indptr.push_back(0);

    std::vector<int> col_indices;

    std::vector<T> line(nlocal);
    for(int row = 0; row < nlocal; ++row)
    {
        ModuleBase::GlobalFunc::ZEROS(line.data(), nlocal);

        if (!options.reduce || pv.global2local_row(row) >= 0)
        {
            auto iter = XR.find(row);
            if (iter != XR.end())
            {
                for (auto &value : iter->second)
                {
                    if (value.first >= static_cast<size_t>(nlocal))
                    {
                        std::cerr << "Sparse column index out of range." << std::endl;
                        ModuleBase::WARNING_QUIT("ModuleIO::save_lat_r",
                                                 "Sparse column index out of range.");
                    }
                    line[value.first] = value.second;
                }
            }
        }

        if (options.reduce)
        {
            Parallel_Reduce::reduce_all(line.data(), nlocal);
        }

        if (!options.reduce || GlobalV::DRANK == 0)
        {
            long long nonzeros_count = 0;
            for (int col = 0; col < nlocal; ++col)
            {
                if (std::abs(line[col]) > options.threshold)
                {
                    if (options.binary)
                    {
                        ofs.write(reinterpret_cast<char*>(&line[col]), sizeof(T));
                    }
                    else
                    {
                        write_data(ofs, line[col], options.precision);
                    }
                    col_indices.push_back(col);
                    nonzeros_count++;
                }

            }
            nonzeros_count += indptr.back();
            indptr.push_back(nonzeros_count);
        }
    }

    if (!options.reduce || GlobalV::DRANK == 0)
    {
        if (options.binary)
        {
            for (int col : col_indices)
            {
                ofs.write(reinterpret_cast<char*>(&col), sizeof(int));
            }
            for (auto &i : indptr)
            {
                ofs.write(reinterpret_cast<char *>(&i), sizeof(long long));
            }
        }
        else
        {
            ofs << std::endl;
            for (int col : col_indices)
            {
                ofs << " " << col;
            }
            ofs << std::endl;
            for (auto &i : indptr)
            {
                ofs << " " << i;
            }
            ofs << std::endl;
        }
    }
}

template void ModuleIO::save_lat_r<double>(std::ofstream& ofs,
    const SparseRBlock<double>& XR,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options);

template void ModuleIO::save_lat_r<std::complex<double>>(std::ofstream& ofs,
    const SparseRBlock<std::complex<double>>& XR,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options);
