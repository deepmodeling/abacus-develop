#include "lat_r_csr.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/parallel_reduce.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/module_out/sparse_matrix.h"

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
    const Parallel_Blkcyclic& layout,
    const SparseWriteOptions& options)
{
    const int nlocal = layout.nrglb();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::save_lat_r",
                                 "Block-cyclic global row size must be positive.");
    }

    std::vector<long long> indptr;
    indptr.reserve(nlocal + 1);
    indptr.push_back(0);

    std::vector<int> col_indices;

    std::vector<T> line(nlocal);
    for(int row = 0; row < nlocal; ++row)
    {
        ModuleBase::GlobalFunc::ZEROS(line.data(), nlocal);

        if (!options.reduce || layout.irloc(row) >= 0)
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
    const Parallel_Blkcyclic& layout,
    const SparseWriteOptions& options);

template void ModuleIO::save_lat_r<std::complex<double>>(std::ofstream& ofs,
    const SparseRBlock<std::complex<double>>& XR,
    const Parallel_Blkcyclic& layout,
    const SparseWriteOptions& options);

void ModuleIO::write_csr_banner(std::ofstream& ofs)
{
    ofs << " #----------------------------------------------------------------------#" << std::endl;
    ofs << " #                               CSR Format                             #" << std::endl;
    ofs << " # The outer loop corresponds to the number of Bravais lattice vectors. #" << std::endl;
    ofs << " # The first line contains the index of the Bravais lattice vector      #" << std::endl;
    ofs << " # (Rx, Ry, Rz), followed by the number of non-zero elements.           #" << std::endl;
    ofs << " # The subsequent lines consist of three blocks of data, which are      #" << std::endl;
    ofs << " # values, column indices, row pointers.                                #" << std::endl;
    ofs << " #----------------------------------------------------------------------#" << std::endl;
    ofs << std::endl;
}

template <typename T>
void ModuleIO::save_lat_r_csr(std::ofstream& ofs,
    const SparseRBlock<T>& XR,
    const Parallel_Blkcyclic& layout,
    const SparseWriteOptions& options)
{
    const int nlocal = layout.nrglb();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::save_lat_r_csr",
                                 "Block-cyclic global row size must be positive.");
    }

    ModuleIO::SparseMatrix<T> sparse_matrix(nlocal, nlocal);
    sparse_matrix.setSparseThreshold(options.threshold);

    std::vector<T> line(nlocal);
    for (int row = 0; row < nlocal; ++row)
    {
        ModuleBase::GlobalFunc::ZEROS(line.data(), nlocal);

        if (!options.reduce || layout.irloc(row) >= 0)
        {
            auto iter = XR.find(row);
            if (iter != XR.end())
            {
                for (const auto& value: iter->second)
                {
                    if (value.first >= static_cast<size_t>(nlocal))
                    {
                        ModuleBase::WARNING_QUIT("ModuleIO::save_lat_r_csr",
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

        if (!options.reduce || ofs.is_open())
        {
            for (int col = 0; col < nlocal; ++col)
            {
                sparse_matrix.insert(row, col, line[col]);
            }
        }
    }

    if (!options.reduce || ofs.is_open())
    {
        sparse_matrix.printToCSR(ofs, options.precision);
    }
}

template void ModuleIO::save_lat_r_csr<double>(std::ofstream& ofs,
    const SparseRBlock<double>& XR,
    const Parallel_Blkcyclic& layout,
    const SparseWriteOptions& options);

template void ModuleIO::save_lat_r_csr<std::complex<double>>(std::ofstream& ofs,
    const SparseRBlock<std::complex<double>>& XR,
    const Parallel_Blkcyclic& layout,
    const SparseWriteOptions& options);
