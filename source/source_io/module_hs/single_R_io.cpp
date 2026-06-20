#include "single_R_io.h"
#include "source_base/parallel_reduce.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"

#include <complex>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <vector>

inline void write_data(std::ofstream& ofs, const double& data)
{
    ofs << " " << std::fixed << std::scientific << std::setprecision(16) << data;
}
inline void write_data(std::ofstream& ofs, const std::complex<double>& data)
{
    ofs << " (" << data.real() << "," << data.imag() << ")";
}

template<typename T>
void ModuleIO::output_single_R(std::ofstream& ofs,
    const SparseRBlock<T>& XR,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options)
{
    const int nlocal = pv.get_global_row_size();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::output_single_R",
                                 "Parallel_Orbitals global row size must be positive.");
    }

    std::vector<long long> indptr;
    indptr.reserve(nlocal + 1);
    indptr.push_back(0);

    std::stringstream tem1;
    tem1 << options.temp_dir << std::to_string(GlobalV::DRANK)
         << "temp_sparse_indices.dat";
    std::ofstream ofs_tem1;
    std::ifstream ifs_tem1;

    if (!options.reduce || GlobalV::DRANK == 0)
    {
        if (options.binary)
        {
            ofs_tem1.open(tem1.str().c_str(), std::ios::binary);
        }
        else
        {
            ofs_tem1.open(tem1.str().c_str());
        }
    }

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
                        ofs_tem1.write(reinterpret_cast<char *>(&col), sizeof(int));
                    }
                    else
                    {
                        write_data(ofs, line[col]);
                        ofs_tem1 << " " << col;
                    }

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
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str(), std::ios::binary);
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs.write(reinterpret_cast<char *>(&i), sizeof(long long));
            }
        }
        else
        {
            ofs << std::endl;
            ofs_tem1 << std::endl;
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str());
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs << " " << i;
            }
            ofs << std::endl;
        }

        std::remove(tem1.str().c_str());
    }
}

template void ModuleIO::output_single_R<double>(std::ofstream& ofs,
    const SparseRBlock<double>& XR,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options);

template void ModuleIO::output_single_R<std::complex<double>>(std::ofstream& ofs,
    const SparseRBlock<std::complex<double>>& XR,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options);
