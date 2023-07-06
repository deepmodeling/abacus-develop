#include "single_R_io.h"
#include "module_base/parallel_reduce.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"

void ModuleIO::output_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, double>> &XR, const double &sparse_threshold, const bool &binary, const Parallel_Orbitals &pv)
{
    double *line = nullptr;
    std::vector<int> indptr;
    indptr.reserve(GlobalV::NLOCAL + 1);
    indptr.push_back(0);

    std::stringstream tem1;
    tem1 << GlobalV::global_out_dir << "temp_sparse_indices.dat";
    std::ofstream ofs_tem1;
    std::ifstream ifs_tem1;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.open(tem1.str().c_str(), std::ios::binary);
        }
        else
        {
            ofs_tem1.open(tem1.str().c_str());
        }
    }

    line = new double[GlobalV::NLOCAL];
    for(int row = 0; row < GlobalV::NLOCAL; ++row)
    {
        // line = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(line, GlobalV::NLOCAL);

        if(pv.global2local_row(row) >= 0)
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

        Parallel_Reduce::reduce_double_all(line, GlobalV::NLOCAL);

        if(GlobalV::DRANK == 0)
        {
            int nonzeros_count = 0;
            for (int col = 0; col < GlobalV::NLOCAL; ++col)
            {
                if (std::abs(line[col]) > sparse_threshold)
                {
                    if (binary)
                    {
                        ofs.write(reinterpret_cast<char *>(&line[col]), sizeof(double));
                        ofs_tem1.write(reinterpret_cast<char *>(&col), sizeof(int));
                    }
                    else
                    {
                        ofs << " " << std::fixed << std::scientific << std::setprecision(8) << line[col];
                        ofs_tem1 << " " << col;
                    }

                    nonzeros_count++;

                }

            }
            nonzeros_count += indptr.back();
            indptr.push_back(nonzeros_count);
        }

        // delete[] line;
        // line = nullptr;

    }

    delete[] line;
    line = nullptr;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str(), std::ios::binary);
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs.write(reinterpret_cast<char *>(&i), sizeof(int));
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

void ModuleIO::output_soc_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, std::complex<double>>> &XR, const double &sparse_threshold, const bool &binary, const Parallel_Orbitals &pv)
{
    std::complex<double> *line = nullptr;
    std::vector<int> indptr;
    indptr.reserve(GlobalV::NLOCAL + 1);
    indptr.push_back(0);

    std::stringstream tem1;
    tem1 << GlobalV::global_out_dir << "temp_sparse_indices.dat";
    std::ofstream ofs_tem1;
    std::ifstream ifs_tem1;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.open(tem1.str().c_str(), std::ios::binary);
        }
        else
        {
            ofs_tem1.open(tem1.str().c_str());
        }
    }

    line = new std::complex<double>[GlobalV::NLOCAL];
    for(int row = 0; row < GlobalV::NLOCAL; ++row)
    {
        // line = new std::complex<double>[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(line, GlobalV::NLOCAL);

        if(pv.global2local_row(row) >= 0)
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

        Parallel_Reduce::reduce_complex_double_all(line, GlobalV::NLOCAL);

        if (GlobalV::DRANK == 0)
        {
            int nonzeros_count = 0;
            for (int col = 0; col < GlobalV::NLOCAL; ++col)
            {
                if (std::abs(line[col]) > sparse_threshold)
                {
                    if (binary)
                    {
                        ofs.write(reinterpret_cast<char *>(&line[col]), sizeof(std::complex<double>));
                        ofs_tem1.write(reinterpret_cast<char *>(&col), sizeof(int));
                    }
                    else
                    {
                        ofs << " (" << std::fixed << std::scientific << std::setprecision(8) << line[col].real() << "," 
                                    << std::fixed << std::scientific << std::setprecision(8) << line[col].imag() << ")";
                        ofs_tem1 << " " << col;
                    }

                    nonzeros_count++;

                }

            }
            nonzeros_count += indptr.back();
            indptr.push_back(nonzeros_count);
        }

        // delete[] line;
        // line = nullptr;

    }

    delete[] line;
    line = nullptr;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str(), std::ios::binary);
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs.write(reinterpret_cast<char *>(&i), sizeof(int));
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
