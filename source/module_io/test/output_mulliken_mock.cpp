#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ModuleIO
{

std::vector<double> read_k(std::string filename, int ik)
{
    std::ifstream skFile(filename);
    // Initialize variables for file parsing
    std::string line;
    int nrow = 0;
    int ncol = 0;

    // find ik before read the following
    while (std::getline(skFile, line))
    {
        if (line.rfind("# k-point", 0) == 0)
        {
            std::istringstream ss(line.substr(9)); // 9 to skip "# k-point"
            int ik_read;
            ss >> ik_read;
            if (ik_read == ik)
            {
                break;
            }
        }
    }

    // Read the file comments, begin parsing when nrow and ncol are found
    while (std::getline(skFile, line))
    {
        if (line.rfind("# nrow", 0) == 0)
        {
            std::istringstream ss(line.substr(7)); // 7 to skip "# nrow "
            ss >> nrow;
        }
        else if (line.rfind("# ncol", 0) == 0)
        {
            std::istringstream ss(line.substr(7)); // 7 to skip "# ncol "
            ss >> ncol;
            break;
        }
    }

    // Initialize a vector of size nrow * ncol with default values of 0.0
    std::vector<double> matrix(nrow * ncol, 0.0);

    // Read the rest of the file and populate the vector
    int index;
    double value;
    while (skFile >> index >> value)
    {
        if (index >= 0 && index < nrow * ncol)
        {
            matrix[index] = value;
        }
    }

    skFile.close();
    return matrix;
}

} // namespace ModuleIO

#include "module_io/output_dmk.h"
#include "module_io/output_sk.h"

namespace ModuleIO
{

template <typename TK>
Output_DMK<TK>::Output_DMK(elecstate::DensityMatrix<TK, double>* p_DM, Parallel_Orbitals* ParaV, int nspin, int nks)
    : p_DM_(p_DM), ParaV_(ParaV), nspin_(nspin), nks_(nks)
{
}

template <typename TK>
TK* Output_DMK<TK>::get_DMK(int ik)
{
    if (this->nspin_ == 1)
    {
        std::vector<double> dmk = read_k("./support/DMK_nspin1", ik);
        /// convert sk to TK
        this->DMK.resize(dmk.size());
        for (size_t i = 0; i < dmk.size(); i++)
        {
            this->DMK[i] = dmk[i];
            // std::cout << "DMK[" << i << "] = " << DMK[i] << std::endl;
        }
    }
    else if (this->nspin_ == 2)
    {
        std::vector<double> dmk = read_k("./support/DMK_nspin2", ik);
        /// convert sk to TK
        this->DMK.resize(dmk.size());
        for (size_t i = 0; i < dmk.size(); i++)
        {
            this->DMK[i] = dmk[i];
            // std::cout << "DMK[" << i << "] = " << DMK[i] << std::endl;
        }
    }
    else if (this->nspin_ == 4)
    {
        std::vector<double> dmk1 = read_k("./support/DMK_nspin2", 0);
        std::vector<double> dmk2 = read_k("./support/DMK_nspin2", 1);
        this->DMK.resize(dmk1.size() * 4, 0.0);
        for (size_t i = 0; i < dmk1.size(); i++)
        {
            /// sk1 is column-major, find ir and irc
            int ir_nspin2 = i % (ParaV_->nrow / 2);
            int ic_nspin2 = i / (ParaV_->nrow / 2);
            int ir_nspin4 = ir_nspin2 * 2;
            int ic_nspin4 = ic_nspin2 * 2;
            int index_nspin4_up = ir_nspin4 + ic_nspin4 * ParaV_->nrow;
            int index_nspin4_dn = ir_nspin4 + 1 + (ic_nspin4 + 1) * ParaV_->nrow;
            this->DMK[index_nspin4_up] = dmk1[i];
            this->DMK[index_nspin4_dn] = dmk2[i];
        }
    }
    return this->DMK.data();
}

template class Output_DMK<double>;
template class Output_DMK<std::complex<double>>;

template <typename TK>
Output_Sk<TK>::Output_Sk(hamilt::Hamilt<TK>* p_hamilt, Parallel_Orbitals* ParaV, int nspin, int nks)
    : p_hamilt_(p_hamilt), ParaV_(ParaV), nspin_(nspin), nks_(nks)
{
}

template <typename TK>
TK* Output_Sk<TK>::get_Sk(int ik)
{
    if (this->nspin_ == 1)
    {
        std::vector<double> sk = read_k("./support/SK_nspin1", ik);
        /// convert sk to TK
        this->SK.resize(sk.size());
        for (size_t i = 0; i < sk.size(); i++)
        {
            this->SK[i] = sk[i];
            // std::cout << "SK[" << i << "] = " << SK[i] << std::endl;
        }
    }
    else if (this->nspin_ == 2)
    {
        std::vector<double> sk = read_k("./support/SK_nspin2", ik);
        /// convert sk to TK
        this->SK.resize(sk.size());
        for (size_t i = 0; i < sk.size(); i++)
        {
            this->SK[i] = sk[i];
            // std::cout << "SK[" << i << "] = " << SK[i] << std::endl;
        }
    }
    else if (this->nspin_ == 4)
    {
        std::vector<double> sk1 = read_k("./support/SK_nspin2", 0);
        std::vector<double> sk2 = read_k("./support/SK_nspin2", 1);
        this->SK.resize(sk1.size() * 4, 0.0);
        for (size_t i = 0; i < sk1.size(); i++)
        {
            /// sk1 is column-major, find ir and irc
            int ir_nspin2 = i % (ParaV_->nrow / 2);
            int ic_nspin2 = i / (ParaV_->nrow / 2);
            int ir_nspin4 = ir_nspin2 * 2;
            int ic_nspin4 = ic_nspin2 * 2;
            int index_nspin4_up = ir_nspin4 + ic_nspin4 * ParaV_->nrow;
            int index_nspin4_dn = ir_nspin4 + 1 + (ic_nspin4 + 1) * ParaV_->nrow;
            this->SK[index_nspin4_up] = sk1[i];
            this->SK[index_nspin4_dn] = sk2[i];
        }
    }
    return this->SK.data();
}

template class Output_Sk<double>;
template class Output_Sk<std::complex<double>>;

} // namespace ModuleIO