#pragma once
#include "module_elecstate/module_dm/density_matrix.h"
namespace LR_Util
{
    template<typename TR>
    void print_HR(const hamilt::HContainer<TR>& HR, const int& nat, const std::string& label, const double& threshold = 1e-10)
    {
        std::cout << label << "\n";
        for (int ia = 0;ia < nat;ia++)
            for (int ja = 0;ja < nat;ja++)
            {
                auto ap = HR.find_pair(ia, ja);
                for (int iR = 0;iR < ap->get_R_size();++iR)
                {
                    std::cout << "atom pair (" << ia << ", " << ja << "),  "
                        << "R=(" << ap->get_R_index(iR)[0] << ", " << ap->get_R_index(iR)[1] << ", " << ap->get_R_index(iR)[2] << "): \n";
                    auto& mat = ap->get_HR_values(iR);
                    std::cout << "rowsize=" << ap->get_row_size() << ", colsize=" << ap->get_col_size() << "\n";
                    for (int i = 0;i < ap->get_row_size();++i)
                    {
                        for (int j = 0;j < ap->get_col_size();++j)
                        {
                            auto& v = mat.get_value(i, j);
                            std::cout << (std::abs(v) > threshold ? v : 0) << " ";
                        }
                        std::cout << "\n";
                    }
                }
            }
    }
    template <typename TK, typename TR>
    void print_DMR(const elecstate::DensityMatrix<TK, TR>& DMR, const int& nat, const std::string& label, const double& threshold = 1e-10)
    {
        std::cout << label << "\n";
        int is = 0;
        for (auto& dr : DMR.get_DMR_vector())
            print_HR(*dr, nat, "DMR[" + std::to_string(is++) + "]", threshold);
    }
    void get_DMR_real_imag_part(const elecstate::DensityMatrix<std::complex<double>, std::complex<double>>& DMR,
        elecstate::DensityMatrix<std::complex<double>, double>& DMR_real,
        const int& nat,
        const char& type = 'R');
    void set_HR_real_imag_part(const hamilt::HContainer<double>& HR_real,
        hamilt::HContainer<std::complex<double>>& HR,
        const int& nat,
        const char& type = 'R');
}