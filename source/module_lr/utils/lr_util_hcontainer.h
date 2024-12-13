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

    template <typename T, typename TR>
    void initialize_HR(hamilt::HContainer<TR>& hR,
                       const UnitCell& ucell,
                       const Grid_Driver& gd,
                       const std::vector<double>& orb_cutoff)
    {
        const auto& pmat = *hR.get_paraV();
        for (int iat1 = 0; iat1 < ucell.nat; iat1++)
        {
            auto tau1 = ucell.get_tau(iat1);
            int T1, I1;
            ucell.iat2iait(iat1, &I1, &T1);
            AdjacentAtomInfo adjs;
            gd.Find_atom(ucell, tau1, T1, I1, &adjs);
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T2 = adjs.ntype[ad];
                const int I2 = adjs.natom[ad];
                int iat2 = ucell.itia2iat(T2, I2);
                if (pmat.get_row_size(iat1) <= 0 || pmat.get_col_size(iat2) <= 0) { continue; }
                const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
                if (ucell.cal_dtau(iat1, iat2, R_index).norm() * ucell.lat0 >= orb_cutoff[T1] + orb_cutoff[T2]) { continue; }
                hamilt::AtomPair<TR> tmp(iat1, iat2, R_index.x, R_index.y, R_index.z, &pmat);
                hR.insert_pair(tmp);
            }
        }
        hR.allocate(nullptr, true);
        // hR.set_paraV(&pmat);
        if (std::is_same<T, double>::value) { hR.fix_gamma(); }
    }
    template <typename T, typename TR>
    void initialize_DMR(elecstate::DensityMatrix<T, TR>& dm,
                        const Parallel_Orbitals& pmat,
                        const UnitCell& ucell,
                        const Grid_Driver& gd,
                        const std::vector<double>& orb_cutoff)
    {
        hamilt::HContainer<TR> hR_tmp(&pmat);
        initialize_HR<T, TR>(hR_tmp, ucell, gd, orb_cutoff);
        dm.init_DMR(hR_tmp);
    }
}