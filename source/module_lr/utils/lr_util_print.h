#pragma once
#include <ATen/tensor.h>
#include "module_psi/psi.h"
namespace LR_Util
{
    template<typename T>
    void print_psi_bandfirst(const psi::Psi<T>& psi, const std::string& label, const int& ib, const double& threshold = 1e-10)
    {
        assert(psi.get_k_first() == 0);
        std::cout << label << ": band " << ib << "\n";
        for (int ik = 0;ik < psi.get_nk();++ik)
        {
            std::cout << "iks=" << ik << "\n";
            for (int i = 0;i < psi.get_nbasis();++i)
            {
                const auto& v = psi(ib, ik, i);
                std::cout << (std::abs(v) > threshold ? v : 0) << " ";
            }
            std::cout << "\n";
        }
    }
    template<typename T>
    void write_psi_bandfirst(const psi::Psi<T>& psi, const std::string& filename, const int& rank, const double& threshold = 1e-10, const int& precision = 8)
    {
        assert(psi.get_k_first() == 0);
        std::ofstream ofs(filename + "_" + std::to_string(rank) + ".dat");
        ofs << std::setprecision(precision) << std::scientific;
        ofs << psi.get_nbands() << " " << psi.get_nk() << " " << psi.get_nbasis() << "\n";
        for (int ib = 0;ib < psi.get_nbands();++ib)
        {
            for (int ik = 0;ik < psi.get_nk();++ik)
            {
                for (int i = 0;i < psi.get_nbasis();++i)
                {
                    const auto& v = psi(ib, ik, i);
                    ofs << (std::abs(v) > threshold ? v : 0) << " ";
                }
                ofs << "\n";
            }
        }
        ofs.close();
    }
    template<typename T>
    psi::Psi<T> read_psi_bandfirst(const std::string& filename, const int& rank)
    {
        std::ifstream ifs(filename + "_" + std::to_string(rank) + ".dat");
        int nbands, nks, nbasis;
        ifs >> nbands >> nks >> nbasis;
        psi::Psi<T> psi(nks, nbands, nbasis, nullptr, false);
        for (int ib = 0;ib < psi.get_nbands();++ib) {
            for (int ik = 0;ik < psi.get_nk();++ik) {
                for (int i = 0;i < psi.get_nbasis();++i) {
                    ifs >> psi(ib, ik, i);
}
}
}
        ifs.close();
        return psi;
    }
    template<typename T >
    void print_psi_kfirst(const psi::Psi<T>& psi, const std::string& label, const double& threshold = 1e-10)
    {
        assert(psi.get_k_first() == 1);
        for (int ik = 0;ik < psi.get_nk();++ik)
        {
            std::cout << label << ": k " << ik << "\n";
            for (int ib = 0;ib < psi.get_nbands();++ib)
            {
                std::cout << "ib=" << ib << ": ";
                for (int i = 0;i < psi.get_nbasis();++i)
                {
                    const auto& v = psi(ik, ib, i);
                    std::cout << (std::abs(v) > threshold ? v : 0) << " ";
                }
                std::cout << "\n";
            }
        }
    }
    template<typename T>
    void print_tensor(const container::Tensor& t, const std::string& label, const Parallel_2D* pmat, const double& threshold = 1e-10)
    {
        std::cout << label << "\n";
        for (int j = 0; j < pmat->get_col_size();++j)
        {
            for (int i = 0;i < pmat->get_row_size();++i)
            {
                const auto& v = t.data<T>()[j * pmat->get_row_size() + i];
                std::cout << (std::abs(v) > threshold ? v : 0) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "\n";
    }
    template<typename T>
    void print_grid_nonzero(T* rho, const int& nrxx, const int& nnz, const std::string& label, const double& threshold = 1e-5)
    {
        std::cout << "first " << nnz << " non-zero elements of " << label << "\n";
        int inz = 0;int i = 0;
        while (inz < nnz && i < nrxx) {
            if (rho[++i] - 0.0 > threshold) { std::cout << rho[i] << " ";++inz; }
};
    }

}