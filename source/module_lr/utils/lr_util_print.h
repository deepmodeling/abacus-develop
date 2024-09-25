#pragma once
#include <ATen/tensor.h>
#include "module_psi/psi.h"
#ifdef __EXX
#include <RI/global/Tensor.h>
#endif
namespace LR_Util
{
    template<typename T>
    constexpr T filter(const T& v, const double& threshold = 1e-10)
    {
        return (std::abs(v) > threshold ? v : 0);
    }

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
                std::cout << filter(psi(ib, ik, i)) << " ";
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
                    ofs << filter(psi(ib, ik, i)) << " ";
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
                    std::cout << filter(psi(ik, ib, i)) << " ";
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
                std::cout << filter(t.data<T>()[j * pmat->get_row_size() + i]) << " ";
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


    using TC = std::array<int, 3>;
    using TAC = std::pair<int, TC>;
    template<typename T>
    using TLRIX = std::map<int, std::map<TAC, std::vector<T>>>;

    template<typename T>
    void print_CsX(const TLRIX<T>& CsX, const int nvirt, const std::string& label, const double& threshold = 1e-10)
    {
        std::cout << label << "\n";
        for (const auto& tmp1 : CsX)
        {
            const int& iat1 = tmp1.first;
            for (const auto& tmp2 : tmp1.second)
            {
                const int& iat2 = tmp2.first.first;
                const auto& R = tmp2.first.second;
                auto& t = tmp2.second;
                const int nocc = t.size() / nvirt;
                std::cout << "iat1=" << iat1 << " iat2=" << iat2 << " R=(" << R[0] << " " << R[1] << " " << R[2] << ")\n";
                for (int io = 0;io < nocc;++io)
                {
                    for (int iv = 0;iv < nvirt;++iv) { std::cout << filter(t[io * nvirt + iv]) << " "; }
                    std::cout << "\n";
                }
                std::cout << "\n";
            }
        }
    }

#ifdef __EXX
    template <typename T>
    using TLRI = std::map<int, std::map<TAC, RI::Tensor<T>>>;
    template<typename T>
    void print_CV(const TLRI<T>& Cs, const std::string& label, const double& threshold = 1e-10)
    {
        std::cout << label << "\n";
        for (const auto& tmp1 : Cs)
        {
            const int& iat1 = tmp1.first;
            for (const auto& tmp2 : tmp1.second)
            {
                const int& iat2 = tmp2.first.first;
                const auto& R = tmp2.first.second;
                auto& t = tmp2.second;
                if (R != TC({ 0, 0, 0 })) {continue;}   // for test
                std::cout << "iat1=" << iat1 << " iat2=" << iat2 << " R=(" << R[0] << " " << R[1] << " " << R[2] << ")\n";
                if (t.shape.size() == 2)
                {
                    for (int iabf1 = 0;iabf1 < t.shape[0];++iabf1)
                    {
                        for (int iabf2 = 0;iabf2 < t.shape[1];++iabf2)
                        {
                            std::cout << filter(t(iabf1, iabf2)) << " ";
                        }
                        std::cout << "\n";
                    }
                }
                else if (t.shape.size() == 3)
                {
                    const int nabf = t.shape[0];
                    const int nw1 = t.shape[1];
                    const int nw2 = t.shape[2];
                    std::cout << "size: " << nabf << " " << nw1 << " " << nw2 << "\n";
                    for (int iabf = 0;iabf < nabf;++iabf)
                    {
                        for (int iw1 = 0;iw1 < nw1;++iw1)
                        {
                            for (int iw2 = 0;iw2 < nw2;++iw2)
                            {
                                std::cout << filter(t(iabf, iw1, iw2)) << " ";
                            }
                            std::cout << "\n";
                        }
                        std::cout << "\n";
                    }
                }
            }
        }
    }
#endif
}