// #include "module_ri/LRI_CV_Tools.h"
#include <map>
#include <cstddef>
#define IZ(x) int x = 0;
#define SZ(x) std::size_t x = 0;

namespace LRI_CV_Tools
{
    using TC = std::array<int, 3>;
    using TAC = std::pair<int, TC>;
    template <typename T>
    using TLRI = std::map<int, std::map<TAC, RI::Tensor<T>>>;;

    template<typename T>
    inline double absmax(const RI::Tensor<T>& t)
    {
        double res = 0;
        for (int i = 0;i < t.get_shape_all();++i)
        {
            double v = std::abs(std::norm(t.ptr()[i]));
            if (v > res) { res = v; }
        }
        return std::sqrt(res);
    }

    template<typename T>
    TLRI<T> read_Cs_ao(const std::string& file_path, const double& threshold)
    {
        IZ(natom) IZ(ncell) IZ(ia1) IZ(ia2) IZ(ic_1) IZ(ic_2) IZ(ic_3)
        SZ(nw1) SZ(nw2) SZ(nabf)
        std::ifstream infile;
        infile.open(file_path);
        infile >> natom >> ncell;   // no use of ncell

        TLRI<T> Cs;
        while (infile.peek() != EOF)
        {
            infile >> ia1 >> ia2 >> ic_1 >> ic_2 >> ic_3 >> nw1 >> nw2 >> nabf;
            const TC& box = { ic_1, ic_2, ic_3 };
            RI::Tensor<T> tensor_cs({ nabf, nw1, nw2 });
            for (std::size_t i = 0; i != nw1; i++) { for (std::size_t j = 0; j != nw2; j++) { for (std::size_t mu = 0; mu != nabf; mu++) { infile >> tensor_cs(mu, i, j); } } }
            // no screening for data-structure consistency
            if (absmax(tensor_cs) >= threshold) { Cs[ia1 - 1][{ia2 - 1, box}] = tensor_cs; }
        }
        infile.close();
        return Cs;
    }

    template<typename T>
    void write_Cs_ao(const TLRI<T>& Cs, const std::string& file_path)
    {
        std::ofstream outfile;
        outfile.open(file_path);
        outfile << Cs.size() << " " << Cs.at(0).size() / Cs.size() << std::endl;    //natom, ncell
        for (auto& it1 : Cs)
        {
            const int& ia1 = it1.first;
            for (auto& it2 : it1.second)
            {
                const int& ia2 = it2.first.first;
                const auto& box = it2.first.second;
                const auto& tensor_cs = it2.second;
                outfile << ia1 + 1 << " " << ia2 + 1 << " " << box[0] << " " << box[1] << " " << box[2] << std::endl;
                const int& nw1 = tensor_cs.shape[1], nw2 = tensor_cs.shape[2], nabf = tensor_cs.shape[0];
                outfile << nw1 << " " << nw2 << " " << nabf << std::endl;
                for (int i = 0; i != nw1; i++)
                {
                    for (int j = 0; j != nw2; j++)
                    {
                        for (int mu = 0; mu != nabf; mu++) { outfile << tensor_cs(mu, i, j) << " "; }
                        outfile << std::endl;
                    }
                }
            }
        }
        outfile.close();
    }

    template<typename T>
    TLRI<T> read_Vs_abf(const std::string& file_path, const double& threshold)
    {
        IZ(natom) IZ(ncell) IZ(ia1) IZ(ia2) IZ(ic_1) IZ(ic_2) IZ(ic_3)
        SZ(nabf1) SZ(nabf2)
        std::ifstream infile;
        infile.open(file_path);
        infile >> natom >> ncell;   // no use of ncell

        TLRI<T> Vs;
        while (infile.peek() != EOF)
        {
            infile >> ia1 >> ia2 >> ic_1 >> ic_2 >> ic_3 >> nabf1 >> nabf2;
            const TC& box = { ic_1, ic_2, ic_3 };
            RI::Tensor<T> tensor_vs({ nabf1, nabf2 });
            for (std::size_t i = 0; i != nabf1; i++) {
                for (std::size_t j = 0; j != nabf2; j++){ infile >> tensor_vs(i, j);}
}
            if (absmax(tensor_vs) >= threshold) { Vs[ia1 - 1][{ia2 - 1, box}] = tensor_vs; }
            // else ++cs_discard;
        }
        infile.close();
        return Vs;
    }

    template <typename T>
    void write_Vs_abf(const TLRI<T>& Vs, const std::string& file_path)
    {
        std::ofstream outfile;
        outfile.open(file_path);
        outfile << Vs.size() << " " << Vs.at(0).size() / Vs.size() << std::endl;    //natom, ncell
        for (const auto& it1 : Vs)
        {
            const int& ia1 = it1.first;
            for (const auto& it2 : it1.second)
            {
                const int& ia2 = it2.first.first;
                const auto& box = it2.first.second;
                const auto& tensor_v = it2.second;
                outfile << ia1 + 1 << " " << ia2 + 1 << " " << box[0] << " " << box[1] << " " << box[2] << std::endl;
                outfile << tensor_v.shape[0] << " " << tensor_v.shape[1] << std::endl;
                for (int i = 0; i != tensor_v.shape[0]; i++)
                {
                    for (int j = 0; j != tensor_v.shape[1]; j++) {outfile << tensor_v(i, j) << " ";}
                    outfile << std::endl;
                }
            }
        }
        outfile.close();
    }
}