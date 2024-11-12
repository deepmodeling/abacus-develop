#pragma once
#include "module_io/restart_exx_csr.h"
#include "module_cell/unitcell.h"
#include "module_io/csr_reader.h"
#include "module_io/write_HS_sparse.h"
#include "module_ri/serialization_cereal.h"
#include <RI/global/Tensor.h>
#include <map>

namespace ModuleIO
{
    template<typename Tdata>
    void read_Hexxs_csr(const std::string& file_name, const UnitCell& ucell,
        const int nspin, const int nbasis,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs)
    {
        ModuleBase::TITLE("Exx_LRI", "read_Hexxs_csr");
        Hexxs.resize(nspin);
        for (int is = 0;is < nspin;++is)
        {
            ModuleIO::csrFileReader<Tdata> csr(file_name + "_" + std::to_string(is) + ".csr");
            int nR = csr.getNumberOfR();
            int nbasis = csr.getMatrixDimension();
            // allocate Hexxs[is]
            for (int iat1 = 0; iat1 < ucell.nat; ++iat1) {
                for (int iat2 = 0;iat2 < ucell.nat;++iat2) {
                    for (int iR = 0;iR < nR;++iR)
                    {
                        const std::vector<int>& R = csr.getRCoordinate(iR);
                        TC dR({ R[0], R[1], R[2] });
                        Hexxs[is][iat1][{iat2, dR}] = RI::Tensor<Tdata>({ static_cast<size_t>(ucell.atoms[ucell.iat2it[iat1]].nw), static_cast<size_t>(ucell.atoms[ucell.iat2it[iat2]].nw) });
                    }
                }
            }
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (int iR = 0;iR < nR;++iR)
            {
                const std::vector<int>& R = csr.getRCoordinate(iR);
                const SparseMatrix<Tdata>& matrix = csr.getMatrix(iR);
                for (auto& ijv : matrix.getElements())
                {
                    const int& npol = ucell.get_npol();
                    const int& i = ijv.first.first * npol;
                    const int& j = ijv.first.second * npol;
                    Hexxs.at(is).at(ucell.iwt2iat[i]).at({ ucell.iwt2iat[j], { R[0], R[1], R[2] } })(ucell.iwt2iw[i] / npol, ucell.iwt2iw[j] / npol) = ijv.second;
                }
            }
        }
    }

    template<typename Tdata>
    void read_Hexxs_cereal(const std::string& file_name, 
        std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs)
    {
        ModuleBase::TITLE("Exx_LRI", "read_Hexxs_cereal");
        ModuleBase::timer::tick("Exx_LRI", "read_Hexxs_cereal");
        std::ifstream ifs(file_name, std::ios::binary);
        cereal::BinaryInputArchive iar(ifs);
        iar(Hexxs);
        ModuleBase::timer::tick("Exx_LRI", "read_Hexxs_cereal");
    }

    template<typename Tdata>
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>>
        calculate_RI_Tensor_sparse(const double& sparse_threshold,
            const std::map<int, std::map<TAC, RI::Tensor<Tdata>>>& Hexxs,
            const UnitCell& ucell)
    {
        ModuleBase::TITLE("Exx_LRI_Interface", "calculate_HContainer_sparse_d");
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>> target;
        for (auto& a1_a2R_data : Hexxs)
        {
            int iat1 = a1_a2R_data.first;
            for (auto& a2R_data : a1_a2R_data.second)
            {
                int iat2 = a2R_data.first.first;
                int nw1 = ucell.atoms[ucell.iat2it[iat1]].nw;
                int nw2 = ucell.atoms[ucell.iat2it[iat2]].nw;
                int start1 = ucell.atoms[ucell.iat2it[iat1]].stapos_wf / ucell.get_npol() + ucell.iat2ia[iat1] * nw1;
                int start2 = ucell.atoms[ucell.iat2it[iat2]].stapos_wf / ucell.get_npol() + ucell.iat2ia[iat2] * nw2;

                const TC& R = a2R_data.first.second;
                auto& matrix = a2R_data.second;
                Abfs::Vector3_Order<int> dR(R[0], R[1], R[2]);
                for (int i = 0;i < nw1;++i) {
                    for (int j = 0;j < nw2;++j) {
                        target[dR][start1 + i][start2 + j] = ((std::abs(matrix(i, j)) > sparse_threshold) ? matrix(i, j) : static_cast<Tdata>(0));
                    }
                }
            }
        }
        return target;
    }
    template<typename Tdata>
    void write_Hexxs_csr(const std::string& file_name, const UnitCell& ucell,
        const std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs)
    {
        ModuleBase::TITLE("Exx_LRI", "write_Hexxs_csr");
        std::set<Abfs::Vector3_Order<int>> all_R_coor;
        double sparse_threshold = 1e-10;
        for (int is = 0;is < Hexxs.size();++is)
        {
            for (const auto& HexxA : Hexxs[is])
            {
                const int iat0 = HexxA.first;
                for (const auto& HexxB : HexxA.second)
                {
                    const int iat1 = HexxB.first.first;
                    const Abfs::Vector3_Order<int> R = RI_Util::array3_to_Vector3(HexxB.first.second);
                    all_R_coor.insert(R);
                }
            }
            ModuleIO::save_sparse(
                calculate_RI_Tensor_sparse(sparse_threshold, Hexxs[is], ucell),
                all_R_coor,
                sparse_threshold,
                false, //binary
                file_name + "_" + std::to_string(is) + ".csr",
                Parallel_Orbitals(),
                "Hexxs_" + std::to_string(is),
                -1,
                false);  //no reduce, one file for each process
        }
    }
}
