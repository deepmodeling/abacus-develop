#include "restart_exx_csr.h"

#include "source_base/global_function.h"
#include "source_cell/unitcell.h"
#include "source_base/module_out/csr_reader.h"
#include "source_base/module_out/sparse_matrix.h"
#include "source_cell/ucell_io.h"
#include "source_io/module_hs/hs_sparse_io.h"
#include "source_lcao/module_ri/serialization_cereal.h"
#include <cereal/types/utility.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/valarray.hpp>
#include <RI/global/Tensor.h>
#include <map>

namespace ModuleIO
{
    template<typename Tdata>
    void read_Hexxs_csr(const std::string& file_name, const UnitCell& ucell,
        const int nspin, const int nbasis,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs)
    {
        ModuleBase::TITLE("ModuleIO", "read_Hexxs_csr");
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
                        Hexxs[is][iat1][{iat2, dR}] = RI::Tensor<Tdata>(
                                {
                                static_cast<size_t>(ucell.atoms[ucell.iat2it[iat1]].nw),
                                static_cast<size_t>(ucell.atoms[ucell.iat2it[iat2]].nw)
                                }
                                );
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
                    Hexxs.at(is).at(ucell.iwt2iat[i]).at(
                            {
                                ucell.iwt2iat[j],
                                { R[0], R[1], R[2] }
                            }
                            )(ucell.iwt2iw[i] / npol, ucell.iwt2iw[j] / npol) = ijv.second;
                }
            }
        }
    }

    template<typename Tdata>
    void read_Hexxs_cereal(const std::string& file_name,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs)
    {
        ModuleBase::TITLE("ModuleIO", "read_Hexxs_cereal");
        ModuleBase::timer::start("Exx_LRI", "read_Hexxs_cereal");
        std::ifstream ifs(file_name, std::ios::binary);
        if(!ifs.is_open())
            { ModuleBase::WARNING_QUIT("read_Hexxs_cereal", file_name+" not found."); }
        cereal::BinaryInputArchive iar(ifs);
        iar(Hexxs);
        ModuleBase::timer::end("Exx_LRI", "read_Hexxs_cereal");
    }

    template<typename Tdata>
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>>
        calculate_RI_Tensor_sparse(const double& sparse_threshold,
            const std::map<int, std::map<TAC, RI::Tensor<Tdata>>>& Hexxs,
            const UnitCell& ucell)
    {
        ModuleBase::TITLE("ModuleIO", "calculate_HContainer_sparse_d");
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
                        if (std::abs(matrix(i, j)) > sparse_threshold)
                        {
                            target[dR][start1 + i][start2 + j] = matrix(i, j);
                        }
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
        ModuleBase::TITLE("ModuleIO", "write_Hexxs_csr");
        const double sparse_threshold = 1e-10;
        int matrix_dimension = 0;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            matrix_dimension += ucell.atoms[it].na * ucell.atoms[it].nw * ucell.get_npol();
        }
        if (matrix_dimension <= 0)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_Hexxs_csr",
                                     "LCAO matrix dimension must be positive.");
        }
        // The Hexx restart files must round-trip through read_Hexxs_csr, whose
        // csrFileReader parses the sparse-file layout shared with the DM(R)
        // restart files: header, lattice block, then CSR blocks per R. The
        // bare save_sparse header cannot be read back.
        for (int is = 0;is < Hexxs.size();++is)
        {
            std::set<Abfs::Vector3_Order<int>> all_R_coor;
            for (const auto& HexxA : Hexxs[is])
            {
                for (const auto& HexxB : HexxA.second)
                {
                    all_R_coor.insert(RI_Util::array3_to_Vector3(HexxB.first.second));
                }
            }
            std::ofstream ofs(file_name + "_" + std::to_string(is) + ".csr");
            if (!ofs.is_open())
            {
                ModuleBase::WARNING_QUIT("ModuleIO::write_Hexxs_csr",
                                         "cannot open " + file_name + "_" + std::to_string(is) + ".csr");
            }
            ofs << " --- Ionic Step 1 ---" << std::endl;
            ofs << " # print exact exchange matrix in real space Hexxs_" << is << "(R)" << std::endl;
            ofs << " " << Hexxs.size() << " # number of spin directions" << std::endl;
            ofs << " " << is + 1 << " # spin index" << std::endl;
            ofs << " " << matrix_dimension << " # number of localized basis" << std::endl;
            ofs << " " << all_R_coor.size() << " # number of Bravais lattice vector R" << std::endl;
            ofs << std::endl;
            ModuleIO::UcellIO::write_ucell(ofs, &ucell);
            ofs << std::endl;

            const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>> sparse_by_R
                = calculate_RI_Tensor_sparse(sparse_threshold, Hexxs[is], ucell);
            for (const auto& R : all_R_coor)
            {
                ModuleIO::SparseMatrix<Tdata> sparse_matrix(matrix_dimension, matrix_dimension);
                sparse_matrix.setSparseThreshold(sparse_threshold);
                const auto itR = sparse_by_R.find(R);
                if (itR != sparse_by_R.end())
                {
                    for (const auto& row_cols : itR->second)
                    {
                        for (const auto& col_val : row_cols.second)
                        {
                            sparse_matrix.insert(row_cols.first, col_val.first, col_val.second);
                        }
                    }
                }
                ofs << " " << R.x << " " << R.y << " " << R.z << " "
                    << sparse_matrix.getNNZ() << std::endl;
                sparse_matrix.printToCSR(ofs, 8);
            }
            ofs.close();
        }
    }

template void read_Hexxs_csr<double>(const std::string&,
                                     const UnitCell&,
                                     const int,
                                     const int,
                                     std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>&);
template void read_Hexxs_csr<std::complex<double>>(const std::string&,
                                                   const UnitCell&,
                                                   const int,
                                                   const int,
                                                   std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>&);
template void read_Hexxs_cereal<double>(const std::string&,
                                        std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>&);
template void read_Hexxs_cereal<std::complex<double>>(const std::string&,
                                                      std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>&);
template std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>
    calculate_RI_Tensor_sparse<double>(const double&,
                                       const std::map<int, std::map<TAC, RI::Tensor<double>>>&,
                                       const UnitCell&);
template std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
    calculate_RI_Tensor_sparse<std::complex<double>>(const double&,
                                                     const std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>&,
                                                     const UnitCell&);
template void write_Hexxs_csr<double>(const std::string&,
                                      const UnitCell&,
                                      const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>&);
template void write_Hexxs_csr<std::complex<double>>(const std::string&,
                                                    const UnitCell&,
                                                    const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>&);
} // namespace ModuleIO
