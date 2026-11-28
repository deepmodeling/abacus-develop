#include "restart_exx_csr.h"

#include "source_base/global_function.h"
#include "source_cell/unitcell.h"
#include "source_base/module_out/csr_reader.h"
#include "source_io/module_hs/hs_sparse_io.h"
#include "source_lcao/module_ri/serialization_cereal.h"
#include <cereal/types/utility.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/valarray.hpp>
#include <RI/global/Tensor.h>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

namespace ModuleIO
{
    std::string exx_csr_filename(const int ispin, const bool numbered, const int istep)
    {
        std::string name = "hrs" + std::to_string(ispin + 1);
        if (numbered && istep >= 0)
        {
            name += "g" + std::to_string(istep + 1);
        }
        name += "_exx.csr";
        return name;
    }

    template<typename Tdata>
    void read_Hexxs_csr(const std::string& file_dir, const UnitCell& ucell,
        const int nspin, const int nbasis,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs,
        double& eexx)
    {
        ModuleBase::TITLE("ModuleIO", "read_Hexxs_csr");
        Hexxs.clear();
        Hexxs.resize(nspin);
        eexx = 0.0;
        for (int is = 0;is < nspin;++is)
        {
            const std::string path = file_dir + exx_csr_filename(is, false, -1);
            // A single-frame file carries exactly one self-describing header.
            // out_app_flag = 1 appends one header per written ionic step, and
            // csrFileReader would silently parse only the first frame, so stop
            // here and tell the user how to produce a readable restart file.
            int frame_count = 0;
            {
                std::ifstream frame_stream(path);
                std::string frame_line;
                while (std::getline(frame_stream, frame_line))
                {
                    if (frame_line.rfind("Ionic_Step", 0) == 0)
                    {
                        ++frame_count;
                    }
                }
            }
            if (frame_count > 1)
            {
                ModuleBase::WARNING_QUIT("read_Hexxs_csr",
                    "Found " + std::to_string(frame_count) + " ionic-step records in " + path
                    + "; init_exx = file expects a single-frame file written with out_app_flag = 0"
                      " (or, for out_app_flag = 2, the numbered hrs{is}g{step}_exx.csr file).");
            }
            ModuleIO::csrFileReader<Tdata> csr(path);
            int nR = csr.getNumberOfR();
            if (csr.getMatrixDimension() != nbasis)
            {
                ModuleBase::WARNING_QUIT("read_Hexxs_csr",
                                         "EXX Hexx(R) matrix dimension does not match the LCAO basis: " + path);
            }
            // Parse E_exx from the description line of the first spin file.
            // It is mandatory: a CSR file written without out_mat_exx carries no
            // E_exx and cannot be used to restart the EXX outer loop.
            if (is == 0)
            {
                const std::string& desc = csr.getDescription();
                const std::string tag = "E_exx = ";
                const size_t pos = desc.find(tag);
                if (pos == std::string::npos)
                {
                    ModuleBase::WARNING_QUIT("read_Hexxs_csr",
                        "E_exx is missing from the CSR description line of " + path
                        + "; the file was probably written without out_mat_exx = 1.");
                }
                std::istringstream iss(desc.substr(pos + tag.size()));
                if (!(iss >> eexx))
                {
                    ModuleBase::WARNING_QUIT("read_Hexxs_csr",
                        "Failed to parse E_exx from CSR description: " + desc);
                }
            }
            const std::string confirmation = "Read HexxR from " + path;
            std::cout << confirmation << std::endl;
            ModuleBase::TITLE("read_Hexxs_csr", confirmation, false);
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
    void write_Hexxs_csr(const std::string& out_dir, const UnitCell& ucell,
        const std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs,
        const double eexx,
        const int precision,
        const int istep,
        const int out_app_flag)
    {
        ModuleBase::TITLE("ModuleIO", "write_Hexxs_csr");
        std::set<Abfs::Vector3_Order<int>> all_R_coor;
        double sparse_threshold = 1e-10;
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
        Parallel_Orbitals pv;
        pv.para_2d().set_serial(matrix_dimension, matrix_dimension);

        // Format E_exx for embedding in the CSR description line
        std::ostringstream eexx_ss;
        eexx_ss << ", E_exx = " << std::setprecision(precision) << eexx << " Ry";

        // out_app_flag = 2 embeds the one-based ionic step in the file name;
        // modes 0 and 1 keep the fixed name (mode 1 appends later steps).
        const bool numbered = (istep >= 0) && (out_app_flag == 2);
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
            ModuleIO::SparseWriteOptions options;
            options.filename = out_dir + exx_csr_filename(is, numbered, istep);
            options.label = "Hexxs_" + std::to_string(is);
            options.threshold = sparse_threshold;
            options.binary = false;
            options.precision = precision;
            options.istep = istep;
            options.out_app_flag = out_app_flag;
            // Hexxs(R) is fully replicated on every rank, so its rows must NOT be
            // summed across ranks; use reduce = false and let only rank 0 write.
            options.reduce = false;
            options.csr_format = true;
            options.ucell = &ucell;
            options.nspin = static_cast<int>(Hexxs.size());
            options.ispin = is;
            options.description_suffix = eexx_ss.str();
            // Only one rank writes the replicated matrix to avoid concurrent
            // writes and duplicated (reduced) data.
            if (GlobalV::DRANK == 0)
            {
                ModuleIO::save_sparse(
                    calculate_RI_Tensor_sparse(sparse_threshold, Hexxs[is], ucell),
                    all_R_coor,
                    pv,
                    options);
            }
        }
    }

template void read_Hexxs_csr<double>(const std::string&,
                                     const UnitCell&,
                                     const int,
                                     const int,
                                     std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>&,
                                     double&);
template void read_Hexxs_csr<std::complex<double>>(const std::string&,
                                                   const UnitCell&,
                                                   const int,
                                                   const int,
                                                   std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>&,
                                                   double&);
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
                                      const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>&,
                                      const double,
                                      const int,
                                      const int,
                                      const int);
template void write_Hexxs_csr<std::complex<double>>(const std::string&,
                                                    const UnitCell&,
                                                    const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>&,
                                                    const double,
                                                    const int,
                                                    const int,
                                                    const int);
} // namespace ModuleIO
