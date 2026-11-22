#ifndef HS_SPARSE_IO_DETAIL_H
#define HS_SPARSE_IO_DETAIL_H

#include "hs_sparse_io.h"
#include "source_base/module_out/sparse_matrix.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"
#include "source_cell/ucell_io.h"
#include "source_io/module_output/provenance_header.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

/**
 * @brief Internal helpers shared by the R-space sparse text writers.
 *        Declared here (instead of duplicated in each anonymous namespace) so
 *        the header layout, R records and file-append rules cannot diverge
 *        between hs_sparse_io.cpp, dhs_sparse_writer.cpp, pos_op_csr.cpp,
 *        hsr_writer.cpp and write_dmr.cpp.
 */
namespace ModuleIO
{
namespace detail
{
/**
 * @brief Count non-zero elements per R-block of a sparse R-matrix.
 * @param smat       Sparse matrix keyed by R-coordinate.
 * @param all_R_coor Ordered set of R-coordinates to scan.
 * @param threshold  Magnitude below which an element is treated as zero.
 * @param reduce     If true, allreduce the per-R counts across MPI ranks.
 * @return Per-R nonzero counts aligned with all_R_coor order.
 */
template <typename Tdata>
std::vector<long long> count_nonzeros_by_R(
    const SparseRMatrix<Tdata>& smat,
    const std::set<RCoordinate>& all_R_coor,
    const double threshold,
    const bool reduce)
{
    std::vector<long long> nonzero_num(all_R_coor.size(), 0);
    int count = 0;
    for (const auto& R_coor: all_R_coor)
    {
        const auto iter = smat.find(R_coor);
        if (iter != smat.end())
        {
            for (const auto& row_loop: iter->second)
            {
                for (const auto& col_value: row_loop.second)
                {
                    if (std::abs(col_value.second) > threshold)
                    {
                        ++nonzero_num[count];
                    }
                }
            }
        }
        ++count;
    }

    if (reduce)
    {
        Parallel_Reduce::reduce_all(nonzero_num.data(), static_cast<int>(nonzero_num.size()));
    }
    return nonzero_num;
}

/**
 * @brief Quit if a sparse-matrix output stream failed to open.
 * @param ofs      Stream to check.
 * @param filename Filename for the error message.
 * @param context Caller tag for the error message.
 */
inline void check_output_file_open(const std::ofstream& ofs,
                                   const std::string& filename,
                                   const std::string& context)
{
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT(context, "Cannot open sparse matrix file: " + filename);
    }
}

/**
 * @brief Decide whether the current ionic step appends to an existing file.
 * @details Mode 1 appends all ionic steps after the first written one; mode 0
 *          and mode 2 always start a fresh stream (mode 2 embeds the step in
 *          the file name). Single-point runs pass istep <= 0 and truncate.
 * @param out_app_flag Output append mode from INPUT.
 * @param istep Zero-based ionic step index.
 * @return True when the stream must be opened in append mode.
 */
inline bool append_ionic_step(const int out_app_flag, const int istep)
{
    return out_app_flag == 1 && istep > 0;
}

/**
 * @brief Build the "<prefix><spin suffix>[g<step>]_nao<ext>" name shared by
 *        every R-space sparse output.
 * @details Centralizing the suffix order keeps H(R), S(R), DM(R), dH/dR,
 *          r(R), T(R), Vxc(R) and L(R) file names consistent and avoids one
 *          writer drifting from the others.
 * @param prefix Property prefix without spin or step, e.g. "dhrx", "rrx".
 * @param spin_suffix Spin token appended to the prefix ("s1", "1") or empty
 *        when the quantity is spin independent.
 * @param numbered Whether to embed the one-based ionic step.
 * @param istep Zero-based ionic step index.
 * @param extension File extension including the leading dot.
 * @return Assembled file name without any directory prefix.
 */
inline std::string make_sparse_R_filename(const std::string& prefix,
                                          const std::string& spin_suffix,
                                          const bool numbered,
                                          const int istep,
                                          const std::string& extension)
{
    if (prefix.empty())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::detail::make_sparse_R_filename",
                                 "File name prefix must not be empty.");
    }
    if (extension.empty() || extension[0] != '.')
    {
        ModuleBase::WARNING_QUIT("ModuleIO::detail::make_sparse_R_filename",
                                 "File extension must start with a dot.");
    }

    std::string name = prefix + spin_suffix;
    if (numbered && istep >= 0)
    {
        name += "g" + std::to_string(istep + 1);
    }
    name += "_nao" + extension;
    return name;
}

/**
 * @brief Metadata shared by every self-describing CSR text writer.
 * @details Callers set all fields explicitly so a new output type cannot
 *          silently inherit the wrong spin or cell metadata.
 */
struct CsrTextHeader
{
    int step;                        ///< One-based ionic step.
    std::string description;         ///< Full " # print ..." matrix description line.
    int nspin;                       ///< Number of spin directions.
    int ispin;                       ///< One-based spin index.
    int nbasis;                      ///< Number of localized basis functions.
    int number_of_R;                 ///< Number of R blocks actually written.
    const UnitCell* ucell;           ///< Unit cell; must not be null.
    std::string representation_note; ///< Optional note written after the cell block.
    bool has_efermi;                 ///< Append the Fermi energy to the spin line.
    double efermi_eV;                ///< Fermi energy in eV.
    bool write_banner;               ///< Write the CSR layout banner.
};

/**
 * @brief Write the common self-describing header of a CSR text file.
 * @param ofs Output text stream.
 * @param header Header fields; ucell must not be null.
 */
inline void write_csr_text_header(std::ofstream& ofs, const CsrTextHeader& header)
{
    if (header.ucell == nullptr)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::detail::write_csr_text_header",
                                 "UnitCell pointer must not be null.");
    }
    if (header.nspin <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::detail::write_csr_text_header",
                                 "Number of spin directions must be positive.");
    }
    if (header.ispin < 1 || header.ispin > header.nspin)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::detail::write_csr_text_header",
                                 "Spin index must be within the spin range.");
    }
    if (header.nbasis <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::detail::write_csr_text_header",
                                 "Number of localized basis functions must be positive.");
    }
    if (header.number_of_R < 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::detail::write_csr_text_header",
                                 "Number of Bravais lattice vectors must not be negative.");
    }

    const std::string time_stamp = ModuleIO::local_timestamp();
    ofs << "Ionic_Step " << header.step << ModuleIO::prov_inline(time_stamp) << std::endl;
    ofs << header.description << std::endl;
    ofs << " " << header.nspin << " # number of spin directions" << std::endl;
    ofs << " " << header.ispin << " # spin index";
    if (header.has_efermi)
    {
        ofs << ", E_Fermi = " << std::setprecision(6) << header.efermi_eV << " eV";
    }
    ofs << std::endl;
    ofs << " " << header.nbasis << " # number of localized basis" << std::endl;
    ofs << " " << header.number_of_R << " # number of Bravais lattice vector R" << std::endl;
    ofs << std::endl;

    ModuleIO::UcellIO::write_ucell(ofs, header.ucell);
    if (!header.representation_note.empty())
    {
        ofs << "# representation: " << header.representation_note << std::endl;
    }
    ofs << std::endl;

    if (header.write_banner)
    {
        ModuleIO::write_csr_banner(ofs);
    }
}

/**
 * @brief Write the three-integer binary CSR header.
 * @param ofs Output binary stream.
 * @param step One-based ionic step.
 * @param nbasis Number of localized basis functions.
 * @param number_of_R Number of R blocks actually written.
 */
inline void write_csr_binary_header(std::ofstream& ofs,
                                    const int step,
                                    const int nbasis,
                                    const int number_of_R)
{
    ofs.write(reinterpret_cast<const char*>(&step), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&nbasis), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&number_of_R), sizeof(int));
}

/**
 * @brief Write the "Rx Ry Rz nnz" record that starts one R block.
 * @param ofs Output stream.
 * @param rx First Bravais lattice vector component.
 * @param ry Second Bravais lattice vector component.
 * @param rz Third Bravais lattice vector component.
 * @param nnz Number of non-zero elements in the block.
 * @param binary Whether the stream uses the binary record layout.
 * @param leading_space Prefix the text record with a space (save_sparse layout).
 */
inline void write_R_record(std::ofstream& ofs,
                           const int rx,
                           const int ry,
                           const int rz,
                           const long long nnz,
                           const bool binary,
                           const bool leading_space)
{
    if (binary)
    {
        if (nnz > static_cast<long long>(std::numeric_limits<int>::max()))
        {
            ModuleBase::WARNING_QUIT("ModuleIO::detail::write_R_record",
                                     "Too many non-zero elements for the binary R record.");
        }
        const int count = static_cast<int>(nnz);
        ofs.write(reinterpret_cast<const char*>(&rx), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&ry), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&rz), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&count), sizeof(int));
    }
    else
    {
        if (leading_space)
        {
            ofs << " ";
        }
        ofs << rx << " " << ry << " " << rz << " " << nnz << std::endl;
    }
}
} // namespace detail
} // namespace ModuleIO

#endif // HS_SPARSE_IO_DETAIL_H
