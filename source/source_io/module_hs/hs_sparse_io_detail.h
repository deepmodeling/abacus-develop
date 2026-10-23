#ifndef HS_SPARSE_IO_DETAIL_H
#define HS_SPARSE_IO_DETAIL_H

#include "hs_sparse_io.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

/**
 * @brief Internal helpers shared by hs_sparse_io.cpp and dhs_sparse_writer.cpp.
 *        Declared here (instead of duplicated in each anonymous namespace) so
 *        the threshold/nonzero logic cannot diverge between the two writers.
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
} // namespace detail
} // namespace ModuleIO

#endif // HS_SPARSE_IO_DETAIL_H
