#include "hs_dense_io.h"

#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/global_function.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <type_traits>
#include <vector>

namespace
{
// Gather row i of a distributed 2D-block square matrix into a dense row.
// Under MPI, each rank fills the columns it owns locally; the caller then
// reduces across ranks so rank 0 holds the complete row. Without MPI the
// matrix is local and the whole row is read directly.
template <typename T>
void gather_row(const T* mat,
                const int dim,
                const int i,
                const bool tri,
                const Parallel_2D& pv,
                const std::string& ks_solver,
                std::vector<T>& line)
{
#ifdef __MPI
    std::fill(line.begin(), line.end(), T(0));
    const int ir = pv.global2local_row(i);
    if (ir >= 0)
    {
        // data collection
        for (int j = (tri ? i : 0); j < dim; ++j)
        {
            const int ic = pv.global2local_col(j);
            if (ic >= 0)
            {
                int iic = 0;
                if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(ks_solver))
                {
                    iic = ir + ic * pv.nrow;
                }
                else
                {
                    iic = ir * pv.ncol + ic;
                }
                line[tri ? j - i : j] = mat[iic];
            }
        }
    }
#else
    (void)pv;
    (void)ks_solver;
    for (int j = (tri ? i : 0); j < dim; ++j)
    {
        line[tri ? j - i : j] = mat[i * dim + j];
    }
#endif
}
} // namespace

// output a square matrix
template <typename T>
void ModuleIO::save_mat(const int istep,
    const T* mat,
    const int dim,
    const bool bit,
    const int precision,
    const bool tri,
    const bool app,
    const std::string& filename,
    const Parallel_2D& pv,
    const int drank,
    const std::string& ks_solver,
    const bool reduce)
{
    ModuleBase::TITLE("ModuleIO", "save_mat");
    ModuleBase::timer::start("ModuleIO", "save_mat");

    const bool gamma_only = std::is_same<T, double>::value;

    // write .dat file
    if (bit)
    {
        FILE* out_matrix = nullptr;
#ifdef __MPI
        if (drank == 0)
        {
#endif
            const char* mode = (app && istep > 0) ? "ab" : "wb";
            out_matrix = fopen(filename.c_str(), mode);
            if (out_matrix == nullptr)
            {
                ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
            }
            fwrite(&dim, sizeof(int), 1, out_matrix);
#ifdef __MPI
        }
#endif

        std::vector<T> line(tri ? dim : dim);
        for (int i = 0; i < dim; ++i)
        {
            const int line_len = tri ? dim - i : dim;
            line.resize(line_len);
            gather_row(mat, dim, i, tri, pv, ks_solver, line);

#ifdef __MPI
            if (reduce)
            {
                Parallel_Reduce::reduce_all(line.data(), line_len);
            }

            if (drank == 0)
            {
                for (int j = (tri ? i : 0); j < dim; ++j)
                {
                    fwrite(&line[tri ? j - i : j], sizeof(T), 1, out_matrix);
                }
            }

            MPI_Barrier(DIAG_WORLD);
#else
            for (int j = (tri ? i : 0); j < dim; ++j)
            {
                fwrite(&line[tri ? j - i : j], sizeof(T), 1, out_matrix);
            }
#endif
        }

#ifdef __MPI
        if (drank == 0)
        {
            fclose(out_matrix);
        }
#else
        fclose(out_matrix);
#endif
    } // end writing .dat file
    else // write .txt file
    {
        std::ofstream out_matrix;
        out_matrix << std::scientific << std::setprecision(precision);
#ifdef __MPI
        if (drank == 0)
        {
            if (app && istep > 0)
            {
                out_matrix.open(filename.c_str(), std::ofstream::app);
            }
            else
            {
                out_matrix.open(filename.c_str());
            }
            if (!out_matrix.is_open())
            {
                ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
            }
            out_matrix << "#------------------------------------------------------------------------" << std::endl;
            out_matrix << "# ionic step " << istep+1 << std::endl; // istep starts from 0
            out_matrix << "# filename " << filename << std::endl;
            out_matrix << "# gamma only " << gamma_only << std::endl;
            out_matrix << "# rows " << dim << std::endl;
            out_matrix << "# columns " << dim << std::endl;
            out_matrix << "#------------------------------------------------------------------------" << std::endl;

        }

        std::vector<T> line(tri ? dim : dim);
        for (int i = 0; i < dim; i++)
        {
            const int line_len = tri ? dim - i : dim;
            line.resize(line_len);
            gather_row(mat, dim, i, tri, pv, ks_solver, line);

            if (reduce)
            {
                Parallel_Reduce::reduce_all(line.data(), line_len);
            }

            if (drank == 0)
            {
                out_matrix << "Row " << i+1 << std::endl;
                size_t count = 0;
                for (int j = (tri ? i : 0); j < dim; j++)
                {
                    out_matrix << " " << line[tri ? j - i : j];
                    ++count;
                    if(count%8==0)
                    {
                        if(j!=dim-1)
                        {
                            out_matrix << std::endl;
                        }
                    }
                }
                out_matrix << std::endl;
            }
        }

        if (drank == 0)
        {
            out_matrix.close();
        }
#else
        if (app)
        {
            out_matrix.open(filename.c_str(), std::ofstream::app);
        }
        else
        {
            out_matrix.open(filename.c_str());
        }
        if (!out_matrix.is_open())
        {
            ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
        }

        out_matrix << dim;
        out_matrix << std::setprecision(precision);
        std::vector<T> line(dim);
        for (int i = 0; i < dim; i++)
        {
            gather_row(mat, dim, i, tri, pv, ks_solver, line);
            for (int j = (tri ? i : 0); j < dim; j++)
            {
                out_matrix << " " << line[tri ? j - i : j];
            }
            out_matrix << std::endl;
        }
        out_matrix.close();
#endif
    }
    ModuleBase::timer::end("ModuleIO", "save_mat");
    return;
}

// Explicit instantiations
template void ModuleIO::save_mat<double>(const int, const double*, const int, const bool,
    const int, const bool, const bool, const std::string&, const Parallel_2D&,
    const int, const std::string&, const bool);
template void ModuleIO::save_mat<std::complex<double>>(const int, const std::complex<double>*, const int,
    const bool, const int, const bool, const bool, const std::string&, const Parallel_2D&,
    const int, const std::string&, const bool);
template void ModuleIO::save_mat<float>(const int, const float*, const int, const bool,
    const int, const bool, const bool, const std::string&, const Parallel_2D&,
    const int, const std::string&, const bool);
template void ModuleIO::save_mat<std::complex<float>>(const int, const std::complex<float>*, const int,
    const bool, const int, const bool, const bool, const std::string&, const Parallel_2D&,
    const int, const std::string&, const bool);
