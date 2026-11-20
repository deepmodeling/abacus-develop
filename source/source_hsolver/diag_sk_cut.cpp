#include "diag_sk_cut.h"

#include "source_base/matrix_block.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/tool_quit.h"
#include "source_hsolver/module_linalg/matrix_block_ops.h"
#include "source_hsolver/module_linalg/overlap_spectrum.h"
#include "source_psi/psi.h"

#ifdef __MPI
#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

namespace hsolver
{
namespace
{
void set_layout(Parallel_Blkcyclic& layout, int rows, int cols, const int* desc)
{
#ifdef __MPI
    if (desc != nullptr && desc[1] >= 0 && desc[4] > 0)
    {
        const int block = desc[4];
        const int context = desc[1];
        layout.bind_descriptor(rows, cols, block, context);
        return;
    }
#endif
    layout.set_serial(rows, cols);
}

template <typename T>
struct Matrix
{
    Parallel_Blkcyclic layout;
    std::vector<T> data;

    Matrix(int rows, int cols, const int* desc)
    {
        set_layout(layout, rows, cols, desc);
        const std::size_t storage = std::max<std::size_t>(1, layout.size());
        data.resize(storage);
    }

    ModuleBase::MatrixBlock<T> view()
    {
        const std::size_t local_rows = layout.nrloc();
        const std::size_t local_cols = layout.ncloc();
        const int* desc = nullptr;
#ifdef __MPI
        if (!layout.is_serial())
        {
            desc = layout.desc();
        }
#endif
        return {data.data(), local_rows, local_cols, desc};
    }
};

template <typename T>
void set_identity(Matrix<T>& matrix)
{
    const auto& layout = matrix.layout;
    const int rows = layout.nrloc();
    for (int col = 0; col < layout.ncloc(); ++col)
    {
        const int global_col = layout.icglb(col);
        const int row = layout.irloc(global_col);
        if (row >= 0)
        {
            matrix.data[row + col * rows] = T(1);
        }
    }
}

template <typename T>
void pack_transform(Matrix<T>& overlap_vectors, const std::vector<double>& weights,
                    int first, Matrix<T>& transform)
{
    const auto& layout = transform.layout;
    const int order = layout.nrglb();
    const int rank = layout.ncglb();
#ifdef __MPI
    if (!layout.is_serial())
    {
        const int first_column = first + 1;
        Cpxgemr2d(order, rank, overlap_vectors.data.data(), 1, first_column,
                  overlap_vectors.layout.desc(), transform.data.data(), 1, 1,
                  transform.layout.desc(), layout.ictxt());
    }
    else
#endif
    {
        const std::size_t first_entry = static_cast<std::size_t>(first) * order;
        const T* retained_columns = overlap_vectors.data.data() + first_entry;
        std::copy_n(retained_columns, layout.size(), transform.data.data());
    }
    const int ld = layout.nrloc();
    for (int col = 0; col < layout.ncloc(); ++col)
    {
        const int global_col = layout.icglb(col);
        const double scale = weights[first + global_col];
        for (int row = 0; row < ld; ++row)
        {
            transform.data[row + col * ld] *= scale;
        }
    }
}

int projected_solve_order(const Parallel_Blkcyclic& layout, const int rank)
{
    const int grid_width = std::max(layout.nrproc(), layout.ncproc());
    const int first_entry_on_last_process = (grid_width - 1) * layout.nb() + 1;
    // Require a populated grid. Distributed ELPA can also stall on a
    // compact 2x2 problem, so retain the original layout in that case.
    const bool compact_layout = rank >= first_entry_on_last_process && (rank > 2 || grid_width == 1);
    return rank == 1 || compact_layout ? rank : layout.nrglb();
}

template <typename T>
double single_projected_eigenvalue(const Matrix<T>& projected)
{
    double value = std::real(projected.data.front());
#ifdef __MPI
    const auto& layout = projected.layout;
    if (!layout.is_serial())
    {
        char scope[] = "All";
        char topology[] = " ";
        const int context = layout.ictxt();
        if (layout.in_this_processor(0, 0))
        {
            Cxgebs2d(context, scope, topology, 1, 1, &value, 1);
        }
        else
        {
            Cxgebr2d(context, scope, topology, 1, 1, &value, 1, 0, 0);
        }
    }
#endif
    return value;
}

template <typename T>
void pad_projected_hamiltonian(Matrix<T>& projected, const int rank)
{
#ifdef __MPI
    const auto& layout = projected.layout;
    const int order = layout.nrglb();
    if (order > rank)
    {
        const int one = 1;
        double work = 0.0;
        const double max_entry = ScalapackConnector::lange(
            'M', rank, rank, projected.data.data(), one, one, layout.desc(), &work);
        // Every physical eigenvalue is <= rank * max_entry. The padding
        // keeps the grid populated and never enters the requested bands.
        const double padding_energy = 2.0 * rank * std::max(1.0, max_entry);
        if (!std::isfinite(padding_energy))
        {
            ModuleBase::WARNING_QUIT("LCAO overlap", "Nonfinite energy bound for the projected Hamiltonian.");
        }
        const int rows = layout.nrloc();
        for (int col = 0; col < layout.ncloc(); ++col)
        {
            const int global_col = layout.icglb(col);
            const int row = layout.irloc(global_col);
            if (global_col >= rank && row >= 0)
            {
                projected.data[row + col * rows] = T(padding_energy);
            }
        }
    }
#endif
}

template <typename T>
void diagonalize_workspace(ModuleBase::MatrixBlock<T>& h,
                          ModuleBase::MatrixBlock<T>& s,
                          T* vectors,
                          double* eigenvalues,
                          const LcaoDiagOptions& options)
{
    const int rows = h.row;
    const int cols = h.col;
    psi::Psi<T> workspace(vectors, 1, cols, rows, rows, true);
    diag_lcao(h, s, workspace, eigenvalues, options);
}

template <typename T>
bool build_transform(const ModuleBase::MatrixBlock<T>& s,
                     double cutoff,
                     const LcaoDiagOptions& options,
                     std::vector<double>& values,
                     Matrix<T>& overlap_vectors,
                     std::vector<double>& weights)
{
    const int order = options.order;
    Matrix<T> overlap(order, order, s.desc);
    Matrix<T> identity(order, order, s.desc);
    auto overlap_view = overlap.view();
    copy_hermitian_upper(order, s, overlap_view);
    set_identity(identity);
    auto identity_view = identity.view();
    identity_view.desc = overlap_view.desc;
    auto spectrum_options = options;
    spectrum_options.bands = order;
    diagonalize_workspace(overlap_view, identity_view, overlap_vectors.data.data(), values.data(), spectrum_options);
    return overlap_inverse_sqrt(values, cutoff, weights);
}
}

template <typename T>
void diag_sk_cut(ModuleBase::MatrixBlock<T>& h,
                const ModuleBase::MatrixBlock<T>& s,
                psi::Psi<T>& psi,
                double* eigenvalues,
                const LcaoDiagOptions& options,
                int kpoint,
                double cutoff,
                std::ostream* log)
{
    const int order = options.order;
    const int bands = options.bands;
    T* vectors = psi.get_pointer();
    Matrix<T> overlap_vectors(order, order, s.desc);
    std::vector<double> values(order);
    std::vector<double> weights;
    const int displayed_k = kpoint + 1;
    auto workspace_options = options;
    workspace_options.reuse_overlap_factor = false;
    if (!build_transform(s, cutoff, workspace_options, values, overlap_vectors, weights))
    {
        const std::string message = "Invalid overlap spectrum for k point " + std::to_string(displayed_k) + ".";
        ModuleBase::WARNING_QUIT("LCAO overlap", message);
        return;
    }
    const auto retained = std::find_if(weights.begin(), weights.end(),
                                       [](double weight) { return weight > 0.0; });
    const int first = retained - weights.begin();
    const int rank = order - first;
    if (log != nullptr)
    {
        *log << " LCAO overlap k point " << displayed_k << " retained " << rank << '/' << order
             << " (min/max eigenvalue " << values.front() / values.back() << ")\n";
    }
    if (rank < bands)
    {
        const std::string message = "At k point " + std::to_string(displayed_k)
                                    + " the retained overlap rank " + std::to_string(rank)
                                    + " is smaller than nbands " + std::to_string(bands) + ".";
        ModuleBase::WARNING_QUIT("LCAO overlap", message);
        return;
    }
    if (rank == order)
    {
        // Real ELPA factors S in place. Keep the physical overlap available for
        // the next spectrum check, even when no modes were discarded this time.
        auto overlap = overlap_vectors.view();
        copy_hermitian_upper(order, s, overlap);
        overlap.desc = h.desc;
        const std::size_t local_bands = psi.get_nbands();
        if (options.method == "genelpa" && local_bands < h.col)
        {
            // GEN-ELPA writes into a square block-cyclic eigenvector buffer.
            Matrix<T> eigenvectors(order, order, h.desc);
            diagonalize_workspace(h, overlap, eigenvectors.data.data(), eigenvalues, workspace_options);
            const std::size_t elements = local_bands * h.row;
            if (elements != 0)
            {
                std::copy_n(eigenvectors.data.data(), elements, vectors);
            }
        }
        else
        {
            diag_lcao(h, overlap, psi, eigenvalues, workspace_options);
        }
        return;
    }

    Matrix<T> transform(order, rank, s.desc);
    pack_transform(overlap_vectors, weights, first, transform);
    std::vector<T>().swap(overlap_vectors.data);
    auto x = transform.view();

    const int solve_order = projected_solve_order(overlap_vectors.layout, rank);
    Matrix<T> projected(solve_order, solve_order, h.desc);
    auto h_projected = projected.view();
    {
        Matrix<T> product(order, rank, h.desc);
        auto hx = product.view();
        multiply_hermitian_upper(order, rank, h, x, hx);
        const char adjoint = std::is_same<T, double>::value ? 'T' : 'C';
        multiply_matrix_blocks(adjoint, 'N', rank, rank, order, x, hx, h_projected);
    }

    if (rank == 1)
    {
        // The 1x1 problem needs no backend eigensolver. Its eigenvector is X.
        eigenvalues[0] = single_projected_eigenvalue(projected);
        const std::size_t elements = transform.layout.size();
        if (elements != 0)
        {
            std::copy_n(transform.data.data(), elements, vectors);
        }
        return;
    }
    pad_projected_hamiltonian(projected, rank);
    Matrix<T> identity(solve_order, solve_order, h.desc);
    set_identity(identity);
    auto s_identity = identity.view();
    Matrix<T> subspace_vectors(solve_order, solve_order, h.desc);
    s_identity.desc = h_projected.desc;
    auto projected_options = workspace_options;
    projected_options.order = solve_order;
    diagonalize_workspace(h_projected, s_identity, subspace_vectors.data.data(), eigenvalues, projected_options);

    Parallel_Blkcyclic output_layout;
    set_layout(output_layout, order, bands, h.desc);
    const std::size_t rows = output_layout.nrloc();
    const std::size_t cols = output_layout.ncloc();
    const int* output_desc = nullptr;
#ifdef __MPI
    if (!output_layout.is_serial())
    {
        output_desc = output_layout.desc();
    }
#endif
    auto y = subspace_vectors.view();
    ModuleBase::MatrixBlock<T> c{vectors, rows, cols, output_desc};
    multiply_matrix_blocks('N', 'N', order, bands, rank, x, y, c);
}

template void diag_sk_cut<double>(ModuleBase::MatrixBlock<double>&,
                                  const ModuleBase::MatrixBlock<double>&,
                                  psi::Psi<double>&,
                                  double*,
                                  const LcaoDiagOptions&,
                                  int,
                                  double,
                                  std::ostream*);
template void diag_sk_cut<std::complex<double>>(ModuleBase::MatrixBlock<std::complex<double>>&,
                                                const ModuleBase::MatrixBlock<std::complex<double>>&,
                                                psi::Psi<std::complex<double>>&,
                                                double*,
                                                const LcaoDiagOptions&,
                                                int,
                                                double,
                                                std::ostream*);
} // namespace hsolver
