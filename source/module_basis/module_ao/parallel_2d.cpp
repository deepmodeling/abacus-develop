#include "parallel_2d.h"

#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"

#include <cassert>
#include <numeric>

bool Parallel_2D::in_this_processor(const int iw1_all, const int iw2_all) const
{
    return global2local_row(iw1_all) != -1 && global2local_col(iw2_all) != -1;
}

int Parallel_2D::get_global_row_size() const
{
#ifdef __MPI
    return desc[2];
#else
    return nrow;
#endif
}

int Parallel_2D::get_global_col_size() const
{
#ifdef __MPI
    return desc[3];
#else
    return ncol;
#endif
}

#ifdef __MPI
void Parallel_2D::_init_proc_grid(const MPI_Comm comm, const bool mode)
{
    // determine the number of rows and columns of the process grid
    // by factorizing n = p * q such that p, q are closest and p <= q
    int num_proc = 0;
    MPI_Comm_size(comm, &num_proc);
    dim0 = static_cast<int>(std::sqrt(num_proc + 0.5));
    while (dim1 = num_proc / dim0, dim0 * dim1 != num_proc)
    {
        --dim0;
    }

    if (mode)
    {
        std::swap(dim0, dim1);
    }

    // create a 2D Cartesian MPI communicator (row-major by default)
    int period[2] = {1, 1};
    int dim[2] = {dim0, dim1};
    const int reorder = 0;
    MPI_Cart_create(comm, 2, dim, period, reorder, &comm_2D);
    MPI_Cart_get(comm_2D, 2, dim, period, coord);

    // initialize the BLACS grid accordingly
    blacs_ctxt = Csys2blacs_handle(comm_2D);
    char order = 'R'; // row-major
    Cblacs_gridinit(&blacs_ctxt, &order, dim0, dim1);

    // TODO Currently MPI and BLACS are made to have the same Cartesian grid.
    // In theory, however, BLACS would split any given communicator to create
    // new ones for its own purpose when initializing the process grid, so it
    // might be unnecessary to create an MPI communicator with Cartesian topology.
    // ***This needs to be verified***
}

void Parallel_2D::_set_dist_info(const int mg, const int ng, const int nb)
{
    this->nb = nb;

    // number of local rows and columns
    const int zero = 0;
    nrow = numroc_(&mg, &nb, &coord[0], &zero, &dim0);
    ncol = numroc_(&ng, &nb, &coord[1], &zero, &dim1);
    nloc = static_cast<int64_t>(nrow) * ncol;

    // initialize the ScaLAPACK descriptor
    int info = 0, lld = std::max(nrow, 1);
    descinit_(desc, &mg, &ng, &nb, &nb, &zero, &zero, &blacs_ctxt, &lld, &info);

    // generate the global-to-local and local-to-global index maps
    local2global_row_.resize(nrow);
    global2local_row_ = std::vector<int>(mg, -1);
    for (int i = 0; i < nrow; ++i)
    {
        local2global_row_[i] = (i / nb * dim0 + coord[0]) * nb + i % nb;
        global2local_row_[local2global_row_[i]] = i;
    }

    local2global_col_.resize(ncol);
    global2local_col_ = std::vector<int>(ng, -1);
    for (int j = 0; j < ncol; ++j)
    {
        local2global_col_[j] = (j / nb * dim1 + coord[1]) * nb + j % nb;
        global2local_col_[local2global_col_[j]] = j;
    }
}

int Parallel_2D::init(const int mg, const int ng, const int nb, const MPI_Comm comm, const bool mode)
{
    _init_proc_grid(comm, mode);
    _set_dist_info(mg, ng, nb);
    return nrow == 0 || ncol == 0;
}

int Parallel_2D::set(const int mg, const int ng, const int nb, const MPI_Comm comm_2D, const int blacs_ctxt)
{
    this->comm_2D = comm_2D;
    this->blacs_ctxt = blacs_ctxt;
    Cblacs_gridinfo(blacs_ctxt, &dim0, &dim1, &coord[0], &coord[1]);
    _set_dist_info(mg, ng, nb);
    return nrow == 0 || ncol == 0;
}
#endif

void Parallel_2D::set_serial(const int mg, const int ng)
{
    assert(mg > 0 && ng > 0);

    nb = 1;
    dim0 = dim1 = 1;
    coord[0] = coord[1] = 0;
    nrow = mg;
    ncol = ng;
    nloc = nrow * ncol;
    local2global_row_.resize(nrow);
    local2global_col_.resize(ncol);
    std::iota(local2global_row_.begin(), local2global_row_.end(), 0);
    std::iota(local2global_col_.begin(), local2global_col_.end(), 0);
    global2local_row_ = local2global_row_;
    global2local_col_ = local2global_col_;
#ifdef __MPI
    comm_2D = MPI_COMM_NULL;
    blacs_ctxt = -1;
#endif
}
