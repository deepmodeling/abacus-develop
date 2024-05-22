#include "parallel_2d.h"
#include "module_base/memory.h"

Parallel_2D::Parallel_2D()
{}
Parallel_2D::~Parallel_2D()
{
    delete[] global2local_row_;
    delete[] global2local_col_;
}

void Parallel_2D::set_proc_dim(const int& dsize, bool mode /*= 0*/)
{
    // default mode = 0: dim0 <= dim1
    this->dim0 = (int)sqrt((double)dsize);
    while (dsize % this->dim0 != 0)
    {
        this->dim0 = this->dim0 - 1;
    }
    assert(this->dim0 > 0);
    this->dim1 = dsize / this->dim0;

    if (mode) { // mode = 1: dim0 >= dim1
        std::swap(this->dim0, this->dim1);
    }
}

bool Parallel_2D::in_this_processor(const int& iw1_all, const int& iw2_all) const
{
    return global2local_row(iw1_all) != -1 && global2local_col(iw2_all) != -1;
}

void Parallel_2D::set_global2local(const int& M_A, const int& N_A,
    const bool& div_2d, std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Parallel_2D", "set_global2local");

    this->init_global2local(M_A, N_A, ofs_running);
    if (!div_2d) // xiaohui add 2013-09-02
    {
        for (int i = 0; i < M_A; i++) this->global2local_row_[i] = i;
        for (int i = 0; i < N_A; i++) this->global2local_col_[i] = i;
        this->nrow = M_A;
        this->ncol = N_A;
        this->nloc = this->nrow * this->ncol;
    }
    else
    {
        // ofs_running << " nrow=" << nrow << std::endl;
        for (int irow = 0; irow < this->nrow; irow++)
        {
            int global_row = this->local2global_row(irow);
            this->global2local_row_[global_row] = irow;
        }

        // ofs_running << " ncol=" << ncol << std::endl;
        for (int icol = 0; icol < this->ncol; icol++)
        {
            int global_col = this->local2global_col(icol);
            this->global2local_col_[global_col] = icol;
        }
    }

    return;
}

void Parallel_2D::init_global2local(const int& M_A, const int& N_A, std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Parallel_2D", "init_global2local");
    assert(M_A > 0);
    assert(N_A > 0);

    delete[] this->global2local_row_;
    delete[] this->global2local_col_;

    ModuleBase::GlobalFunc::OUT(ofs_running, "global2local_row dimension", M_A);
    ModuleBase::GlobalFunc::OUT(ofs_running, "global2local_col dimension", N_A);

    this->global2local_row_ = new int[M_A];
    this->global2local_col_ = new int[N_A];

    for (int i = 0; i < M_A; i++) this->global2local_row_[i] = -1;
    for (int i = 0; i < N_A; i++) this->global2local_col_[i] = -1;

    ModuleBase::Memory::record("trace_row_col", sizeof(int) * M_A);
    ModuleBase::Memory::record("trace_row_col", sizeof(int) * N_A);
    return;
}

#ifdef __MPI
extern "C"
{
#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"
}

// FIXME In theory BLACS would split the given communicator to get some new
// ones for its own purpose when initializing the process grid, so there might
// be unnecessary to create a Cartesian MPI communicator in advance.
// ***This needs to be verified***

void Parallel_2D::mpi_create_cart(const MPI_Comm& diag_world)
{
    ModuleBase::TITLE("Parallel_2D", "mpi_create_cart");
#ifdef __DEBUG
    assert(this->comm_2D != MPI_COMM_NULL);
    assert(this->dim0 > 0 && this->dim1 > 0);
#endif
    // the matrix is divided as ( dim0 * dim1 )
    int period[2] = { 1,1 };
    int dim[2] = { this->dim0, this->dim1 };
    int reorder = 0;
    MPI_Cart_create(diag_world, 2, dim, period, reorder, &this->comm_2D);
    return;
}

void Parallel_2D::set_desc(const int& gr, const int& gc, const int& lld, bool first_time)
{
    ModuleBase::TITLE("Parallel_2D", "set_desc");
#ifdef __DEBUG
    assert(this->comm_2D != MPI_COMM_NULL);
    assert(gr > 0 && gc > 0 && lld > 0);
    assert(this->nb > 0 && this->dim0 > 0 && this->dim1 > 0);
#endif
    if (first_time)
    {
        int myprow=0, mypcol=0;
        char order = 'R'; // row major process grid

        blacs_ctxt = Csys2blacs_handle(comm_2D);
        Cblacs_gridinit(&blacs_ctxt, &order, dim0, dim1);
        Cblacs_gridinfo(blacs_ctxt, &dim0, &dim1, &myprow, &mypcol);
    }
    int ISRC = 0;
    int info = 0;
    descinit_(desc, &gr, &gc, &this->nb, &this->nb, &ISRC, &ISRC, &this->blacs_ctxt, &lld, &info);
}

int Parallel_2D::set_local2global(
    const int& M_A,
    const int& N_A,
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("Parallel_2D", "set_local2global");
#ifdef __DEBUG
    assert(M_A > 0 && N_A > 0);
    assert(this->nb > 0);
#endif

    int dim[2];
    int period[2];

    // (0) every processor get it's id on the 2D comm
    // : ( coord[0], coord[1] )
    MPI_Cart_get(this->comm_2D, 2, dim, period, coord);
    assert(dim[0] == this->dim0);
    assert(dim[1] == this->dim1);

    // local number of row and columns
    const int zero = 0;
    nrow = numroc_(&M_A, &nb, &coord[0], &zero, &dim0);
    ncol = numroc_(&N_A, &nb, &coord[1], &zero, &dim1);
    nloc = nrow * ncol;

    // mohan add 2010-09-12
    if (nrow == 0 || ncol == 0)
    {
        ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
        ofs_warning << " but, the number of row and column blocks are "
                    << M_A / nb + (M_A % nb != 0) << "*" << N_A / nb + (N_A % nb != 0)
                    << std::endl;
        if (nb > 1)
        {
            return 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("Parallel_2D::set_local2global", "some processor has no blocks, try a smaller 'nb2d' parameter or reduce the number of mpi processes.");
        }
    }

    local2global_row_.resize(nrow);
    for (int i = 0; i < nrow; ++i)
    {
        local2global_row_[i] = (i / nb) * dim0 * nb + coord[0] * nb + i % nb;
    }

    local2global_col_.resize(ncol);
    for (int j = 0; j < ncol; ++j)
    {
        local2global_col_[j] = (j / nb) * dim1 * nb + coord[1] * nb + j % nb;
    }

    return 0;
}
#else
void Parallel_2D::set_serial(const int& M_A, const int& N_A)
{
    ModuleBase::TITLE("Parallel_2D", "set_serial");
    this->nrow = M_A;
    this->ncol = N_A;
    this->nloc = this->nrow * this->ncol;
    this->local2global_row_.resize(this->nrow);
    this->local2global_col_.resize(this->ncol);
    for (int i = 0; i < this->nrow; i++) this->local2global_row_[i] = i;
    for (int i = 0; i < this->ncol; i++) this->local2global_col_[i] = i;
}
#endif
