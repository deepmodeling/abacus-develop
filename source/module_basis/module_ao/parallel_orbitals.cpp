#include "parallel_orbitals.h"

#include "module_base/memory.h"
#include "module_basis/module_ao/ORB_control.h"
#ifdef __MPI
extern "C"
{
#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"
}
#endif
Parallel_2D::Parallel_2D()
{}
Parallel_2D::~Parallel_2D()
{
    delete[] trace_loc_row;
    delete[] trace_loc_col;
}

Parallel_Orbitals::Parallel_Orbitals()
{
    loc_sizes = nullptr;

    testpb = 0; // mohan add 2011-03-16

    // in multi-k, 2D-block-division variables for FT (R<->k)
    nnr = 1;
    nlocdim = nullptr;
    nlocstart = nullptr;
}

Parallel_Orbitals::~Parallel_Orbitals()
{
    delete[] loc_sizes;    
    delete[] nlocdim;
    delete[] nlocstart;
}

void Parallel_2D::set_proc_dim(const int& dsize, bool mode /*= 0*/)
{
    if (mode) //dim0 >= dim1
    {
        this->dim1 = (int)sqrt((double)dsize);
        while (dsize % this->dim1 != 0)
        {
            this->dim1 = this->dim1 - 1;
        }
        assert(this->dim1 > 0);
        this->dim0 = dsize / this->dim1;
    }
    else    //dim0 <= dim1
    {
        this->dim0 = (int)sqrt((double)dsize);
        while (dsize % this->dim0 != 0)
        {
            this->dim0 = this->dim0 - 1;
        }
        assert(this->dim0 > 0);
        this->dim1 = dsize / this->dim0;
    }
}

#ifdef __MPI

void Parallel_2D::mpi_create_cart()
{
    ModuleBase::TITLE("Parallel_2D", "mpi_create_cart");
    // the matrix is divided as ( dim0 * dim1 )
    int period[2] = { 1,1 };
    int dim[2] = { this->dim0, this->dim1 };
    int reorder = 0;
    MPI_Cart_create(DIAG_WORLD, 2, dim, period, reorder, &this->comm_2D);
    return;
}
#endif

bool Parallel_2D::in_this_processor(const int& iw1_all, const int& iw2_all) const
{
    if (trace_loc_row[iw1_all] == -1)
        return false;
    else if (trace_loc_col[iw2_all] == -1)
        return false;
    return true;
}

void Parallel_2D::set_global2local(const int& M_A, const int& N_A,
    bool& div_2d, std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Parallel_Orbitals", "set_global2local");

    this->init_global2local(M_A, N_A, ofs_running);
    if (!div_2d) // xiaohui add 2013-09-02
    {
        std::cout << " common settings for trace_loc_row and trace_loc_col " << std::endl;
        for (int i = 0; i < M_A; i++) this->trace_loc_row[i] = i;
        for (int i = 0; i < N_A; i++) this->trace_loc_col[i] = i;
        this->nrow = M_A;
        this->ncol = N_A;
        this->nloc = this->nrow * this->ncol;
    }
    else
    {
        // ofs_running << " nrow=" << nrow << std::endl;
        for (int irow = 0; irow < this->nrow; irow++)
        {
            int global_row = this->row_set[irow];
            this->trace_loc_row[global_row] = irow;
            // ofs_running << " global_row=" << global_row
            // << " trace_loc_row=" << this->trace_loc_row[global_row] << std::endl;
        }

        // ofs_running << " ncol=" << ncol << std::endl;
        for (int icol = 0; icol < this->ncol; icol++)
        {
            int global_col = this->col_set[icol];
            this->trace_loc_col[global_col] = icol;
            // ofs_running << " global_col=" << global_col
            // << " trace_loc_col=" << this->trace_loc_row[global_col] << std::endl;
        }
    }

    return;
}

void Parallel_2D::init_global2local(const int& M_A, const int& N_A, std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Parallel_Orbitals", "init_global2local");
    assert(M_A > 0);
    assert(N_A > 0);

    delete[] this->trace_loc_row;
    delete[] this->trace_loc_col;

    ModuleBase::GlobalFunc::OUT(ofs_running, "trace_loc_row dimension", M_A);
    ModuleBase::GlobalFunc::OUT(ofs_running, "trace_loc_col dimension", N_A);

    this->trace_loc_row = new int[M_A];
    this->trace_loc_col = new int[N_A];

    for (int i = 0; i < M_A; i++) this->trace_loc_row[i] = -1;
    for (int i = 0; i < N_A; i++) this->trace_loc_col[i] = -1;

    ModuleBase::Memory::record("trace_row_col", sizeof(int) * M_A);
    ModuleBase::Memory::record("trace_row_col", sizeof(int) * N_A);
    return;
}

#ifdef __MPI
int Parallel_2D::set_local2global(
    const int& M_A,
    const int& N_A,
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("Parallel_2D", "set_local2global");

    int dim[2];
    int period[2];
    int j, end_id, block;
    int row_b, col_b;

    // (0) every processor get it's id on the 2D comm
    // : ( coord[0], coord[1] )
    MPI_Cart_get(this->comm_2D, 2, dim, period, coord);
    assert(dim[0] == this->dim0);
    assert(dim[1] == this->dim1);

    // (1.1) how many blocks at least
    // eg. M_A = 6400, nb = 64;
    // so block = 10;
    block = M_A / nb;

    // (1.2) If data remain, add 1.
    if (block * nb < M_A)
    {
        block++;
    }

    if (this->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Total Row Blocks Number", block);

    // mohan add 2010-09-12
    if (dim[0] > block)
    {
        ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
        ofs_warning << " but, the number of row blocks is " << block << std::endl;
        if (nb > 1)
        {
            return 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("Parallel_Orbitals::set_local2global", "some processor has no row blocks, try a smaller 'nb2d' parameter.");
        }
    }

    // (2.1) row_b : how many blocks for this processor. (at least)
    row_b = block / dim[0];

    // (2.2) row_b : how many blocks in this processor.
    // if there are blocks remain, some processors add 1.
    if (coord[0] < block % dim[0])
    {
        row_b++;
    }

    if (this->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Local Row Block Number", row_b);

    // (3) end_id indicates the last block belong to
    // which processor.
    if (block % dim[0] == 0)
    {
        end_id = dim[0] - 1;
    }
    else
    {
        end_id = block % dim[0] - 1;
    }

    if (this->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Ending Row Block in processor", end_id);

    // (4) this->nrow : how many rows in this processors :
    // the one owns the last block is different.
    if (coord[0] == end_id)
    {
        this->nrow = (row_b - 1) * nb + (M_A - (block - 1) * nb);
    }
    else
    {
        this->nrow = row_b * nb;
    }

    if (this->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Local rows (including nb)", this->nrow);

    // (5) row_set, it's a global index :
    // save explicitly : every row in this processor
    // belongs to which row in the global matrix.
    this->row_set.resize(this->nrow);
    j = 0;
    for (int i = 0; i < row_b; i++)
    {
        for (int k = 0; k < nb && (coord[0] * nb + i * nb * dim[0] + k < M_A); k++, j++)
        {
            this->row_set[j] = coord[0] * nb + i * nb * dim[0] + k;
            // ofs_running << " j=" << j << " row_set=" << this->row_set[j] << std::endl;
        }
    }

    // the same procedures for columns.
    if (this->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Total Col Blocks Number", block);

    if (dim[1] > block)
    {
        ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
        ofs_warning << " but, the number of column blocks is " << block << std::endl;
        if (nb > 1)
        {
            return 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("Parallel_Orbitals::set_local2global", "some processor has no column blocks.");
        }
    }

    col_b = block / dim[1];
    if (coord[1] < block % dim[1])
    {
        col_b++;
    }

    if (this->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Local Row Block Number", col_b);

    if (block % dim[1] == 0)
    {
        end_id = dim[1] - 1;
    }
    else
    {
        end_id = block % dim[1] - 1;
    }

    if (this->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Ending Row Block in processor", end_id);

    if (coord[1] == end_id)
    {
        this->ncol = (col_b - 1) * nb + (M_A - (block - 1) * nb);
    }
    else
    {
        this->ncol = col_b * nb;
    }

    if (this->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Local columns (including nb)", this->ncol);

    //set nloc
    this->nloc = this->nrow * this->ncol;

    this->col_set.resize(this->ncol);

    j = 0;
    for (int i = 0; i < col_b; i++)
    {
        for (int k = 0; k < nb && (coord[1] * nb + i * nb * dim[1] + k < M_A); k++, j++)
        {
            this->col_set[j] = coord[1] * nb + i * nb * dim[1] + k;
        }
    }
}

int Parallel_Orbitals::set_nloc_wfc_Eij(
    const int& N_A,
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("Parallel_Orbitals", "set_local2global");
    // for wavefuncton , calculate nbands_loc
    int end_id;
    int block = N_A / nb;
    if (block * nb < N_A)
    {
        block++;
    }
    if (dim1 > block)
    {
        ofs_warning << " cpu 2D distribution : " << dim0 << "*" << dim1 << std::endl;
        ofs_warning << " but, the number of bands-row-block is " << block << std::endl;
        if (nb > 1)
        {
            return 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("Parallel_Orbitals::set_local2global", "some processor has no bands-row-blocks.");
        }
    }
    int col_b_bands = block / dim1;
    if (coord[1] < block % dim1)
    {
        col_b_bands++;
    }
    if (block % dim1 == 0)
    {
        end_id = dim1 - 1;
    }
    else
    {
        end_id = block % dim1 - 1;
    }
    if (coord[1] == end_id)
    {
        this->ncol_bands = (col_b_bands - 1) * nb + (N_A - (block - 1) * nb);
    }
    else
    {
        this->ncol_bands = col_b_bands * nb;
    }
    this->nloc_wfc = this->ncol_bands * this->nrow;

    this->nloc_Eij = this->ncol_bands * this->ncol_bands;

    return 0;
}
#endif


void Parallel_Orbitals::set_atomic_trace(const int* iat2iwt, const int &nat, const int &nlocal)
{
    this->atom_begin_col.resize(nat);
    this->atom_begin_row.resize(nat);
    for(int iat=0;iat<nat-1;iat++)
    {
        this->atom_begin_col[iat] = -1;
        this->atom_begin_row[iat] = -1;
        int irow = iat2iwt[iat];
        int icol = iat2iwt[iat];
        const int max = (iat == nat-1) ? (nlocal - irow): (iat2iwt[iat+1] - irow);
        //find the first row index of atom iat
        for(int i=0;i<max;i++)
        {
            if(this->trace_loc_row[irow]!=-1)
            {
                this->atom_begin_row[iat] = irow;
                break;
            }
            irow++;
        }
        //find the first col index of atom iat
        for(int i=0;i<max;i++)
        {
            if(this->trace_loc_col[icol]!=-1)
            {
                this->atom_begin_col[iat] = icol;
                break;
            }
            icol++;
        }
    }
}

// Get the number of columns of the parallel orbital matrix
int Parallel_Orbitals::get_col_size()const
{
    return this->ncol;
}
// Get the number of rows of the parallel orbital matrix
int Parallel_Orbitals::get_row_size()const
{
    return this->nrow;
}
// Get the number of columns of the orbital matrix of the iat-th atom
int Parallel_Orbitals::get_col_size(int iat) const
{
    int size = this->atom_begin_col[iat];
    // If the iat-th atom does not have an orbital matrix, return 0
    if(size == -1)
    {
        return 0;
    }
    iat += 1;
    // Traverse the orbital matrices of the atom and calculate the number of columns
    while(this->atom_begin_col[iat] <= this->ncol)
    {
        if(this->atom_begin_col[iat] != -1)
        {
            size = this->atom_begin_col[iat] - size;
            return size;
        }
        iat++;
    }
    // If the orbital matrix is not found after all atoms are traversed, throw an exception
    throw std::string("error in get_col_size(iat)");
}
// Get the number of rows of the orbital matrix of the iat-th atom
int Parallel_Orbitals::get_row_size(int iat) const
{
    int size = this->atom_begin_row[iat];
    if(size == -1)
    {
        return 0;
    }
    iat += 1;
    while(this->atom_begin_row[iat] <= this->ncol)
    {
        if(this->atom_begin_row[iat] != -1)
        {
            size = this->atom_begin_row[iat] - size;
            return size;
        }
        iat++;
    }
    // If the orbital matrix is not found after all atoms are traversed, throw an exception
    throw std::string("error in get_col_size(iat)");
}
