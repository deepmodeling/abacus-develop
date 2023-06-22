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

void Parallel_2D::set_proc_dim(const int& dsize)
{
    this->dim0 = (int)sqrt((double)dsize); // mohan update 2012/01/13
    while (dsize % this->dim0 != 0)
    {
        this->dim0 = this->dim0 - 1;
    }
    assert(this->dim0 > 0);
    this->dim1 = dsize / this->dim0;
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
