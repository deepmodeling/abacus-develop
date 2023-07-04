#include "parallel_orbitals.h"

#ifdef __MPI
extern "C"
{
#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"
}
#endif

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
            if (this->global2local_row[irow] != -1)
            {
                this->atom_begin_row[iat] = irow;
                break;
            }
            irow++;
        }
        //find the first col index of atom iat
        for(int i=0;i<max;i++)
        {
            if (this->global2local_col[icol] != -1)
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

#ifdef __MPI
void Parallel_Orbitals::set_desc_wfc_Eij(const int& nbasis, const int& nbands, const int& lld)
{
    ModuleBase::TITLE("Parallel_2D", "set_desc_wfc_Eij");
#ifdef __DEBUG
    assert(this->comm_2D != MPI_COMM_NULL);
    assert(nbasis > 0 && nbands > 0 && lld > 0);
    assert(this->nb > 0 && this->dim0 > 0 && this->dim1 > 0);
#endif
    int ISRC = 0;
    int info = 0;
    descinit_(desc_wfc, &nbasis, &nbands, &this->nb, &this->nb, &ISRC, &ISRC, &this->blacs_ctxt, &lld, &info);
    descinit_(desc_wfc1, &nbands, &nbasis, &this->nb, &this->nb, &ISRC, &ISRC, &this->blacs_ctxt, &lld, &info);
    descinit_(desc_Eij, &nbands, &nbands, &this->nb, &this->nb, &ISRC, &ISRC, &this->blacs_ctxt, &lld, &info);
}
int Parallel_Orbitals::set_nloc_wfc_Eij(
    const int& N_A,
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("Parallel_Orbitals", "set_nloc_wfc_Eij");
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
            ModuleBase::WARNING_QUIT("Parallel_Orbitals::set_nloc_wfc_Eij", "some processor has no bands-row-blocks.");
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