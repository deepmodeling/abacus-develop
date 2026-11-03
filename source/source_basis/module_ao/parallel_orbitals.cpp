#include "parallel_orbitals.h"

#include "source_base/global_function.h"

#include <algorithm>

Parallel_Orbitals::Parallel_Orbitals()
{
    // in multi-k, 2D-block-division variables for FT (R<->k)
    this->nnr = 1;
    this->ncol_bands = 0;
    this->nrow_bands=0;
    this->nloc_wfc=0;
    this->lastband_in_proc=0;
    this->lastband_number=0;
    this->nbands = 0;

}

Parallel_Orbitals::~Parallel_Orbitals()
{
}

int Parallel_Orbitals::get_wfc_global_nbands() const
{
#ifdef __MPI
    return this->desc_wfc[3];
#else    
    return this->ncol_bands;
#endif
}

int Parallel_Orbitals::get_wfc_global_nbasis() const
{
#ifdef __MPI
    return this->desc_wfc[2];
#else
    return this->nrow_bands;
#endif
}

int Parallel_Orbitals::get_nbands() const
{
    return this->nbands;
}

void Parallel_Orbitals::set_atomic_trace(const int* iat2iwt, const int &nat, const int &nlocal)
{
    ModuleBase::TITLE("Parallel_Orbitals", "set_atomic_trace");
    this->iat2iwt_ = iat2iwt;
    int nat_plus_1 = nat + 1;
    this->atom_begin_col.resize(nat_plus_1);
    this->atom_begin_row.resize(nat_plus_1);
    for(int iat=0;iat<nat;iat++)
    {
        this->atom_begin_col[iat] = -1;
        this->atom_begin_row[iat] = -1;
        int irow = iat2iwt[iat];
        int icol = iat2iwt[iat];
        const int nw_global = (iat == nat-1) ? (nlocal - irow): (iat2iwt[iat+1] - irow);
        //find the first local row index of atom iat
        for(int i=0;i<nw_global;i++)
        {
            if (this->para_2d_.irloc(irow) != -1)
            {
                this->atom_begin_row[iat] = this->para_2d_.irloc(irow);
                break;
            }
            irow++;
        }
        //find the first local col index of atom iat
        for(int i=0;i<nw_global;i++)
        {
            if (this->para_2d_.icloc(icol) != -1)
            {
                this->atom_begin_col[iat] = this->para_2d_.icloc(icol);
                break;
            }
            icol++;
        }
    }
    this->atom_begin_row[nat] = this->para_2d_.nrloc();
    this->atom_begin_col[nat] = this->para_2d_.ncloc();
}

// Get the number of columns of the orbital matrix of the iat-th atom
int Parallel_Orbitals::get_ncol_atom(int iat) const
{
    int size = this->atom_begin_col[iat];
    // If the iat-th atom does not have an orbital matrix, return 0
    if(size == -1)
    {
        return 0;
    }
    iat += 1;
    // Traverse the orbital matrices of the atom and calculate the number of columns
    while(this->atom_begin_col[iat] <= this->para_2d_.ncloc())
    {
        if(this->atom_begin_col[iat] != -1)
        {
            size = this->atom_begin_col[iat] - size;
            return size;
        }
        iat++;
    }
    // If the orbital matrix is not found after all atoms are traversed, throw an exception
    throw std::string("error in get_ncol_atom(iat)");
}
// Get the number of rows of the orbital matrix of the iat-th atom
int Parallel_Orbitals::get_nrow_atom(int iat) const
{
    int size = this->atom_begin_row[iat];
    if(size == -1)
    {
        return 0;
    }
    iat += 1;
    while(this->atom_begin_row[iat] <= this->para_2d_.nrloc())
    {
        if(this->atom_begin_row[iat] != -1)
        {
            size = this->atom_begin_row[iat] - size;
            return size;
        }
        iat++;
    }
    // If the orbital matrix is not found after all atoms are traversed, throw an exception
    throw std::string("error in get_nrow_atom(iat)");
}

bool Parallel_Orbitals::is_invalid_atom_pair(int iat1, int iat2) const
{
    return get_nrow_atom(iat1) <= 0 || get_ncol_atom(iat2) <= 0;
}

// Get the global indexes of the rows of the parallel orbital matrix
std::vector<int> Parallel_Orbitals::get_indexes_row() const
{
    std::vector<int> indexes(this->para_2d_.nrloc());
    for(int i = 0; i < this->para_2d_.nrloc(); i++)
    {
#ifdef __MPI
        indexes[i] = this->para_2d_.irglb(i);
#else
        indexes[i] = i;
#endif
    }
    return indexes;
}
// Get the global indexes of the columns of the parallel orbital matrix
std::vector<int> Parallel_Orbitals::get_indexes_col() const
{
    std::vector<int> indexes(this->para_2d_.ncloc());
    for(int i = 0; i < this->para_2d_.ncloc(); i++)
    {
#ifdef __MPI
        indexes[i] = this->para_2d_.icglb(i);
#else
        indexes[i] = i;
#endif
    }
    return indexes;
}
// Get the global indexes of the rows of the orbital matrix of the iat-th atom
std::vector<int> Parallel_Orbitals::get_indexes_row(int iat) const
{
    int size = this->get_nrow_atom(iat);
    if(size == 0)
    {
        return std::vector<int>();
    }
    std::vector<int> indexes(size);
    int irow = this->atom_begin_row[iat];
    int begin = this->iat2iwt_[iat];
    for(int i = 0; i < size; ++i)
    {
#ifdef __MPI
        indexes[i] = this->para_2d_.irglb(irow + i) - begin;
#else
        indexes[i] = i;
#endif
    }
    return indexes;
}
// Get the global indexes of the columns of the orbital matrix of the iat-th atom
std::vector<int> Parallel_Orbitals::get_indexes_col(int iat) const
{
    int size = this->get_ncol_atom(iat);
    if(size == 0)
    {
        return std::vector<int>();
    }
    std::vector<int> indexes(size);
    int icol = this->atom_begin_col[iat];
    int begin = this->iat2iwt_[iat];
    for(int i = 0; i < size; ++i)
    {
#ifdef __MPI
        indexes[i] = this->para_2d_.icglb(icol + i) - begin;
#else
        indexes[i] = i;
#endif
    }
    return indexes;
}

#ifdef __MPI
void Parallel_Orbitals::set_wfc_and_mo_layouts(const int& nbasis, const int& nmo)
{
    ModuleBase::TITLE("BlockCyclicDistributionInformation", "set_wfc_and_mo_layouts");
#ifdef __DEBUG
    assert(nbasis > 0 && nmo > 0);
    assert(this->para_2d_.nb() > 0 && this->para_2d_.nrproc() > 0
           && this->para_2d_.ncproc() > 0);
#endif
    const int block_size = this->para_2d_.nb();
    const int blacs_context = this->para_2d_.ictxt();

    this->wfc_layout_.set(nbasis, nmo, block_size, blacs_context);
    this->mo_repr_layout_.set(nmo, nmo, block_size, blacs_context);

    std::copy(this->wfc_layout_.desc(), this->wfc_layout_.desc() + 9, this->desc_wfc);
}
int Parallel_Orbitals::set_wfc_and_mo_sizes(
    const int& nmo,
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("Parallel_Orbitals", "set_wfc_and_mo_sizes");
    // for wavefuncton , calculate nbands_loc
    this->nbands = nmo;
    int end_id = 0;
    int block = nmo / this->para_2d_.nb();
    if (block * this->para_2d_.nb() < nmo)
    {
        block++;
    }
    if (this->para_2d_.ncproc() > block)
    {
        ofs_warning << " cpu 2D distribution : " << this->para_2d_.nrproc() << "*"
                    << this->para_2d_.ncproc() << std::endl;
        ofs_warning << " but, the number of bands-row-block is " << block << std::endl;
        if (this->para_2d_.nb() > 1)
        {
            return 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("Parallel_Orbitals::set_wfc_and_mo_sizes",
                "The number of columns of the 2D process grid exceeds the number of bands. "
                "Try launching the calculation with fewer MPI processes."
            );
        }
    }
    int col_b_bands = block / this->para_2d_.ncproc();
    if (this->para_2d_.icproc() < block % this->para_2d_.ncproc())
    {
        col_b_bands++;
    }
    if (block % this->para_2d_.ncproc() == 0)
    {
        end_id = this->para_2d_.ncproc() - 1;
    }
    else
    {
        end_id = block % this->para_2d_.ncproc() - 1;
    }
    if (this->para_2d_.icproc() == end_id)
    {
        this->ncol_bands = (col_b_bands - 1) * this->para_2d_.nb()
                           + (nmo - (block - 1) * this->para_2d_.nb());
    }
    else
    {
        this->ncol_bands = col_b_bands * this->para_2d_.nb();
    }
    this->nrow_bands = this->para_2d_.nrloc();
    this->nloc_wfc = this->ncol_bands * this->para_2d_.nrloc();

    return 0;
}
#endif
