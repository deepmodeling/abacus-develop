#include "parallel_orbitals.h"

Parallel_Orbitals::Parallel_Orbitals()
{
    loc_sizes = nullptr;
    trace_loc_row = nullptr;
    trace_loc_col = nullptr;

    testpb = 0; // mohan add 2011-03-16
    // default value of nb is 1,
    // but can change to larger value from input.
    nb = 1;

    // in multi-k, 2D-block-division variables for FT (R<->k)
    nnr = 1;
    nlocdim = nullptr;
    nlocstart = nullptr;
}

Parallel_Orbitals::~Parallel_Orbitals()
{
    delete[] trace_loc_row;
    delete[] trace_loc_col;
    delete[] loc_sizes;    
    delete[] nlocdim;
    delete[] nlocstart;
}

bool Parallel_Orbitals::in_this_processor(const int& iw1_all, const int& iw2_all) const
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