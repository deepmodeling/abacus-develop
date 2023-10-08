#include "./transfer.h"

#include <set>

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#ifdef __MPI
#include <mpi.h>

#include <algorithm>

namespace hamilt
{

// ------------------------------------------------
// HTransPara
// ------------------------------------------------

template <typename T>
HTransPara<T>::HTransPara(int n_processes, HContainer<T>* hr_in)
{
    this->hr = hr_in;
    this->ap_indexes.resize(n_processes);
    this->size_values.resize(n_processes);
    this->paraV = hr_in->get_atom_pair(0).get_paraV();
    this->atom_i_index.resize(n_processes);
}

template <typename T>
HTransPara<T>::~HTransPara()
{
}

// cal_orb_indexes
template <typename T>
void HTransPara<T>::cal_orb_indexes(int irank, std::vector<int>* orb_indexes)
{
    // calculate the size of orb_indexes
    int size_orb_indexes = 1;
#ifdef __DEBUG
    assert(orb_indexes != nullptr);
    assert(this->atom_i_index.size() > 0);
#endif
    // loop atoms and cal total size
    for (int i = 0; i < this->atom_i_index[irank].size(); ++i)
    {
        int atom = this->atom_i_index[irank][i];
        size_orb_indexes += 3;
        size_orb_indexes += this->paraV->get_row_size(atom);
        size_orb_indexes += this->paraV->get_col_size(atom);
    }
    orb_indexes->resize(size_orb_indexes);
    int* data = orb_indexes->data();
    // size of atom
    *data++ = this->atom_i_index[irank].size();
    for (int i = 0; i < this->atom_i_index[irank].size(); ++i)
    {
        int atom = this->atom_i_index[irank][i];
        // atom index
        *data++ = atom;
        // size of row for this atom
        *data = this->paraV->get_row_size(atom);
        if (*data++ > 0)
        {
            // indexes of row for this atom
            auto row_indexes = this->paraV->get_indexes_row(atom);
            for (int k = 0; k < row_indexes.size(); k++)
            {
                *data++ = row_indexes[k];
            }
        }
        // size of col for this atom
        *data = this->paraV->get_col_size(atom);
        if (*data++ > 0)
        {
            // indexes of col for this atom
            auto col_indexes = this->paraV->get_indexes_col(atom);
            for (int k = 0; k < col_indexes.size(); k++)
            {
                *data++ = col_indexes[k];
            }
        }
    }
#ifdef __DEBUG
    assert(data - orb_indexes->data() == size_orb_indexes);
#endif
    return;
}

// receive_ap_indexes
template <typename T>
void HTransPara<T>::receive_ap_indexes(int irank, const int* ap_indexes_in, const long& size_ap_indexes_in)
{
    // sender and receiver are not same process
    if (ap_indexes_in == nullptr)
    {
        long size_ap_indexes = 0;
        MPI_Recv(&size_ap_indexes, 1, MPI_LONG, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        this->ap_indexes[irank].resize(size_ap_indexes);
        MPI_Recv(this->ap_indexes[irank].data(), size_ap_indexes, MPI_INT, irank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // sender and receiver are same process
    else
    {
        this->ap_indexes[irank].assign(ap_indexes_in, ap_indexes_in + size_ap_indexes_in);
    }
    // calculate the size of values
    this->size_values[irank] = 0;
    const int number_atom = this->ap_indexes[irank][0];
    const int* ap_data = this->ap_indexes[irank].data() + 1;
    std::set<int> atom_set;
    for (int i = 0; i < number_atom; ++i)
    {
        const int atom_i = *ap_data++;
        atom_set.insert(atom_i);
        const int number_atom_j = *ap_data++;
        const int size_row = this->paraV->get_row_size(atom_i);
        for (int j = 0; j < number_atom_j; ++j)
        {
            const int atom_j = *ap_data++;
            atom_set.insert(atom_j);
            const int size_col = this->paraV->get_col_size(atom_j);
            const int number_R = *ap_data++;
            if (size_row > 0 && size_col > 0)
            {
                this->size_values[irank] += size_row * size_col * number_R;
            }
            ap_data += number_R * 3;
        }
    }
    // revert atom_set to atom_i_index[irank]
    this->atom_i_index[irank].resize(atom_set.size());
    std::copy(atom_set.begin(), atom_set.end(), this->atom_i_index[irank].begin());
#ifdef __DEBUG
    assert(ap_data - this->ap_indexes[irank].data() == this->ap_indexes[irank].size());
#endif
    return;
}

// send_orb_indexes
template <typename T>
void HTransPara<T>::send_orb_indexes(int irank, MPI_Request* request)
{
    std::vector<int> orb_indexes;
    this->cal_orb_indexes(irank, &orb_indexes);
    long size_orb_indexes = orb_indexes.size();
    if (request != nullptr)
    {
        MPI_Isend(&size_orb_indexes, 1, MPI_LONG, irank, 0, MPI_COMM_WORLD, request);
        MPI_Isend(orb_indexes.data(), orb_indexes.size(), MPI_INT, irank, 1, MPI_COMM_WORLD, request);
    }
    else
    {
        MPI_Send(&size_orb_indexes, 1, MPI_LONG, irank, 0, MPI_COMM_WORLD);
        MPI_Send(orb_indexes.data(), orb_indexes.size(), MPI_INT, irank, 1, MPI_COMM_WORLD);
    }
}

template <typename T>
void HTransPara<T>::send_data(int irank, MPI_Request* request)
{
    std::vector<T> values(this->size_values[irank]);
    this->pack_data(irank, values.data());
    if (request != nullptr)
    {
        MPI_Isend(values.data(), values.size(), MPITraits<T>::datatype(), irank, 2, MPI_COMM_WORLD, request);
    }
    else
    {
        MPI_Send(values.data(), values.size(), MPITraits<T>::datatype(), irank, 2, MPI_COMM_WORLD);
    }
}

template <typename T>
void HTransPara<T>::receive_data(int irank, const T* values)
{
    // sender and receiver are not same process
    if (values == nullptr)
    {
        std::vector<T> values_tmp(this->size_values[irank]);
        MPI_Recv(values_tmp.data(),
                 values_tmp.size(),
                 MPITraits<T>::datatype(),
                 irank,
                 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        this->unpack_data(irank, values_tmp.data());
    }
    else
    {
        this->unpack_data(irank, values);
    }
}

template <typename T>
void HTransPara<T>::pack_data(int irank, T* values)
{
#ifdef __DEBUG
    assert(values != nullptr);
    assert(this->size_values[irank] != 0);
#endif
    const int number_atom = this->ap_indexes[irank][0];
    const int* ap_data = this->ap_indexes[irank].data() + 1;

    T* value_data = values;
    for (int i = 0; i < number_atom; ++i)
    {
        const int atom_i = *ap_data++;
        const int number_atom_j = *ap_data++;
        const int size_row = this->paraV->get_row_size(atom_i);
        for (int j = 0; j < number_atom_j; ++j)
        {
            const int atom_j = *ap_data++;
            const int size_col = this->paraV->get_col_size(atom_j);
            const int number_R = *ap_data++;
            for (int k = 0; k < number_R; ++k)
            {
                int r_index[3];
                r_index[0] = *ap_data++;
                r_index[1] = *ap_data++;
                r_index[2] = *ap_data++;
                if (size_row > 0 && size_col > 0)
                {
                    const T* matrix_pointer = this->hr->data(atom_i, atom_j, r_index);
                    ModuleBase::GlobalFunc::COPYARRAY(matrix_pointer, value_data, size_row * size_col);
                    value_data += size_row * size_col;
                }
            }
        }
    }
#ifdef __DEBUG
    assert(value_data - values == this->size_values[irank]);
#endif
    return;
}

template <typename T>
void HTransPara<T>::unpack_data(int irank, const T* values)
{
#ifdef __DEBUG
    assert(values != nullptr);
    assert(this->size_values[irank] != 0);
#endif
    const T alpha(1.0);
    const int number_atom = this->ap_indexes[irank][0];
    // loop AtomPairs and unpack values
    int* ap_data = this->ap_indexes[irank].data() + 1;

    const T* value_data = values;

    for (int i = 0; i < number_atom; ++i)
    {
        const int atom_i = *ap_data++;
        const int size_row = this->paraV->get_row_size(atom_i);
        const int size_j = *ap_data++;
        for (int j = 0; j < size_j; ++j)
        {
            const int atom_j = *ap_data++;
            const int size_col = this->paraV->get_col_size(atom_j);
            const int number_R = *ap_data++;
            for (int k = 0; k < number_R; ++k)
            {
                int r_index[3];
                r_index[0] = *ap_data++;
                r_index[1] = *ap_data++;
                r_index[2] = *ap_data++;
                if (size_row > 0 && size_col > 0)
                {
                    T* matrix_pointer = this->hr->data(atom_i, atom_j, r_index);
                    BlasConnector::axpy(size_row * size_col, alpha, value_data, 1, matrix_pointer, 1);
                    value_data += size_row * size_col;
                }
            }
        }
    }
#ifdef __DEBUG
    //assert(value_data - values == this->data_size[irank]);
    assert(ap_data - this->ap_indexes[irank].data() == this->ap_indexes[irank].size());
#endif
}

template <typename T>
long HTransPara<T>::get_max_size() const
{
    return *std::max_element(this->size_values.begin(), this->size_values.end());
}

template <typename T>
void HTransPara<T>::get_value_size(int* out) const
{
    for (int i = 0; i < this->size_values.size(); ++i)
    {
        out[i] = this->size_values[i];
    }
    return;
}

// ------------------------------------------------
// HTransSerial
// ------------------------------------------------

template <typename T>
HTransSerial<T>::HTransSerial(int n_processes, HContainer<T>* hr_in)
{
    this->hr = hr_in;
    this->orb_indexes.resize(n_processes);
    this->size_values.resize(n_processes);
    this->orb_col_indexes.resize(n_processes);
    this->orb_row_indexes.resize(n_processes);
}

template <typename T>
HTransSerial<T>::~HTransSerial()
{
}

template <typename T>
void HTransSerial<T>::cal_ap_indexes(int irank, std::vector<int>* ap_indexes)
{
    // calculate the size of ap_indexes
    long size_ap_indexes = this->hr->size_atom_pairs() * 2; // count of atom_j and size_r
    for (int i = 0; i < this->hr->size_atom_pairs(); i++)
    {
        size_ap_indexes += this->hr->get_atom_pair(i).get_R_size() * 3; // count of rx, ry, rz
    }
    auto& sparse_ap = this->hr->get_sparse_ap();
    int size_atom = 0;
    for (int i = 0; i < sparse_ap.size(); i++)
    {
        if (sparse_ap[i].size() > 0)
        {
            size_atom++;
        }
    }
    size_ap_indexes += size_atom * 2 + 1; // count of atom_i and size_j and size_i
    ap_indexes->resize(size_ap_indexes);
    int* data = ap_indexes->data();
    // size of atom
    *data++ = size_atom;
    auto& sparse_ap_index = this->hr->get_sparse_ap_index();
    for (int atom = 0; atom < sparse_ap.size(); ++atom)
    {
        if (sparse_ap[atom].size() > 0)
        {
            // atom index
            *data++ = atom;
            // size of atom_j
            *data++ = sparse_ap[atom].size();
            // loop of atom_j
            for (int j = 0; j < sparse_ap[atom].size(); j++)
            {
                // atom_j index
                *data++ = sparse_ap[atom][j];
                hamilt::AtomPair<T>& atom_pair = this->hr->get_atom_pair(sparse_ap_index[atom][j]);
                // size of R
                *data++ = atom_pair.get_R_size();
                // loop of R
                for (int k = 0; k < atom_pair.get_R_size(); k++)
                {
                    // rx
                    *data++ = atom_pair.get_R_index(k)[0];
                    // ry
                    *data++ = atom_pair.get_R_index(k)[1];
                    // rz
                    *data++ = atom_pair.get_R_index(k)[2];
                }
            }
        }
    }
#ifdef __DEBUG
    assert(data - ap_indexes->data() == size_ap_indexes);
#endif
    // size of atom
    return;
}

template <typename T>
void HTransSerial<T>::send_ap_indexes(int irank, MPI_Request* request)
{
    std::vector<int> ap_indexes;
    this->cal_ap_indexes(irank, &ap_indexes);
    long size_ap_indexes = ap_indexes.size();
    if (request != nullptr)
    {
        MPI_Isend(&size_ap_indexes, 1, MPI_LONG, irank, 0, MPI_COMM_WORLD, request);
        MPI_Isend(ap_indexes.data(), ap_indexes.size(), MPI_INT, irank, 1, MPI_COMM_WORLD, request);
    }
    else
    {
        MPI_Send(&size_ap_indexes, 1, MPI_LONG, irank, 0, MPI_COMM_WORLD);
        MPI_Send(ap_indexes.data(), ap_indexes.size(), MPI_INT, irank, 1, MPI_COMM_WORLD);
    }
}

// receive_orb_indexes
template <typename T>
void HTransSerial<T>::receive_orb_indexes(int irank, const int* orb_indexes_in, const long& size_orb_indexes_in)
{
    // sender and receiver are not same process
    if (orb_indexes_in == nullptr)
    {
        long size_orb_indexes = 0;
        MPI_Recv(&size_orb_indexes, 1, MPI_LONG, irank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        this->orb_indexes[irank].resize(size_orb_indexes);
        MPI_Recv(this->orb_indexes[irank].data(),
                 size_orb_indexes,
                 MPI_INT,
                 irank,
                 1,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
    // sender and receiver are same process
    else
    {
        this->orb_indexes[irank].assign(orb_indexes_in, orb_indexes_in + size_orb_indexes_in);
    }
    // calculate the index of row and col for each atom
    this->size_values[irank] = 0;
    const int number_atom = this->orb_indexes[irank][0];
    const int* orb_data = this->orb_indexes[irank].data() + 1;
    for (int i = 0; i < number_atom; ++i)
    {
        const int atom_i = *orb_data++;
        this->orb_row_indexes[irank][atom_i] = orb_data - orb_indexes[irank].data();
        orb_data += *orb_data + 1;
        this->orb_col_indexes[irank][atom_i] = orb_data - orb_indexes[irank].data();
        orb_data += *orb_data + 1;
    }
#ifdef __DEBUG
    assert(orb_data - this->orb_indexes[irank].data() == this->orb_indexes[irank].size());
#endif
    // calculate the size of values
    for (int iap = 0; iap < this->hr->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<T>& atom_pair = this->hr->get_atom_pair(iap);
        const int atom_i = atom_pair.get_atom_i();
        const int atom_j = atom_pair.get_atom_j();
        const int size_row = this->orb_indexes[irank][this->orb_row_indexes[irank][atom_i]];
        const int size_col = this->orb_indexes[irank][this->orb_col_indexes[irank][atom_j]];
        const int number_R = atom_pair.get_R_size();
        if (size_row > 0 && size_col > 0)
        {
            this->size_values[irank] += size_row * size_col * number_R;
        }
    }
}

template <typename T>
void HTransSerial<T>::send_data(int irank, MPI_Request* request)
{
    std::vector<T> values(this->size_values[irank]);
    this->pack_data(irank, values.data());
    if (request != nullptr)
    {
        MPI_Isend(values.data(), values.size(), MPITraits<T>::datatype(), irank, 2, MPI_COMM_WORLD, request);
    }
    else
    {
        MPI_Send(values.data(), values.size(), MPITraits<T>::datatype(), irank, 2, MPI_COMM_WORLD);
    }
}

template <typename T>
void HTransSerial<T>::receive_data(int irank, const T* values)
{
    // sender and receiver are not same process
    if (values == nullptr)
    {
        std::vector<T> values_tmp(this->size_values[irank]);
        MPI_Recv(values_tmp.data(),
                 values_tmp.size(),
                 MPITraits<T>::datatype(),
                 irank,
                 2,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        this->unpack_data(irank, values_tmp.data());
    }
    else
    {
        this->unpack_data(irank, values);
    }
}

template <typename T>
void HTransSerial<T>::pack_data(int irank, T* values)
{
#ifdef __DEBUG
    assert(values != nullptr);
    assert(this->size_values[irank] != 0);
#endif
    auto& sparse_ap = this->hr->get_sparse_ap();
    auto& sparse_ap_index = this->hr->get_sparse_ap_index();

#ifdef _OPENMP
    // calculate the index of each atom
    std::vector<T*> value_atoms(sparse_ap.size(), values);
    long size_begin = 0;
    for (int i = 0; i < sparse_ap.size(); ++i)
    {
        value_atoms[i] += size_begin;
        const int atom_i = i;
        if(sparse_ap[i].size() == 0) continue;
        const int size_row = this->orb_indexes[irank][this->orb_row_indexes[irank][atom_i]];
        for (int j = 0; j < sparse_ap[i].size(); ++j)
        {
            const int atom_j = sparse_ap[i][j];
            const int size_col = this->orb_indexes[irank][this->orb_col_indexes[irank][atom_j]];
            const int number_R = this->hr->get_atom_pair(sparse_ap_index[i][j]).get_R_size();
            size_begin += size_row * size_col * number_R;
        }
    }

#pragma omp parallel for
#else
    T* value_data = values;
#endif
    for (int i = 0; i < sparse_ap.size(); ++i)
    {
#ifdef _OPENMP
        T* value_data = value_atoms[i];
#endif
        if (sparse_ap[i].size() == 0)
        {
            continue;
        }
        const int atom_i = i;
        const int size_row = this->orb_indexes[irank][this->orb_row_indexes[irank][atom_i]];
        if (size_row == 0)
        {
            continue;
        }
        const int* row_index = this->orb_indexes[irank].data() + this->orb_row_indexes[irank][atom_i] + 1;
        for (int j = 0; j < sparse_ap[i].size(); ++j)
        {
            const int atom_j = sparse_ap[i][j];
            const int size_col = this->orb_indexes[irank][this->orb_col_indexes[irank][atom_j]];
            if (size_col == 0)
            {
                continue;
            }
            const int* col_index = this->orb_indexes[irank].data() + this->orb_col_indexes[irank][atom_j] + 1;
            const hamilt::AtomPair<T>& tmp_ap = this->hr->get_atom_pair(sparse_ap_index[i][j]);
            const int number_R = tmp_ap.get_R_size();
            for (int k = 0; k < number_R; ++k)
            {
                const hamilt::BaseMatrix<T>& matrix = tmp_ap.get_HR_values(k);
                for (int irow = 0; irow < size_row; ++irow)
                {
                    const int mu = row_index[irow];
                    for (int icol = 0; icol < size_col; ++icol)
                    {
                        const int nu = col_index[icol];
                        *value_data++ = matrix.get_value(mu, nu);
                    }
                }
            }
        }
    }
#ifdef __DEBUG
    //assert(value_data - values == this->size_values[irank]);
#endif
    return;
}

template <typename T>
void HTransSerial<T>::unpack_data(int irank, const T* values)
{
#ifdef __DEBUG
    assert(values != nullptr);
    assert(this->size_values[irank] != 0);
#endif
    auto& sparse_ap = this->hr->get_sparse_ap();
    auto& sparse_ap_index = this->hr->get_sparse_ap_index();

#ifdef _OPENMP
    // calculate the index of each atom
    std::vector<const T*> value_atoms(sparse_ap.size(), values);
    long size_begin = 0;
    for (int i = 0; i < sparse_ap.size(); ++i)
    {
        value_atoms[i] += size_begin;
        const int atom_i = i;
        if(sparse_ap[i].size() == 0) continue;
        const int size_row = this->orb_indexes[irank][this->orb_row_indexes[irank][atom_i]];
        for (int j = 0; j < sparse_ap[i].size(); ++j)
        {
            const int atom_j = sparse_ap[i][j];
            const int size_col = this->orb_indexes[irank][this->orb_col_indexes[irank][atom_j]];
            const int number_R = this->hr->get_atom_pair(sparse_ap_index[i][j]).get_R_size();
            size_begin += size_row * size_col * number_R;
        }
    }

#pragma omp parallel for
#else
    const T* value_data = values;
#endif
    for (int i = 0; i < sparse_ap.size(); ++i)
    {
#ifdef _OPENMP
        const T* value_data = value_atoms[i];
#endif
        if (sparse_ap[i].size() == 0)
        {
            continue;
        }
        const int atom_i = i;
        const int size_row = this->orb_indexes[irank][this->orb_row_indexes[irank][atom_i]];
        if (size_row == 0)
        {
            continue;
        }
        const int* row_index = this->orb_indexes[irank].data() + this->orb_row_indexes[irank][atom_i] + 1;
        for (int j = 0; j < sparse_ap[i].size(); ++j)
        {
            const int atom_j = sparse_ap[i][j];
            const int size_col = this->orb_indexes[irank][this->orb_col_indexes[irank][atom_j]];
            if (size_col == 0)
            {
                continue;
            }
            const int* col_index = this->orb_indexes[irank].data() + this->orb_col_indexes[irank][atom_j] + 1;
            const hamilt::AtomPair<T>& tmp_ap = this->hr->get_atom_pair(sparse_ap_index[i][j]);
            const int number_R = tmp_ap.get_R_size();
            for (int k = 0; k < number_R; ++k)
            {
                const hamilt::BaseMatrix<T>& matrix = tmp_ap.get_HR_values(k);
                for (int irow = 0; irow < size_row; ++irow)
                {
                    const int mu = row_index[irow];
                    for (int icol = 0; icol < size_col; ++icol)
                    {
                        const int nu = col_index[icol];
                        matrix.get_value(mu, nu) = *value_data++;
                    }
                }
            }
        }
    }
#ifdef __DEBUG
    //assert(value_data - values == this->size_values[irank]);
#endif
    return;
}

template <typename T>
long HTransSerial<T>::get_max_size() const
{
    return *std::max_element(this->size_values.begin(), this->size_values.end());
}
template <typename T>
void HTransSerial<T>::get_value_size(int* out) const
{
    for (int i = 0; i < this->size_values.size(); ++i)
    {
        out[i] = this->size_values[i];
    }
    return;
}

template class HTransPara<double>;
template class HTransPara<std::complex<double>>;
template class HTransSerial<double>;
template class HTransSerial<std::complex<double>>;

} // end namespace hamilt

#endif