#include "atom_pair.h"
#include <complex>
#include "module_base/blas_connector.h"

namespace hamilt
{

//----------------------------------------------------
// atom pair class
//----------------------------------------------------

// destructor
template <typename T>
AtomPair<T>::~AtomPair()
{
    this->values.clear();
}

template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_, const int& atom_j_, const Parallel_Orbitals* paraV_, T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_), paraV(paraV_)
{
    assert(this->paraV != nullptr);
    this->row_ap = this->paraV->atom_begin_row[atom_i];
    this->col_ap = this->paraV->atom_begin_col[atom_j];
    if (this->row_ap == -1 || this->col_ap == -1)
    {
        throw std::string("Atom-pair not belong this process");
    }
    this->row_size = this->paraV->get_row_size(atom_i);
    this->col_size = this->paraV->get_col_size(atom_j);
    this->ldc = this->paraV->get_col_size();
    this->R_index.resize(3, 0);
    this->current_R = 0;
    if (existed_matrix != nullptr)
    {
        BaseMatrix<T> tmp(row_size, col_size, (existed_matrix + row_ap * ldc + col_ap));
        tmp.set_ldc(this->ldc);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
        this->ldc = col_size;
    }
}

template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_,
                      const int& atom_j_,
                      const int& rx,
                      const int& ry,
                      const int& rz,
                      const Parallel_Orbitals* paraV_,
                      T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_), paraV(paraV_)
{
    assert(this->paraV != nullptr);
    this->row_ap = this->paraV->atom_begin_row[atom_i];
    this->col_ap = this->paraV->atom_begin_col[atom_j];
    if (this->row_ap == -1 || this->col_ap == -1)
    {
        throw std::string("Atom-pair not belong this process");
    }
    this->row_size = this->paraV->get_row_size(atom_i);
    this->col_size = this->paraV->get_col_size(atom_j);
    this->ldc = this->paraV->get_col_size();
    this->R_index.resize(3, 0);
    this->current_R = 0;
    this->R_index[0] = rx;
    this->R_index[1] = ry;
    this->R_index[2] = rz;
    if (existed_matrix != nullptr)
    {
        BaseMatrix<T> tmp(row_size, col_size, (existed_matrix + row_ap * ldc + col_ap));
        tmp.set_ldc(this->ldc);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
        this->ldc = col_size;
    }
}
// direct save whole matrix of atom-pair
template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_,
                      const int& atom_j_,
                      const int* row_atom_begin,
                      const int* col_atom_begin,
                      const int& natom,
                      T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_)
{
    assert(row_atom_begin != nullptr && col_atom_begin != nullptr);
    this->row_ap = row_atom_begin[atom_i];
    this->col_ap = col_atom_begin[atom_j];
    this->row_size = row_atom_begin[atom_i + 1] - row_atom_begin[atom_i];
    this->col_size = col_atom_begin[atom_j + 1] - col_atom_begin[atom_j];
    this->R_index.resize(3, 0);
    this->current_R = 0;
    if (existed_matrix != nullptr)
    {
        this->ldc = row_atom_begin[natom] - row_atom_begin[0];
        BaseMatrix<T> tmp(row_size, col_size, (existed_matrix + row_ap * ldc + col_ap));
        tmp.set_ldc(this->ldc);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
        this->ldc = col_size;
    }
}
//
template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_,
                      const int& atom_j_,
                      const int& rx,
                      const int& ry,
                      const int& rz,
                      const int* row_atom_begin,
                      const int* col_atom_begin,
                      const int& natom,
                      T* existed_matrix)
    : atom_i(atom_i_), atom_j(atom_j_)
{
    assert(row_atom_begin != nullptr && col_atom_begin != nullptr);
    this->row_ap = row_atom_begin[atom_i];
    this->col_ap = col_atom_begin[atom_j];
    this->row_size = row_atom_begin[atom_i + 1] - row_atom_begin[atom_i];
    this->col_size = col_atom_begin[atom_j + 1] - col_atom_begin[atom_j];
    this->R_index.resize(3, 0);
    this->current_R = 0;
    if (existed_matrix != nullptr)
    {
        this->ldc = row_atom_begin[natom] - row_atom_begin[0];
        BaseMatrix<T> tmp(row_size, col_size, (existed_matrix + row_ap * ldc + col_ap));
        tmp.set_ldc(this->ldc);
        this->values.push_back(tmp);
    }
    else
    {
        BaseMatrix<T> tmp(row_size, col_size);
        this->values.push_back(tmp);
        this->ldc = col_size;
    }
}

template <typename T>
AtomPair<T>::AtomPair(const int& atom_i_, const int& atom_j_) : atom_i(atom_i_), atom_j(atom_j_)
{
}

// copy constructor
template <typename T>
AtomPair<T>::AtomPair(const AtomPair<T>& other)
    : R_index(other.R_index),
      values(other.values),
      paraV(other.paraV),
      current_R(other.current_R),
      atom_i(other.atom_i),
      atom_j(other.atom_j),
      row_ap(other.row_ap),
      col_ap(other.col_ap),
      row_size(other.row_size),
      col_size(other.col_size),
      ldc(other.ldc)
{
}

//allocate
template <typename T>
void AtomPair<T>::allocate(bool is_zero)
{
    for(int value=0;value<this->values.size();++value)
    {
        this->values[value].allocate(is_zero);
    }
}

// set_zero
template <typename T>
void AtomPair<T>::set_zero()
{
    for(auto& value : values)
    {
        if (value.get_pointer() != nullptr)
        {
        value.set_zero();
        }
    }
}

// The copy assignment operator
template <typename T>
AtomPair<T>& AtomPair<T>::operator=(const AtomPair<T>& other)
{
    if (this != &other)
    {
        R_index = other.R_index;
        values = other.values;
        paraV = other.paraV;
        current_R = other.current_R;
        atom_i = other.atom_i;
        atom_j = other.atom_j;
        row_ap = other.row_ap;
        col_ap = other.col_ap;
        row_size = other.row_size;
        col_size = other.col_size;
        ldc = other.ldc;
    }
    return *this;
}

// move constructor
template <typename T>
AtomPair<T>::AtomPair(AtomPair<T>&& other) noexcept
    : R_index(std::move(other.R_index)),
      values(std::move(other.values)),
      paraV(other.paraV),
      current_R(other.current_R),
      atom_i(other.atom_i),
      atom_j(other.atom_j),
      row_ap(other.row_ap),
      col_ap(other.col_ap),
      row_size(other.row_size),
      col_size(other.col_size),
      ldc(other.ldc)
{
    other.paraV = nullptr;
}

// move assignment operator
template <typename T>
AtomPair<T>& AtomPair<T>::operator=(AtomPair<T>&& other) noexcept
{
    if (this != &other)
    {
        R_index = std::move(other.R_index);
        values = std::move(other.values);
        paraV = other.paraV;
        other.paraV = nullptr;
        current_R = other.current_R;
        atom_i = other.atom_i;
        atom_j = other.atom_j;
        row_ap = other.row_ap;
        col_ap = other.col_ap;
        row_size = other.row_size;
        col_size = other.col_size;
        ldc = other.ldc;
    }
    return *this;
}

template <typename T>
bool AtomPair<T>::operator<(const AtomPair<T>& other) const
{
    if (atom_i < other.atom_i)
    {
        return true;
    }
    else if (atom_i == other.atom_i)
    {
        return atom_j < other.atom_j;
    }
    else
    {
        return false;
    }
}

// get col_size
template <typename T>
int AtomPair<T>::get_col_size() const
{
    return this->col_size;
}

// get row_size
template <typename T>
int AtomPair<T>::get_row_size() const
{
    return this->row_size;
}

// get atom_i
template <typename T>
int AtomPair<T>::get_atom_i() const
{
    return this->atom_i;
}

// get atom_j
template <typename T>
int AtomPair<T>::get_atom_j() const
{
    return this->atom_j;
}

// set size
template <typename T>
void AtomPair<T>::set_size(const int& col_size_in, const int& row_size_in)
{
    this->col_size = col_size_in;
    this->row_size = row_size_in;
}

// get size
template <typename T>
int AtomPair<T>::get_size() const
{
    return this->col_size * this->row_size;
}

// get paraV for check
template <typename T>
const Parallel_Orbitals* AtomPair<T>::get_paraV() const
{
    return this->paraV;
}

// identify
template <typename T>
bool AtomPair<T>::identify(const AtomPair<T>& other) const
{
    if (this->atom_i == other.atom_i && this->atom_j == other.atom_j)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// identify
template <typename T>
bool AtomPair<T>::identify(const int& atom_i_, const int& atom_j_) const
{
    if (this->atom_i == atom_i_ && this->atom_j == atom_j_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// get_HR_values for no-const AtomPair
template <typename T>
BaseMatrix<T>& AtomPair<T>::get_HR_values(int rx_in, int ry_in, int rz_in)
{
    // find existed R index
    if (this->find_R(rx_in, ry_in, rz_in))
    {
        // if found, return this->values[current_R]
        return this->values[current_R];
    }
    // if not found, add a new BaseMatrix for this R index
    R_index.push_back(rx_in);
    R_index.push_back(ry_in);
    R_index.push_back(rz_in);
    values.push_back(BaseMatrix<T>(this->row_size, this->col_size));
    values.back().allocate(true);
    // return the last BaseMatrix reference in values
    return this->values.back();
}

// get_HR_values for const AtomPair
template <typename T>
const BaseMatrix<T>& AtomPair<T>::get_HR_values(int rx_in, int ry_in, int rz_in) const
{
    // if current_R is -1, R index has not been fixed, find existed R index
    if (this->find_R(rx_in, ry_in, rz_in))
    {
        // if found, return this->values[current_R]
        return this->values[current_R];
    }
    // if not found, throw a error message
    else
    {
        throw std::string("AtomPair::get_HR_values: R index not found");
    }
}

// get_HR_values with index
template <typename T>
BaseMatrix<T>& AtomPair<T>::get_HR_values(const int& index) const
{
    if (index >= this->values.size())
    {
        throw std::string("AtomPair::get_HR_values: index out of range");
    }
    return const_cast<BaseMatrix<T>&>(this->values[index]);
}

// find_R
template <typename T>
bool AtomPair<T>::find_R(const int& rx_in, const int& ry_in, const int& rz_in) const
{
    for (int i = 0; i < this->R_index.size(); i += 3)
    {
        if (R_index[i] == rx_in && R_index[i + 1] == ry_in && R_index[i + 2] == rz_in)
        {
            this->current_R = i / 3;
            return true;
        }
    }
    return false;
}

// find_matrix
template <typename T>
const BaseMatrix<T>* AtomPair<T>::find_matrix(const int& rx_in, const int& ry_in, const int& rz_in) const
{
    if(!this->find_R(rx_in, ry_in, rz_in))
    {
        return nullptr;
    }
    else
    {
        return &(this->values[this->current_R]);
    }
}

// find_matrix
template <typename T>
BaseMatrix<T>* AtomPair<T>::find_matrix(const int& rx_in, const int& ry_in, const int& rz_in)
{
    if(!this->find_R(rx_in, ry_in, rz_in))
    {
        return nullptr;
    }
    else
    {
        return &(this->values[this->current_R]);
    }
}

template <typename T>
void AtomPair<T>::convert_add(const BaseMatrix<T>& target, int rx_in, int ry_in, int rz_in)
{
    BaseMatrix<T>& matrix = this->get_HR_values(rx_in, ry_in, rz_in);
    // memory type is 2d-block
    // for 2d-block memory type, the data of input matrix is expected storing in a linear array
    // so we can use pointer to access data, and will be separate to 2d-block in this function
    matrix.add_array(target.get_pointer());
}

// function merge
template <typename T>
void AtomPair<T>::merge(const AtomPair<T>& other, bool skip_R)
{
    if (other.atom_i != atom_i || other.atom_j != atom_j)
    {
        throw std::string("AtomPair::merge: atom pair not match");
    }
    int rx = 0, ry = 0, rz = 0;
    for (int i = 0; i < other.R_index.size(); i += 3)
    {
        if (!skip_R)
        {
            rx = other.R_index[i];
            ry = other.R_index[i + 1];
            rz = other.R_index[i + 2];
        }
        const BaseMatrix<T>& matrix_tmp = other.get_HR_values(i / 3);
        //if not found, push_back this BaseMatrix to this->values
        if (!this->find_R(rx, ry, rz))
        {
            this->R_index.push_back(rx);
            this->R_index.push_back(ry);
            this->R_index.push_back(rz);
            this->values.push_back(matrix_tmp);
        }
        //if found but not allocated, skip this BaseMatrix values
        else if (this->values[current_R].get_pointer() == nullptr || matrix_tmp.get_pointer() == nullptr)
        {
            continue;
        }
        // if found and allocated, add data
        else
        {
            this->values[current_R].add_array(matrix_tmp.get_pointer());
        }
    }
}

// merge_to_gamma
template <typename T>
void AtomPair<T>::merge_to_gamma()
{
    // reset R_index to (0, 0, 0)
    this->R_index.clear();
    this->R_index.resize(3, 0);
    // merge all values to first BaseMatrix
    BaseMatrix<T> tmp(this->row_size, this->col_size);
    tmp.allocate(true);
    for (int i = 0; i < this->values.size(); i++)
    {
        tmp.add_array(this->values[i].get_pointer());
    }
    this->values.clear();
    this->values.push_back(tmp);

    this->ldc = this->col_size;

    this->current_R = 0;
}

template <typename T>
void AtomPair<T>::add_to_matrix(std::complex<T>* hk,
                                const int ld_hk,
                                const std::complex<T>& kphase,
                                const int hk_type) const
{
    const BaseMatrix<T>& matrix = values[current_R];
    T* hr_tmp = matrix.get_pointer();
    std::complex<T>* hk_tmp = hk;
    T* hk_real_pointer = nullptr;
    T* hk_imag_pointer = nullptr;
    const int ld_hk_2 = ld_hk * 2;
    // row major
    if (hk_type == 0)
    {
        hk_tmp += this->row_ap * ld_hk + this->col_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            hk_real_pointer = (T*)hk_tmp;
            hk_imag_pointer = hk_real_pointer+1;
            BlasConnector::axpy(this->col_size, kphase.real(), hr_tmp, 1, hk_real_pointer, 2);
            BlasConnector::axpy(this->col_size, kphase.imag(), hr_tmp, 1, hk_imag_pointer, 2);
            hk_tmp += ld_hk;
            hr_tmp += this->ldc;
        }
    }
    // column major
    else if (hk_type == 1)
    {
        hk_tmp += this->col_ap * ld_hk + this->row_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            hk_real_pointer = (T*)hk_tmp;
            hk_imag_pointer = hk_real_pointer+1;
            BlasConnector::axpy(this->col_size, kphase.real(), hr_tmp, 1, hk_real_pointer, ld_hk_2);
            BlasConnector::axpy(this->col_size, kphase.imag(), hr_tmp, 1, hk_imag_pointer, ld_hk_2);
            hk_tmp ++;
            hr_tmp += this->ldc;
        }
    }
}

// add_to_matrix
template <typename T>
void AtomPair<T>::add_to_matrix(T* hk, const int ld_hk, const T& kphase, const int hk_type) const
{
    const BaseMatrix<T>& matrix = values[current_R];
    T* hr_tmp = matrix.get_pointer();
    T* hk_tmp = hk;
    // row major
    if (hk_type == 0)
    {
        hk_tmp += this->row_ap * ld_hk + this->col_ap;
        for (int mu = 0; mu < this->row_size; mu++)
        {
            BlasConnector::axpy(this->col_size, kphase, hr_tmp, 1, hk_tmp, 1);
            /*for (int nu = 0; nu < this->col_size; nu++)
            {
                hk_tmp[nu] += matrix.get_value(mu, nu) * kphase;
            }*/
            hk_tmp += ld_hk;
            hr_tmp += this->ldc;
        }
    }
    // column major
    else if (hk_type == 1)
    {
        hk_tmp += this->col_ap * ld_hk + this->row_ap;
        /*for (int nu = 0; nu < this->col_size; nu++)
        {
            BlasConnector::axpy(this->row_size, kphase, hr_tmp, this->ldc, hk_tmp, 1);
            for (int mu = 0; mu < this->row_size; mu++)
            {
                hk_tmp[mu] += matrix.get_value(mu, nu) * kphase;
            }
            hk_tmp += ld_hk;
            hr_tmp++;
        }*/
        for (int mu = 0; mu < this->row_size; mu++)
        {
            BlasConnector::axpy(this->col_size, kphase, hr_tmp, 1, hk_tmp, ld_hk);
            /*for (int mu = 0; mu < this->row_size; mu++)
            {
                hk_tmp[mu] += matrix.get_value(mu, nu) * kphase;
            }*/
            ++hk_tmp;
            hr_tmp += this->ldc;
        }
    }
}

// add_to_array
template <typename T>
void AtomPair<T>::add_to_array(T* array, const T& kphase) const
{
    const BaseMatrix<T>& matrix = values[current_R];
    for (int i = 0; i < this->row_size * this->col_size; i++)
    {
        array[i] += matrix.get_pointer()[i] * kphase;
    }
}

// add_to_array
template <typename T>
void AtomPair<T>::add_to_array(std::complex<T>* array, const std::complex<T>& kphase) const
{
    const BaseMatrix<T>& matrix = values[current_R];
    for (int i = 0; i < this->row_size * this->col_size; i++)
    {
        array[i] += matrix.get_pointer()[i] * kphase;
    }
}

template <typename T>
T& AtomPair<T>::get_matrix_value(const size_t& i_row_global, const size_t& j_col_global) const
{
    int i_row_local = this->paraV == nullptr ? i_row_global : this->paraV->global2local_row(i_row_global);
    int j_col_local = this->paraV == nullptr ? j_col_global : this->paraV->global2local_col(j_col_global);
#ifdef __DEBUG
    assert(i_row_local != -1 && j_col_local != -1);
    assert(current_R < this->values.size());
    assert(current_R >= 0);
#endif
    size_t i_row_in = i_row_local - row_ap;
    size_t j_col_in = j_col_local - col_ap;
    return this->values[current_R].get_value(i_row_in, j_col_in);
}

// interface for get (rx, ry, rz) of index-th R-index in this->R_index, the return should be int[3]
template <typename T>
int* AtomPair<T>::get_R_index(const int& index) const
{
    if (index >= R_index.size() / 3 || index < 0)
    {
        std::cout << "Error: index out of range in get_R_index" << std::endl;
        return nullptr;
    }
    else
    {
        // return the (int*) pointer of R_index[index*3]
        int* ptr = const_cast<int*>(&(R_index[index * 3]));
        return ptr;
    }
}

template <typename T>
int* AtomPair<T>::get_R_index() const
{
    // return the (int*) pointer of R_index[index*3]
    int* ptr = const_cast<int*>(&(R_index[current_R * 3]));
    return ptr;
}

// get_value
template <typename T>
T& AtomPair<T>::get_value(const int& i) const
{
#ifdef __DEBUG
    assert(i < this->row_size * this->col_size);
    assert(i >= 0);
    assert(current_R < this->values.size());
    assert(current_R >= 0);
#endif
    return this->values[current_R].get_pointer()[i];
}

// get_value
template <typename T>
T& AtomPair<T>::get_value(const int& row, const int& col) const
{
#ifdef __DEBUG
    assert(row < this->row_size && row >= 0);
    assert(col < this->col_size && col >= 0);
    assert(current_R < this->values.size());
    assert(current_R >= 0);
#endif
    return this->values[current_R].get_value(row, col);
}

// get_pointer
template <typename T>
T* AtomPair<T>::get_pointer(int ir) const
{
    if(ir<0)
    { 
        ir = current_R;
    }
#ifdef __DEBUG
    assert(current_R < this->values.size());
    assert(current_R >= 0);
#endif
    return this->values[ir].get_pointer();
}

// get_R_size
template <typename T>
size_t AtomPair<T>::get_R_size() const
{
#ifdef __DEBUG
    assert(this->R_index.size() / 3 == this->values.size());
    assert(this->R_index.size() % 3 == 0);
#endif
    return this->R_index.size() / 3;
}

// get_memory_size
template <typename T>
size_t AtomPair<T>::get_memory_size() const
{
    size_t memory_size = sizeof(*this);
    for (int i = 0; i < this->values.size(); i++)
    {
        memory_size += this->values[i].get_memory_size();
    }
    return memory_size;
}

// T of AtomPair can be double or complex<double>
template class AtomPair<double>;
template class AtomPair<std::complex<double>>;

} // namespace hamilt