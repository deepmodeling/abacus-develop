#ifndef ATOM_PAIR_H
#define ATOM_PAIR_H

// #include "module_cell/atom_spec.h"
#include "base_matrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"

#include <vector>

namespace hamilt
{
/**
Class: AtomPair

    the simplest way to use this class is:
    {
        AtomPair<T> atom_pair(atom_i, atom_j);
        atom_pair.set_size(row_size, col_size);
        const int rx = 0, ry = 0, rz = 0;
        auto tmp_matrix = atom_pair.get_HR_values(rx, ry, rz);
        //1. save array to tmp_matrix
        std::vector<T> local_matrix_ij = ...;
        tmp_matrix.add_array(local_matrix_ij.data());
        //2. get pointer of tmp_matrix and save array to it
        T* tmp_matrix_pointer = tmp_matrix.get_pointer();
        for(int orb_i = 0;orb_i<atom_pair.get_row_size();orb_i++)
        {
            for(int orb_j = 0;orb_j<atom_pair.get_col_size();orb_j++)
            {
                *tmp_matrix_pointer++ = ...;
            }
        }
    }
*/

template <typename T>
class AtomPair
{
  public:
    // Constructor of class AtomPair
    // Only for 2d-block MPI parallel case
    // This constructor used for initialize a atom-pair local Hamiltonian with only center cell
    // which is used for constructing HK (k space Hamiltonian) objects, (gamma_only case)
    AtomPair(const int& atom_i_,              // atomic index of atom i, used to identify atom
             const int& atom_j_,              // atomic index of atom j, used to identify atom
             const Parallel_Orbitals* paraV_, // information for 2d-block parallel
             T* existed_matrix
             = nullptr // if nullptr, new memory will be allocated, otherwise this class is a data wrapper
    );
    // Constructor of class AtomPair
    // Only for 2d-block MPI parallel case
    // This constructor used for initialize a atom-pair local Hamiltonian with non-zero cell indexes,
    // which is used for constructing HR (real space Hamiltonian) objects.
    AtomPair(const int& atom_i_,              // atomic index of atom i, used to identify atom
             const int& atom_j_,              // atomic index of atom j, used to identify atom
             const int& rx,                   // x coordinate of cell
             const int& ry,                   // y coordinate of cell
             const int& rz,                   // z coordinate of cell
             const Parallel_Orbitals* paraV_, // information for 2d-block parallel
             T* existed_array
             = nullptr // if nullptr, new memory will be allocated, otherwise this class is a data wrapper
    );
    // This constructor used for initialize a atom-pair local Hamiltonian with only center cell
    // which is used for constructing HK (k space Hamiltonian) objects, (gamma_only case)
    AtomPair(const int& atom_i,         // atomic index of atom i, used to identify atom
             const int& atom_j,         // atomic index of atom j, used to identify atom
             const int* row_atom_begin, // array, contains starting indexes in Hamiltonian matrix of atom i
             const int* col_atom_begin, // array, contains starting indexes in Hamiltonian matrix of atom j
             const int& natom,
             T* existed_matrix = nullptr);

    // This constructor used for initialize a atom-pair local Hamiltonian with non-zero cell indexes,
    // which is used for constructing HR (real space Hamiltonian) objects.
    AtomPair(const int& atom_i,         // atomic index of atom i, used to identify atom
             const int& atom_j,         // atomic index of atom j, used to identify atom
             const int& rx,             // x coordinate of cell
             const int& ry,             // y coordinate of cell
             const int& rz,             // z coordinate of cell
             const int* row_atom_begin, // array, contains starting indexes in Hamiltonian matrix of atom i
             const int* col_atom_begin, // array, contains starting indexes in Hamiltonian matrix of atom j
             const int& natom,
             T* existed_matrix = nullptr);

    // copy constructor
    AtomPair(const AtomPair<T>& other);
    // move constructor
    AtomPair(AtomPair&& other) noexcept;

    // simple constructor, only set atom_i and atom_j
    AtomPair(const int& atom_i_, // atomic index of atom i, used to identify atom
             const int& atom_j_  // atomic index of atom j, used to identify atom
    );
    // Destructor of class AtomPair
    ~AtomPair();

    /**
     * @brief allocate memory for all the BaseMatrix
    */
    void allocate(bool if_zero = false);

    /**
     * @brief set values in every BaseMatrix to zero
    */
    void set_zero();

    /**
     * @brief get col_size for this AtomPair
    */
    int get_col_size() const;
    /**
     * @brief get row_size for this AtomPair
    */
    int get_row_size() const;
    /**
     * @brief get atom_i and atom_j for this AtomPair
    */
    int get_atom_i() const;
    int get_atom_j() const;
    /**
     * @brief set col_size and row_size
    */
    void set_size(const int& col_size_in, const int& row_size_in);
    /**
     * @brief get size = col_size * row_size
     * @return int
    */
    int get_size() const;

    /**
     * @brief get Parallel_Orbitals pointer of this AtomPair for checking 2d-block parallel
     * @return const Parallel_Orbitals*
    */
    const Parallel_Orbitals* get_paraV() const;

    /// use atom_i and atom_j to identify the atom-pair
    bool identify(const AtomPair<T>& other) const;
    bool identify(const int& atom_i_, const int& atom_j_) const;

    /**
     * @brief get target BaseMatrix of target cell
     * for const AtomPair, it will return a const BaseMatrix<T> object,
     * and if not found, it will throw a error message
     * for non-const AtomPair, it will return a BaseMatrix<T> object,
     * and if not found, it will insert a new one and return it
     * 
     * @param rx_in x coordinate of cell
     * @param ry_in y coordinate of cell
     * @param rz_in z coordinate of cell
     * @return BaseMatrix<T>&
     */
    BaseMatrix<T>& get_HR_values(int rx_in, int ry_in, int rz_in);
    const BaseMatrix<T>& get_HR_values(int rx_in, int ry_in, int rz_in) const;
    
    /**
     * @brief get target BaseMatrix of index of this->values
     * it will return a BaseMatrix<T> object,
     * and if not found, it will throw a error message
     * 
     * @param index index of this->values
     * @return BaseMatrix<T>&
     */
    BaseMatrix<T>& get_HR_values(const int& index) const;

    // interface for get (rx, ry, rz) of index-th R-index in this->R_index, the return should be int[3]
    int* get_R_index(const int& index) const;
    // interface for get (rx, ry, rz) of current_R, the return should be int[3]
    int* get_R_index() const;
    // interface for search (rx, ry, rz) in this->R_index, if found, current_R would be set to index
    int find_R(const int& rx_in, const int& ry_in, const int& rz_in) const;
    // interface for search (rx, ry, rz) in this->R_index, if found, current_R would be set to index
    // and return BaseMatrix<T>* of this->values[index]
    const BaseMatrix<T>* find_matrix(const int& rx_in, const int& ry_in, const int& rz_in) const;
    BaseMatrix<T>* find_matrix(const int& rx_in, const int& ry_in, const int& rz_in);

    // this interface will call get_value in this->values
    // these four interface can be used only when R-index has been choosed (current_R >= 0)
    T& get_value(const int& i) const;
    T& get_value(const int& row, const int& col) const;
    T& get_matrix_value(const size_t& i_row_global, const size_t& j_col_global) const;
    T* get_pointer(int ir=-1) const;

    // add another BaseMatrix<T> to this->values with specific R index.
    void convert_add(const BaseMatrix<T>& target, int rx_in, int ry_in, int rz_in);
    /**
     * @brief merge another AtomPair to this AtomPair
     *
     * @param other Another AtomPair
     */
    void merge(const AtomPair<T>& other, bool skip_R = false);

    /**
     * @brief merge all values in this AtomPair to one BaseMatrix with R-index (0, 0, 0)
     * in this case, H_gamma = sum_{R} H_R will be saved in this->values[0]
     */
    void merge_to_gamma();

    /**
     * @brief Add this->value[current_R] * kphase as a block matrix of hk.
     *
     * For row major dense matrix (hk_type == 0): value[current_R][i*col_size+j] -> hk[(row_ap+i) * ld_hk + col_ap + j]
     * For column major dense matrix (hk_type == 1): value[current_R][i*col_size+j] -> hk[row_ap + i + (col_ap+j) *
     * ld_hk] For sparse matrix (hk_type == 2): not implemented yet
     *
     * @param hk Pointer to the target matrix.
     * @param ld_hk Leading dimension of the target matrix.
     * @param kphase Complex scalar to be multiplied with the block matrix.
     * @param hk_type The type of matrix layout (default: 0).
     */
    void add_to_matrix(std::complex<T>* hk,
                       const int ld_hk,
                       const std::complex<T>& kphase,
                       const int hk_type = 0) const;

    /**
     * @brief Add this->value[current_R] * kphase as a block matrix of hk.
     * for non-collinear spin case only
     */
    void add_to_matrix(T* hk, const int ld_hk, const T& kphase, const int hk_type = 0) const;

    /**
     * @brief Add this->value[current_R] * kphase to an array.
     * T = double or float
     */
    void add_to_array(std::complex<T>* target_array, const std::complex<T>& kphase) const;
    /**
     * @brief Add this->value[current_R] * kphase to an array.
     * for non-collinear spin case only (T = complex<double> or complex<float>)
     */
    void add_to_array(T* target_array, const T& kphase) const;

    // comparation function, used for sorting
    bool operator<(const AtomPair& other) const;

    // The copy assignment operator
    AtomPair& operator=(const AtomPair& other);
    // move assignment operator
    AtomPair& operator=(AtomPair&& other) noexcept;

    // interface for getting the size of this->R_index
    size_t get_R_size() const;

    /**
     * @brief get total memory size of AtomPair
    */
    size_t get_memory_size() const;

  private:
    // it contains 3 index of cell, size of R_index is three times of values.
    std::vector<int> R_index;

    // it contains containers for accessing matrix of this atom-pair
    std::vector<BaseMatrix<T>> values;

    // only for 2d-block
    const Parallel_Orbitals* paraV = nullptr;

    // the default R index is (0, 0, 0)
    // if current_R > 0, it means R index has been fixed
    // if current_R == 0, it means R index refers to the first cell
    // if current_R == 0 with gamma_only, it means R index refers to the center cell
    // !!!!!!!!!!! BE CAREFUL, current_R IS NOT THREADING-SAFE !!!!!!!!!!!!!!!!!!!!!
    mutable int current_R = 0;

    // index for identifying atom I and J for this atom-pair
    int atom_i = -1;
    int atom_j = -1;
    // start index of row for this Atom-Pair
    int row_ap = -1;
    // start index of col for this Atom-pair
    int col_ap = -1;
    int row_size = 0;
    int col_size = 0;
    int ldc = -1; // leading_dimention_colomn
};

} // namespace hamilt

#endif