#ifndef HCONTAINER_H
#define HCONTAINER_H

#include <vector>
#include <set>

#include "atom_pair.h"
#include "module_cell/unitcell.h"

namespace hamilt
{

/**
 * class HContainer
 * used to store a matrix for atom-pair local Hamiltonian with specific R-index
 * ----------------------------------------
 *   <Psi_{mu_I,R}|H|Psi_{nu_J,0}>
 * ----------------------------------------
 * template T can be double or complex<double>
 * 
 * examples for using this class:
 * 1. initialize a HContainer
 *    a. use unitcell to initialize atom_pairs
 *     ```
 *       // ucell is a UnitCell object
 *       // in this method, all empty atom-pairs will be initialized in HR
 *       // atom-pairs are sorted by matrix of (i, j)
 *       HContainer<double> HR(ucell);
 *     ```
 *    b. use insert_pair() to insert atom-pair
 *     ```
 *       // HR is a HContainer object, from simple constructor
 *       HContainer<double> HR;
 *       // there are many methods to init AtomPair, this case is step by step method
 *       AtomPair<double> atom_ij(0, 1);
 *       // set size of AtomPair by set_size() 
 *       atom_ij.set_size(2, 2);
 *       // allocate local matrix memory with R-index in AtomPair
 *       auto local_matrix = atom_ij.get_HR_values(0, 0, 0);
 *       // save array to local_matrix
 *       std::vector<double> local_matrix_ij = ...;
 *       local_matrix.add_array(local_matrix_ij.data());
 *       // insert atom_ij into HR
 *       HR.insert_pair(atom_ij);
 *     ```
 *   c. use Parallel_Orbital to initialize atom_pairs and HContainer
 *     ```
 *       // paraV is a Parallel_Orbital object, 2d-block-recycle parallel case
 *       HCotainer<double> HR(paraV);
 *       // initialize a new AtomPair with atom paraV
 *       AtomPair<double> atom_ij(0, 1, paraV);
 *       // insert atom_ij into HR
 *       HR.insert_pair(atom_ij);
 *     ```
 *   d. use data_pointer of existed matrix to initialize HContainer, which means HContainer is a wrapper
 *     ```
 *       // data_pointer is a pointer of existed matrix
 *       HContainer<double> HR(paraV, data_pointer);
 *       // initialize a new AtomPair with data_pointer
 *       AtomPair<double> atom_ij(0, 1, paraV, data_pointer);
 *       // insert atom_ij into HR
 *       HR.insert_pair(atom_ij);
 *     ```
 * 2. get target AtomPair with index of atom I and J, or with index in atom_pairs
 *    a. use interface find_pair() to get pointer of target AtomPair
 *     ```
 *       // HR is a HContainer object
 *       AtomPair<double>* atom_ij = HR.find_pair(0, 1);
 *       // check if atom_ij is found
 *       if (atom_ij != nullptr)
 *       {
 *          // do something
 *       }
 *     ```
 *    b. use interface get_atom_pair() to get reference of target AtomPair
 *     ```
 *       // HR is a HContainer object
 *       AtomPair<double>& atom_ij_1 = HR.get_atom_pair(0, 1);
 *       AtomPair<double>& atom_ij_2 = HR.get_atom_pair(1);//suppose 0,1 is the second atom-pair in atom_pairs
 *       // check if atom_ij_1 and atom_ij_2 are the same
 *       if(atom_ij_1.identify(atom_ij_2)) 
 *       {
 *           std::cout<<"atom_ij_1 and atom_ij_2 are the same"<<std::endl;
 *       }
 *     ```
 * 3. get data pointer of target local matrix <Psi_{mu_I,R}|H|Psi_{nu_J,0}>
 *    a. use interface data() with atom_i and atom_j and R index
 *     ```
 *       // HR is a HContainer object
 *       // suppose atom_i = 0, atom_j = 1, int[3] target_R = {0, 0, 0}
 *       double* target_data = HR.data(0, 1, target_R);
 *     ```
 *    b. fix_R and use data() with atom_i and atom_j without R index
 *     ```
 *       HR.fix_R(0, 0, 0);
 *       double* target_data = HR.data(0, 1);
 *       HR.unfix_R();
 *     ```
 * 4. use for-loop to do something with atom-pairs with specific R index
 *    a. loop R-index first and then loop atom-pairs
 *     ```
 *       // HR is a const HContainer object, which has been initialized
 *       int rx, ry, rz;
 *       // call size_R_loop() to get number of different R indexes, 
 *       // be careful, it will cost some time to loop all atom-pairs to gather R indexes
 *       int size_for_loop_R = HR.size_R_loop();
 *       for (int iR = 0; iR < size_for_loop_R ; iR++)
 *       {
 *           // call loop_R() to get R coordinates (rx, ry, rz)
 *           HR.loop_R(iR, rx, ry, rz);
 *           // call fix_R() to save atom-pairs with R index (rx, ry, rz) into tmp_atom_pairs
 *           HR.fix_R(rx, ry, rz);
 *           // loop fixed atom-pairs
 *           for (int i = 0; i < HR.size_atom_pairs(); i++)
 *           {
 *              // get pointer of target atom-pair
 *              double* data_pointer = HR.data(i);
 *              // or get reference of target atom-pair
 *              AtomPair<double>& atom_ijR = HR.get_atom_pair(i);
 *              // do something with atom_ijR or data_pointer
 *              ...
 *           }
 *       }
 *       // call unfix_R() to clear tmp_atom_pairs
 *       HR.unfix_R();
 *     ```
 *    b. loop atom-pairs first and then loop R-index
 *     ```
 *       // HR is a const HContainer object, which has been initialized
 *       // loop atom-pairs
 *       for (int i = 0; i < HR.size_atom_pairs(); i++)
 *       {
 *           // get reference of target atom-pair
 *           AtomPair<double>& atom_ij = HR.get_atom_pair(i);
 *           // loop R-index
 *           for (int iR = 0; iR < atom_ij.size_R(); iR++)
 *           {
 *               const int* r_index = atom_ij.get_R_index(iR);
 *               auto tmp_matrix = atom_ij.get_HR_values(r_index[0], r_index[1], r_index[2]);              
 *               // do something with tmp_matrix
 *               ...
 *           }
 *       }
 *     ```
 *    c. loop atom-pairs with gamma_only case
 *     ```
 *       ...
 *       HR.fix_gamma();
 *       // HR is a const HContainer object, which has been initialized and fixed to gamma       
 *       // loop atom-pairs directly without R index
 *       for (int i = 0; i < HR.size_atom_pairs(); i++)
 *       {
 *           // get data pointer of target atom-pair
 *           double* data_pointer = HR.get_pointer(i);
 *           // do something with data_pointer
 *           ...
 *       }
 *     ```
 * 
*/
template <typename T>
class HContainer
{
  public:
    // Destructor of class HContainer
    ~HContainer();

    // copy constructor
    HContainer(const HContainer<T>& HR_in);

    // move constructor
    HContainer(HContainer<T>&& HR_in);

    // simple constructor
    HContainer(int natom);

    // use unitcell to initialize atom_pairs
    HContainer(const UnitCell& ucell_);

    /**
     * @brief use 2d-block-recycle parallel case to initialize atom_pairs, mainly used now.
     * pass a data pointer to HContainer, which means HContainer is a wrapper
     * it will not allocate memory for atom_pairs
     * this case will forbit inserting empty atom_pair
    */
    HContainer(const Parallel_Orbitals* paraV, T* data_pointer = nullptr);


    /**
     * @brief a AtomPair object can be inserted into HContainer, two steps:
     * 1, find AtomPair with atom index atom_i and atom_j
     * 2.1, if found, add to exist AtomPair,
     * 2.2, if not found, insert new one and sort.
     *
     * @param atom_ij AtomPair object
     */
    void insert_pair(const AtomPair<T>& atom_ij);

    /**
     * @brief find AtomPair with atom index atom_i and atom_j
     * This interface can be used to find AtomPair,
     * if found, return pointer will be the exist one,
     * if not found, return pointer will be nullptr.
     *
     * @param i atom index of atom i
     * @param j atom index of atom j
     * @return AtomPair<T>*
     */
    AtomPair<T>* find_pair(int i, int j) const;

    /**
     * @brief return a reference of AtomPair with index of atom I and J in atom_pairs
     *
     * @param i index of atom i
     * @param j index of atom j
     */
    AtomPair<T>& get_atom_pair(int i, int j) const;
    /**
     * @brief return a reference of AtomPair with index in
     * atom_pairs (R is not fixed)
     * tmp_atom_pairs (R is fixed)
     *
     * @param index index of atom-pair
     */
    AtomPair<T>& get_atom_pair(int index) const;

    /**
     * @brief operator() for accessing value with indexes
     * this interface is not implemented now, because it is too expensive to access data
     *
     * @param atom_i index of atom i
     * @param atom_j index of atom j
     * @param rx_in index of R in x direction
     * @param ry_in index of R in y direction
     * @param rz_in index of R in z direction
     * @param mu index of orbital mu
     * @param nu index of orbital nu
     * @return T&
     */
    //T& operator()(int atom_i, int atom_j, int rx_in, int ry_in, int rz_in, int mu, int nu) const;

    // save atom-pair pointers into this->tmp_atom_pairs for selected R index
    /**
     * @brief save atom-pair pointers into this->tmp_atom_pairs for selected R index
     *
     * @param rx_in index of R in x direction
     * @param ry_in index of R in y direction
     * @param rz_in index of R in z direction
     * @return true if success
     */
    bool fix_R(int rx_in, int ry_in, int rz_in) const;

    /**
     * @brief set current_R to -1, which means R is not fixed
     * clear this->tmp_atom_pairs
     */
    void unfix_R() const;

    /**
     * @brief restrict R indexes of all atom-pair to 0, 0, 0
     * add BaseMatrix<T> with non-zero R index to this->atom_pairs[i].values[0]
     * set gamma_only = true
     * in this mode:
     *   1. fix_R() can not be used
     *   2. there is no interface to set gamma_only = false, user should create a new HContainer if needed
     *   3. get_size_for_loop_R() and loop_R() can not be used
     *   4. get_AP_size() can be used
     *   5. data(i, j) can be used to get pointer of target atom-pair with R = 0, 0, 0 , data(i,j,R) can not be used
     *   6. insert_pair() can be used, but the R index will be ignored
     *   7. get_atom_pair(), find_atom_pair() can be used
     *   8. operator() can be used, but the R index will be ignored
     */
    void fix_gamma();

    // interface for call a R loop for HContainer
    // it can return a new R-index with (rx,ry,rz) for each loop
    // if index==0, a new loop of R will be initialized
    /**
     * @brief interface for call a R loop for HContainer
     * it can return a new R-index with (rx,ry,rz) for each loop
     * if index==0, a new loop of R will be initialized
     *
     * @param index index of R loop
     * @param rx index of R in x direction, would be set in the function
     * @param ry index of R in y direction, would be set in the function
     * @param rz index of R in z direction, would be set in the function
     */
    void loop_R(const size_t& index, int& rx, int& ry, int& rz) const;

    /**
     * @brief calculate number of R index which has counted AtomPairs
     */
    size_t size_R_loop() const;

    /**
     * @brief calculate number of AtomPairs for current R index
     */
    size_t size_atom_pairs() const;

    /**
     * @brief get data pointer of AtomPair with index of I, J
     *
     * @param i index of atom i
     * @param j index of atom j
     * @return T* pointer of data
     */
    T* data(int i, int j) const;

    /**
     * @brief get data pointer of AtomPair with index of I, J, Rx, Ry, Rz
     *
     * @param i index of atom i
     * @param j index of atom j
     * @param R int[3] of R index
     * @return T* pointer of data
     */
    T* data(int i, int j, int* R) const;

    /**
     * @brief get current_R
    */
    int get_current_R() const;

    /**
     * @brief judge if gamma_only
    */
    bool is_gamma_only() const;

    /**
     * @brief get total memory bites of HContainer
    */
    size_t get_memory_size() const;

  private:
    // i-j atom pairs, sorted by matrix of (atom_i, atom_j)
    std::vector<AtomPair<T>> atom_pairs;

    // sparse table for (atom_i, atom_j)->index of atom_pairs
    std::vector<std::vector<int>> sparse_ap;
    std::vector<std::vector<int>> sparse_ap_index;

    /**
     * @brief temporary atom-pair lists to loop selected R index
     */
    mutable std::vector<const AtomPair<T>*> tmp_atom_pairs;
    // it contains 3 index of cell, size of R_index is three times of values.
    mutable std::vector<int> tmp_R_index;
    // current index of R in tmp_atom_pairs, -1 means not initialized
    mutable int current_R = -1;
    /**
     * @brief find index of R in tmp_R_index, used when current_R is fixed
     *
     * @param rx_in index of R in x direction
     * @param ry_in index of R in y direction
     * @param rz_in index of R in z direction
     * @return int index of R in tmp_R_index
     */
    int find_R(const int& rx_in, const int& ry_in, const int& rz_in) const;

    bool gamma_only = false;

    /**
     * @brief if wrapper_pointer is not nullptr, this HContainer is a wrapper
    */
    T* wrapper_pointer = nullptr;

    /**
     * @brief pointer of Parallel_Orbitals, which is used to get atom-pair information
    */
    const Parallel_Orbitals* paraV = nullptr;

    // int multiple = 1;
    // int current_multiple = 0;
};

} // namespace hamilt

#endif