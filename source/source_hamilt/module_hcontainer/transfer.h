#ifndef HCONTAINERTRANSFER_H
#define HCONTAINERTRANSFER_H

#include <unordered_map>
#include <vector>

#include "./hcontainer.h"

#ifdef __MPI
namespace hamilt
{

// MPI tags for the HContainer serial<->parallel transfer (transfer.cpp) and
// Read_HContainer::read() (read_hcontainer.cpp). They must match between the
// send and the recv side, so they live here to be shared. Deliberately odd,
// high values: the transfer previously used fixed tags 0/1/2 on
// MPI_COMM_WORLD, which collide with unrelated p2p traffic (write_eig_occ,
// output_log, func_transfer, ...) on any other call path and cause
// MPI_ERR_TRUNCATE races or deadlocks. Moving to these uncommon tags removes
// the need for the masking MPI_Barrier calls in op_exx_lcao.cpp.
constexpr int kHctTagSize = 14431;   // payload size (MPI_LONG)
constexpr int kHctTagIndex = 14471;  // atom/orbital index array (MPI_INT)
constexpr int kHctTagValue = 14533;  // matrix values (T)

template <typename T>
class HTransPara
{
  public:
    HTransPara(int n_processes, HContainer<T>* hr_in);
    ~HTransPara();

    /**
     * @brief calculate Orbital indexes and will be send to irank
     * the format of the return value is:
     * [size_i, atom_i, number_orb_row, 0, 1, 2, 3, 8, 9, 10, 11, number_orb_col, ... atom_i, number_orb_row, 4, 5, 6,
     * 7, atom_i, ... number_orb_col, ...] i refers to the ith atom in this->ap_indexes[irank] the function is called in
     * plan_indexes
     * @param irank
     * @param orb_indexes
     */
    void cal_orb_indexes(int irank, std::vector<int>* orb_indexes = nullptr);

    /**
     * @brief receive AtomPair_indexes from the ith rank
     * save to this->ap_indexes[irank]
     * @param irank
     */
    void receive_ap_indexes(int irank, const int* ap_indexes_in = nullptr, const long& size_ap_indexes_in = 0);

    /**
     * @brief pack data in this->hr, and send to ith rank
     * @param irank
     */
    void send_orb_indexes(int irank, MPI_Request* request = nullptr);

    /**
     * @brief pack data in this->hr, and send to ith rank
     * @param irank
     */
    void send_data(int irank, MPI_Request* request = nullptr);

    /**
     * @brief receive data from ith rank, save them to this->hr
     * @param irank
     */
    void receive_data(int irank, const T* values = nullptr);

    /**
     * @brief pack BaseMatrix-values for ith rank
     * @param irank
     */
    void pack_data(int irank, T* values = nullptr);

    long get_max_size() const;
    void get_value_size(int* out) const;

  private:
    std::vector<std::vector<int>> ap_indexes;
    HContainer<T>* hr = nullptr;

    // temporary variables
    std::vector<std::vector<int>> atom_i_index;

    const Parallel_Orbitals* paraV = nullptr;

    // unpack BaseMatrix-values from ith rank
    void unpack_data(int irank, const T* values);

    // size of data of all BaseMatrixes
    std::vector<long> size_values;
};

template <typename T>
class HTransSerial
{
  public:
    HTransSerial(int n_processes, HContainer<T>* hr_in);
    ~HTransSerial();

    /**
     * @brief calculate AtomPair indexes and will be send to irank
     * called in plan_indexes
     * @param irank
     */
    void cal_ap_indexes(int irank, std::vector<int>* ap_indexes = nullptr);

    /**
     * @brief calculate AtomPair_indexes of hr_in and send to the ith rank
     * @param irank
     */
    void send_ap_indexes(int irank, MPI_Request* request = nullptr);

    /**
     * @brief receive Orbital_indexes from the ith rank
     * save to this->orb_indexes[irank]
     */
    void receive_orb_indexes(int irank, const int* orb_indexes_in = nullptr, const long& size_orb_indexes_in = 0);

    /**
     * @brief pack data in this->hr, and send to ith rank
     * @param irank
     */
    void send_data(int irank, MPI_Request* request = nullptr);

    /**
     * @brief receive data from ith rank, save them to this->hr
     * @param irank
     */
    void receive_data(int irank, const T* values = nullptr);

    /**
     * @brief pack BaseMatrix-values for ith rank
     * @param irank
     */
    void pack_data(int irank, T* values = nullptr);

    long get_max_size() const;
    void get_value_size(int* out) const;

  private:
    std::vector<std::vector<int>> orb_indexes;
    HContainer<T>* hr = nullptr;

    /**
     * @brief unpack BaseMatrix-values from ith rank
     * @param values
     */
    void unpack_data(int irank, const T* values);

    // temporary variables
    std::vector<std::unordered_map<int, int>> orb_col_indexes;
    std::vector<std::unordered_map<int, int>> orb_row_indexes;

    // size of data of all BaseMatrixes
    std::vector<long> size_values;
};

} // namespace hamilt

/**
 * @brief Struct to get MPI_datatype
 */
template <typename T>
struct MPITraits;

template <>
struct MPITraits<int>
{
    static MPI_Datatype datatype()
    {
        return MPI_INT;
    }
};

template <>
struct MPITraits<double>
{
    static MPI_Datatype datatype()
    {
        return MPI_DOUBLE;
    }
};

template <>
struct MPITraits<std::complex<double>>
{
    static MPI_Datatype datatype()
    {
        return MPI_DOUBLE_COMPLEX;
    }
};

#endif // __MPI

#endif // HCONTAINERTRANSFER_H