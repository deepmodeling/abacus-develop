#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include <string>

#include "module_cell/klist.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

namespace elecstate
{
/**
 * @brief DensityMatrix Class
 * Now only support T = double
 */
template <typename T>
class DensityMatrix
{
  public:
    /**
     * @brief Destructor of class DensityMatrix
     */
    ~DensityMatrix();

    /**
     * @brief Constructor of class DensityMatrix
     * @param _kv pointer of K_Vectors object
     * @param _paraV pointer of Parallel_Orbitals object
     */
    DensityMatrix(const K_Vectors* _kv, const Parallel_Orbitals* _paraV);

    /**
     * @brief initialize density matrix DMR from UnitCell
     * @param GridD_in pointer of Grid_Driver object (used to find ajacent atoms)
     * @param ucell pointer of UnitCell object
     */
    void init_DMR(Grid_Driver* GridD_in, const UnitCell* ucell);

    /**
     * @brief initialize density matrix DMR from another HContainer
     * @param _DMR_in pointer of another HContainer object
     */
    void init_DMR(const hamilt::HContainer<T>& _DMR_in);

    /**
     * @brief set _DMK element directly
     * @param ik k-point index
     * @param i row index
     * @param j column index
     * @param value value to be set
     */
    void set_DMK(const int ik, const int i, const int j, const T value);

    /**
     * @brief read *.dmk into density matrix dm(k) with all k-points, only support serial version now
     * @param directory directory of *.dmk files
     */
    void set_DMK_files(const std::string directory);

    /**
     * @brief write density matrix dm(ik) into *.dmk
     * @param directory directory of *.dmk files
     * @param ik k-point index
     */
    void write_DMK(const std::string directory, const int ik);

    /**
     * @brief read *.dmk into density matrix dm(ik)
     * @param directory directory of *.dmk files
     * @param ik k-point index
     */
    void read_DMK(const std::string directory, const int ik);

    /**
     * @brief write density matrix dm(k) into *.dmk with all k-points
     * @param directory directory of *.dmk files
     */
    void output_DMK(const std::string directory);

    /**
     * @brief get a matrix element of density matrix dm(k)
     * @param ik k-point index
     * @param i row index
     * @param j column index
     * @return T a matrix element of density matrix dm(k)
     */
    T get_DMK(const int ik, const int i, const int j) const;

    /**
     * @brief get total number of k-points of density matrix dm(k)
     */
    int get_DMK_nks() const;

    /**
     * @brief get number of rows of density matrix dm(k)
     */
    int get_DMK_nrow() const;

    /**
     * @brief get number of columns of density matrix dm(k)
     */
    int get_DMK_ncol() const;

    /**
     * @brief get pointer of DMR
     */
    hamilt::HContainer<T>* get_DMR_pointer();

    /**
     * @brief calculate density matrix DMR from dm(k)
     */
    void cal_DMR();

    /**
     * @brief calculate density matrix DMR from dm(k) using blas::axpy
     */
    void cal_DMR_blas();

  private:
    /**
     * @brief HContainer for density matrix in real space
     */
    hamilt::HContainer<T>* _DMR = nullptr;

    /**
     * @brief density matrix in k space, which is a vector[ik]
     */
    std::vector<ModuleBase::ComplexMatrix> _DMK;

    /**
     * @brief K_Vectors object, which is used to get k-point information
     */
    const K_Vectors* _kv;

    /**
     * @brief Parallel_Orbitals object, which contain all information of 2D block cyclic distribution
     */
    const Parallel_Orbitals* _paraV = nullptr;
};

} // namespace elecstate

#endif
