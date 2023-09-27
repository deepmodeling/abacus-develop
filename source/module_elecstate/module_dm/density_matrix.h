#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include <string>

#include "module_cell/klist.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

namespace elecstate
{
/**
 * @brief DensityMatrix Class
 * <TK,TR> = <double,double> for Gamma-only calculation
 * <TK,TR> = <std::complex<double>,double> for multi-k calculation
 */
template <typename TK, typename TR>
class DensityMatrix
{
  public:
    /**
     * @brief Destructor of class DensityMatrix
     */
    ~DensityMatrix();

    /**
     * @brief Constructor of class DensityMatrix for multi-k calculation
     * @param _kv pointer of K_Vectors object
     * @param _paraV pointer of Parallel_Orbitals object
     * @param nspin spin setting (1 - none spin; 2 - spin; 4 - SOC)
     */
    DensityMatrix(const K_Vectors* _kv, const Parallel_Orbitals* _paraV, const int nspin);

    /**
     * @brief Constructor of class DensityMatrix for gamma-only calculation, where kvector is not required
     * @param _paraV pointer of Parallel_Orbitals object
     * @param nspin spin setting (1 - none spin; 2 - spin; 4 - SOC)
     */
    DensityMatrix(const Parallel_Orbitals* _paraV, const int nspin);

    /**
     * @brief initialize density matrix DMR from UnitCell
     * @param GridD_in pointer of Grid_Driver object (used to find ajacent atoms)
     * @param ucell pointer of UnitCell object
     */
    void init_DMR(Grid_Driver* GridD_in, const UnitCell* ucell);

    /**
     * @brief initialize density matrix DMR from UnitCell and RA
     * @param ra pointer of Record_adj object (used to find ajacent atoms)
     * @param ucell pointer of UnitCell object
     */
    void init_DMR(Record_adj& ra, const UnitCell* ucell);

    /**
     * @brief initialize density matrix DMR from another HContainer
     * now only support HContainer<double>
     * @param _DMR_in pointer of another HContainer object
     */
    void init_DMR(const hamilt::HContainer<TR>& _DMR_in);

    /// @brief initialize density matrix DMR from another HContainer
    /// this is a temprory function for NSPIN=4 case 
    /// since copy HContainer from another HContainer with different TR is not supported yet
    /// would be refactor in the future
    /// @param _DMR_in 
    void init_DMR(const hamilt::HContainer<std::complex<double>>& _DMR_in);

    /**
     * @brief set _DMK element directly
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @param ik k-point index
     * @param i row index
     * @param j column index
     * @param value value to be set
     */
    void set_DMK(const int ispin, const int ik, const int i, const int j, const TK value);

    /**
     * @brief set _DMK element to zero
    */
    void set_DMK_zero();
    
    /**
     * @brief get a matrix element of density matrix dm(k)
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @param ik k-point index
     * @param i row index
     * @param j column index
     * @return T a matrix element of density matrix dm(k)
     */
    TK get_DMK(const int ispin, const int ik, const int i, const int j) const;

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
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @return HContainer<TR>* pointer of DMR
     */
    hamilt::HContainer<TR>* get_DMR_pointer(const int ispin) const;

    /**
     * @brief get pointer vector of DMR
     * @return HContainer<TR>* vector of DMR
     */
    std::vector<hamilt::HContainer<TR>*> get_DMR_vector() const;

    /**
     * @brief get pointer of DMK
     * @param ik k-point index, which is the index of _DMK
     * @return TK* pointer of DMK
     */
    TK* get_DMK_pointer(const int ik) const;

    /**
     * @brief get pointer vector of DMK
    */
    std::vector<std::vector<TK>> get_DMK_vector() const;

    /**
     * @brief set _DMK using a input TK* pointer
     * please make sure the size of TK* is correct
    */
    void set_DMK_pointer(const int ik, TK* DMK_in);

    /**
     * @brief get pointer of paraV
     */
    const Parallel_Orbitals* get_paraV_pointer() const;
    const K_Vectors* get_kv_pointer() const;

    /**
     * @brief calculate density matrix DMR from dm(k) using blas::axpy
     */
    void cal_DMR();

    /**
     * @brief calculate density matrix DMR from dm(k) using base_matrix->add_element()
     */
    void cal_DMR_test(); // for reference during development

    /**
     * @brief merge density matrix DMR with different spin
     */
    void sum_DMR_spin();

    /**
     * @brief write density matrix dm(ik) into *.dmk
     * @param directory directory of *.dmk files
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @param ik k-point index
     */
    void write_DMK(const std::string directory, const int ispin, const int ik);

    /**
     * @brief read *.dmk into density matrix dm(ik)
     * @param directory directory of *.dmk files
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @param ik k-point index
     */
    void read_DMK(const std::string directory, const int ispin, const int ik);
    
    std::vector<ModuleBase::ComplexMatrix> EDMK; // for TD-DFT

  private:
    /**
     * @brief HContainer for density matrix in real space for 2D parallelization
     * vector.size() = 1 for non-polarization and SOC
     * vector.size() = 2 for spin-polarization
     */
    std::vector<hamilt::HContainer<TR>*> _DMR;

    /**
     * @brief HContainer for density matrix in real space for gird parallelization
     * vector.size() = 1 for non-polarization and SOC
     * vector.size() = 2 for spin-polarization
     */
    std::vector<hamilt::HContainer<TR>*> _DMR_grid;

    /**
     * @brief density matrix in k space, which is a vector[ik]
     * DMK should be a [_nspin][_nks][i][j] matrix,
     * whose size is _nspin * _nks * _paraV->get_nrow() * _paraV->get_ncol()
     */
    // std::vector<ModuleBase::ComplexMatrix> _DMK;
    std::vector<std::vector<TK>> _DMK;

    /**
     * @brief K_Vectors object, which is used to get k-point information
     */
    const K_Vectors* _kv;

    /**
     * @brief Parallel_Orbitals object, which contain all information of 2D block cyclic distribution
     */
    const Parallel_Orbitals* _paraV = nullptr;

    /**
     * @brief spin-polarization index (1 - none spin and SOC ; 2 - spin polarization)
     * Attention: this is not as same as GlovalV::NSPIN
     * _nspin means the number of isolated spin-polarization states
     */
    int _nspin = 1;

    /**
     * @brief real number of k-points
     * _nks is not equal to _kv->get_nks() when spin-polarization is considered
     * _nks = kv->_nks / nspin
     */
    int _nks = 0;
};

} // namespace elecstate

#endif
