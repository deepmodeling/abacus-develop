#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include <string>

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
template<typename T> struct ShiftRealComplex
{
    using type = void;
};

template<>
struct ShiftRealComplex<double> 
{
	using type = std::complex<double>;
};

template<>
struct ShiftRealComplex<std::complex<double>> 
{
	using type = double;
};

template <typename TK, typename TR>
class DensityMatrix
{
	using TRShift = typename ShiftRealComplex<TR>::type;

	public:
	/**
	 * @brief Destructor of class DensityMatrix
	 */
	~DensityMatrix();

    /**
     * @brief Constructor of class DensityMatrix for multi-k calculation
     * @param _paraV pointer of Parallel_Orbitals object
     * @param nspin number of spin of the density matrix, set by user according to global nspin
     *  (usually {nspin_global -> nspin_dm} = {1->1, 2->2, 4->1}, but sometimes 2->1 like in LR-TDDFT)
     * @param kvec_d direct coordinates of kpoints
     * @param nk number of k-points, not always equal to K_Vectors::get_nks()/nspin_dm.
     *               it will be set to kvec_d.size() if the value is invalid
     */
    DensityMatrix(const Parallel_Orbitals* _paraV, const int nspin, const std::vector<ModuleBase::Vector3<double>>& kvec_d, const int nk);

    /**
     * @brief Constructor of class DensityMatrix for gamma-only calculation, where kvector is not required
     * @param _paraV pointer of Parallel_Orbitals object
     * @param nspin number of spin of the density matrix, set by user according to global nspin
     *  (usually {nspin_global -> nspin_dm} = {1->1, 2->2, 4->1}, but sometimes 2->1 like in LR-TDDFT)
     */
    DensityMatrix(const Parallel_Orbitals* _paraV, const int nspin);

    /**
     * @brief initialize density matrix DMR from UnitCell
     * @param GridD_in pointer of Grid_Driver object (used to find ajacent atoms)
     * @param ucell pointer of UnitCell object
     */
    void init_DMR(const Grid_Driver* GridD_in, const UnitCell* ucell);

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
    // the old input type ``:HContainer<complex<double>` causes redefination error if TR = complex<double>
    void init_DMR(const hamilt::HContainer<TRShift>& _DMR_in);

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
    int get_DMK_size() const;

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
    const std::vector<hamilt::HContainer<TR>*>& get_DMR_vector() const {return this->_DMR;}
    std::vector<hamilt::HContainer<TR>*>& get_DMR_vector() {return this->_DMR;}

    const std::vector<std::vector<TR>>& get_DMR_save() const {return this->_DMR_save;}
    std::vector<std::vector<TR>>& get_DMR_save() {return this->_DMR_save;}

    /**
     * @brief get pointer of DMK
     * @param ik k-point index, which is the index of _DMK
     * @return TK* pointer of DMK
     */
    TK* get_DMK_pointer(const int ik) const;

    /**
     * @brief get pointer vector of DMK
    */
    const std::vector<std::vector<TK>>& get_DMK_vector() const {return this->_DMK;}
    std::vector<std::vector<TK>>& get_DMK_vector() {return this->_DMK;}

    /**
     * @brief set _DMK using a input TK* pointer
     * please make sure the size of TK* is correct
    */
    void set_DMK_pointer(const int ik, TK* DMK_in);

    /**
     * @brief get pointer of paraV
     */
    const Parallel_Orbitals* get_paraV_pointer() const {return this->_paraV;}

    const std::vector<ModuleBase::Vector3<double>>& get_kvec_d() const { return this->_kvec_d; }

    /**
     * @brief calculate density matrix DMR from dm(k) using blas::axpy
     * if ik_in < 0, calculate all k-points
     * if ik_in >= 0, calculate only one k-point without summing over k-points
     */
    void cal_DMR(const int ik_in = -1);

    /**
     * @brief calculate complex density matrix DMR with both real and imaginary part for noncollinear-spin calculation
     * the stored dm(k) has been used to calculate the passin DMR
     * @param dmR_out pointer of HContainer object to store the calculated complex DMR
     */
    void cal_DMR_full(hamilt::HContainer<std::complex<double>>* dmR_out) const;

    /**
     * @brief (Only nspin=2) switch DMR to total density matrix or magnetization density matrix
     * @param mode 0 - original density matrix; 1 - total density matrix; 2 - magnetization density matrix
     */
    void switch_dmr(const int mode);

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

    /**
     * @brief save _DMR into _DMR_save
     */
    void save_DMR();
    
    std::vector<ModuleBase::ComplexMatrix> EDMK; // for TD-DFT

#ifdef __PEXSI
    /**
     * @brief EDM storage for PEXSI
     * used in MD calculation
     */
    std::vector<TK*> pexsi_EDM;
#endif

  private:
    /**
     * @brief HContainer for density matrix in real space for 2D parallelization
     * vector.size() = 1 for non-polarization and SOC
     * vector.size() = 2 for spin-polarization
     */
    std::vector<hamilt::HContainer<TR>*> _DMR;
    std::vector<std::vector<TR>> _DMR_save;

    /**
     * @brief HContainer for density matrix in real space for gird parallelization
     * vector.size() = 1 for non-polarization and SOC
     * vector.size() = 2 for spin-polarization
     */
    std::vector<hamilt::HContainer<TR>*> _DMR_grid;

    /**
     * @brief density matrix in k space, which is a vector[ik]
     * DMK should be a [_nspin][_nk][i][j] matrix,
     * whose size is _nspin * _nk * _paraV->get_nrow() * _paraV->get_ncol()
     */
    // std::vector<ModuleBase::ComplexMatrix> _DMK;
    std::vector<std::vector<TK>> _DMK;

    /**
     * @brief K_Vectors object, which is used to get k-point information
     */
    const std::vector<ModuleBase::Vector3<double>> _kvec_d;

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
     * _nk is not equal to _kv->get_nks() when spin-polarization is considered
     * _nk = kv->get_nks() / nspin when nspin=2
     */
    int _nk = 0;

    /// temporary pointers for switch DMR, only used with nspin=2
    std::vector<TR> dmr_origin_;
    TR* dmr_tmp_ = nullptr;

};

} // namespace elecstate

#endif
