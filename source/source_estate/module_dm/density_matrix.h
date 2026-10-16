#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include <complex>
#include <map>
#include <string>
#include <vector>

#include "source_base/vector3.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/record_adj.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

namespace module_dm
{
/**
 * @brief map a real/complex type to the opposite one
 * ShiftRealComplex<double>::type = std::complex<double>
 * ShiftRealComplex<std::complex<double>>::type = double
 */
template <typename T> struct ShiftRealComplex
{
    using type = void;
};

template <>
struct ShiftRealComplex<double>
{
    using type = std::complex<double>;
};

template <>
struct ShiftRealComplex<std::complex<double>>
{
    using type = double;
};

/**
 * @brief DensityMatrix Class
 * <TK,TR> = <double,double> for Gamma-only calculation
 * <TK,TR> = <std::complex<double>,double> for multi-k calculation
 */
    template <typename TK, typename TR> class DensityMatrix;

// DensityMatrix<complex<double>,TR>::cal_dmr() is illegal in C++, so module_dm is used instead.
    template <typename TK, typename TR_in, typename TR_out>
    extern void cal_dmr(
        DensityMatrix<TK, TR_in> &dm,
        std::vector<hamilt::HContainer<TR_out>*> &dmR_out,
        const int ik_in);

    template <typename TK, typename TR_in, typename TR_out>
    extern void cal_dmr_td(
        DensityMatrix<TK, TR_in> &dm,
        std::vector<hamilt::HContainer<TR_out>*> &dmR_out,
        const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
        const ModuleBase::Vector3<double> At,
        const int ik_in);

    template <typename TK, typename TR_in, typename TR_out>
    extern void cal_dmr_full(
        const DensityMatrix<TK, TR_in> &dm,
        hamilt::HContainer<TR_out>* dmR_out,
        const int ik_in);

    /**
     * @brief shared inner loop of cal_dmr / cal_dmr_td: for each spin channel,
     * zero the DMR HContainer and accumulate kphase * DMK into DMR blocks.
     * Pass an empty phase_hybrid map for the non-TD (cal_dmr) case.
     */
    template <typename TK, typename TR_in, typename TR_out>
    extern void accumulate_dmr(
        DensityMatrix<TK, TR_in> &dm,
        std::vector<hamilt::HContainer<TR_out>*> &dmR_out,
        const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
        const int ik_in,
        const char* func_name);

    template <typename TR>
    extern void exp_mul_dmk(const std::complex<double> kphase,
                                const std::vector<std::complex<double>>& dmk_row,
                                TR* dmr_mat);

    template <typename TR>
    extern void xyz_to_updown(const std::complex<double> spin_block[4],
                                  const int icol,
                                  const int spin_stride[4],
                                  TR* dmr_mat);

    /**
     * @brief geometry of one atom-pair sub-block within the global DMK matrix
     * row0/col0: global index of the block's top-left element in the 2D block-cyclic DMK
     * nrows/ncols: orbital dimensions of the two atoms
     */
    struct DmrBlock
    {
        int row0;
        int col0;
        int nrows;
        int ncols;
        int size() const { return nrows * ncols; }
    };

    /// @brief extract the block geometry for atom pair (iat1, iat2) from the parallel orbitals layout
    DmrBlock get_dmr_block(const Parallel_Orbitals* pv, const int iat1, const int iat2);

    /**
     * @brief precompute k-phase factors e^{ikR} and collect DMR block pointers for one atom pair
     * @param atom_pair the atom pair whose R-vectors and matrices are used
     * @param kvec_d direct coordinates of k-points
     * @param nk number of k-points
     * @param phase_hybrid additional hybrid-gauge phase per R (empty map = no extra phase)
     * @param kphase_vec output: kphase_vec[ik][iR]
     * @param dmr_mats output: dmr_mats[iR] points to the DMR block for R-vector iR
     */
    template <typename TK, typename TR>
    extern void build_kphase(hamilt::AtomPair<TR>& atom_pair,
                             const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                             const int nk,
                             const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
                             std::vector<std::vector<TK>>& kphase_vec,
                             std::vector<TR*>& dmr_mats);

    /// @brief transpose a col-major DMK sub-block into row-major order
    template <typename TK>
    extern void transpose_dmk_block(const TK* dmk_col_major,
                                    const int ld_hk,
                                    const DmrBlock& block,
                                    TK* dmk_row);

    /**
     * @brief nspin=1/2: accumulate Re(kphase * DMK) into DMR blocks
     *
     * Formula: DMR_ij(R) += Re[ e^{ik·R} * DMK_ij(k) ]
     * If ik_in >= 0, only that k-point contributes; if ik_in < 0, sum over all k-points.
     */
    template <typename TK, typename TR>
    extern void add_dmr_real(const DensityMatrix<TK, TR>& dm,
                             const DmrBlock& block,
                             const int ik_begin,
                             const std::vector<std::vector<TK>>& kphase_vec,
                             const int ld_hk,
                             const int ik_in,
                             std::vector<TR*>& dmr_mats);

    /**
     * @brief nspin==4 (SOC): accumulate k-phase * DMK into a per-R complex buffer,
     * then transform 2x2 spin blocks from (upup, updown, downup, downdown) to
     * (rho_0, rho_x, rho_y, rho_z) via xyz_to_updown.
     *
     * Formula:
     *   S_ij(R) = sum_k e^{ik·R} * DMK_ij(k)
     *   rho_0 = rho_upup + rho_downdown
     *   rho_x = rho_updown + rho_downup
     *   rho_y = Im(rho_updown) - Im(rho_downup)   (sign for conjugated stored DM)
     *   rho_z = rho_upup - rho_downdown
     * Each orbital corresponds to a 2x2 spin block, so rows/cols step by 2.
     * If ik_in >= 0, only that k-point contributes; if ik_in < 0, sum over all k-points.
     */
    template <typename TK, typename TR>
    extern void add_dmr_soc(const DensityMatrix<TK, TR>& dm,
                            const DmrBlock& block,
                            const int ik_begin,
                            const std::vector<std::vector<TK>>& kphase_vec,
                            const int ld_hk,
                            const int ik_in,
                            const int col_stride,
                            std::vector<TR*>& dmr_mats);

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
     * @param pv pointer of Parallel_Orbitals object
     * @param spin_mult spin multiplicity used to size the DM: 1 for input nspin 1 or 4
     *  (non-collinear k points are not doubled), 2 for input nspin 2 (LSDA up/down).
     *  This is NOT the physical nspin (1/2/4); it matches K_Vectors::spin_mult.
     * @param kvec_d direct coordinates of kpoints
     * @param nk number of k-points, not always equal to K_Vectors::get_nks()/spin_mult.
     *               it will be set to kvec_d.size() if the value is invalid
     * @param nspin the global physical nspin from INPUT (1/2/4); defaults to spin_mult for
     *               non-SOC cases where they coincide. Pass 4 explicitly for SOC/noncollinear
     *               calculations so that cal_dmr selects the spin-resolved (Pauli) branch.
     */
    DensityMatrix(const Parallel_Orbitals* pv,
            const int spin_mult,
            const std::vector<ModuleBase::Vector3<double>>& kvec_d,
            const int nk,
            const int nspin = 0);

    /**
     * @brief Constructor of class DensityMatrix for gamma-only calculation, where kvector is not required
     * @param pv pointer of Parallel_Orbitals object
     * @param spin_mult spin multiplicity of the density matrix (1 or 2); NOT the physical nspin
     * @param nspin the global physical nspin from INPUT (1/2/4); defaults to spin_mult.
     */
    DensityMatrix(const Parallel_Orbitals* pv, const int spin_mult, const int nspin = 0);

    /**
     * @brief initialize density matrix DMR from UnitCell
     * @param GridD_in pointer of Grid_Driver object (used to find ajacent atoms)
     * @param ucell pointer of UnitCell object
     */
    void init_dmr(const Grid_Driver* GridD_in, const UnitCell* ucell);

    /**
     * @brief initialize density matrix DMR from UnitCell and RA
     * @param ra pointer of Record_adj object (used to find ajacent atoms)
     * @param ucell pointer of UnitCell object
     */
    void init_dmr(Record_adj& ra, const UnitCell* ucell);

    /**
     * @brief initialize density matrix DMR from another HContainer
     * now only support HContainer<double>
     * @param _DMR_in pointer of another HContainer object
     */
    void init_dmr(const hamilt::HContainer<TR>& _DMR_in);

    /// @brief initialize density matrix DMR from another HContainer
    /// this is a temprory function for NSPIN=4 case 
    /// since copy HContainer from another HContainer with different TR is not supported yet
    /// would be refactor in the future
    /// @param _DMR_in 
    // the old input type ``:HContainer<complex<double>` causes redefination error if TR = complex<double>
    void init_dmr(const hamilt::HContainer<TRShift>& _DMR_in);

    /**
     * @brief set dmk element directly
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @param ik k-point index
     * @param i row index
     * @param j column index
     * @param value value to be set
     */
    void set_dmk(const int ispin, const int ik, const int i, const int j, const TK value);

    /**
     * @brief set dmk element to zero
    */
    void set_dmk_zero();
    
    /**
     * @brief get a matrix element of density matrix dm(k)
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @param ik k-point index
     * @param i row index
     * @param j column index
     * @return T a matrix element of density matrix dm(k)
     */
    TK get_dmk(const int ispin, const int ik, const int i, const int j) const;

    /**
     * @brief get total number of k-points of density matrix dm(k)
     */
    int get_dmk_nks() const;
    int get_dmk_size() const;

    /**
     * @brief get number of rows of density matrix dm(k)
     */
    int get_dmk_nrow() const;

    /**
     * @brief get number of columns of density matrix dm(k)
     */
    int get_dmk_ncol() const;

    /**
     * @brief get pointer of DMR
     * @param ispin spin index (1 - spin up (support SOC) or 2 - spin down)
     * @return HContainer<TR>* pointer of DMR
     */
    hamilt::HContainer<TR>* get_dmr_ptr(const int ispin) const;

    /**
     * @brief check whether the stored DMR is a valid density matrix calculated from DMK
     * init_dmr() resets the flag and cal_dmr()/cal_dmr_td() set it, so a freshly
     * allocated, zeroed or file-read DMR is reported as not ready until the first
     * wavefunction-derived calculation
     * @return true if DMR is ready for Hamiltonian construction
     */
    bool is_dmr_ready() const
    {
        return this->_dmr_ready;
    }

    /**
     * @brief get pointer vector of DMR
     * @return HContainer<TR>* vector of DMR
     */
    const std::vector<hamilt::HContainer<TR>*>& get_dmr_vec() const
    {
        return this->dmr;
    }
    std::vector<hamilt::HContainer<TR>*>& get_dmr_vec()
    {
        return this->dmr;
    }

    const std::vector<std::vector<TR>>& get_dmr_save() const
    {
        return this->dmr_save;
    }
    std::vector<std::vector<TR>>& get_dmr_save()
    {
        return this->dmr_save;
    }

    /**
     * @brief get pointer of DMK
     * @param ik k-point index, which is the index of dmk
     * @return TK* pointer of DMK
     */
    TK* get_dmk_ptr(const int ik) const;

    /**
     * @brief get pointer vector of DMK
    */
    const std::vector<std::vector<TK>>& get_dmk_vec() const
    {
        return this->dmk;
    }
    std::vector<std::vector<TK>>& get_dmk_vec()
    {
        return this->dmk;
    }

    /**
     * @brief set dmk using a input TK* pointer
     * please make sure the size of TK* is correct
    */
    void set_dmk_ptr(const int ik, TK* DMK_in);

    /**
     * @brief calculate density matrix DMR from dm(k) using blas::axpy
     * @param ik_in
     * if ik_in < 0, calculate all k-points
     * if ik_in >= 0, calculate only one k-point without summing over k-points
     */
    void cal_dmr(const int ik_in);

    /**
     * @brief calculate density matrix DMR with additional vector potential phase, used for hybrid gauge tddft
     * @param ik_in
     * if ik_in < 0, calculate all k-points
     * if ik_in >= 0, calculate only one k-point
     */
    void cal_dmr_td(const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
                    const ModuleBase::Vector3<double> At,
                    const int ik_in);

    /**
     * @brief calculate complex density matrix DMR with both real and imaginary part for noncollinear-spin calculation
     * the stored dm(k) has been used to calculate the passin DMR
     * @param dmR_out pointer of HContainer object to store the calculated complex DMR
     * @param ik_in
     * if ik_in < 0, calculate all k-points
     * if ik_in >= 0, calculate only one k-point
     */
    void cal_dmr_full(hamilt::HContainer<std::complex<double>>* dmR_out, const int ik_in) const;

    /**
     * @brief (Only nspin=2) switch DMR to total density matrix or magnetization density matrix
     * @param mode 0 - original density matrix; 1 - total density matrix; 2 - magnetization density matrix
     */
    void switch_dmr(const int mode);

    /**
     * @brief save dmr into dmr_save
     */
    void save_dmr();
    
    std::vector<ModuleBase::ComplexMatrix> edmk; // for TD-DFT

#ifdef __PEXSI
    /**
     * @brief EDM storage for PEXSI
     * used in MD calculation
     */
    std::vector<TK*> edm_pexsi;
#endif

  private:
    /**
     * @brief delete all HContainer objects in dmr and clear the vector
     */
    void clear_dmr();

    /**
     * @brief HContainer for density matrix in real space for 2D parallelization
     * vector.size() = 1 for non-polarization and SOC
     * vector.size() = 2 for spin-polarization
     */
    std::vector<hamilt::HContainer<TR>*> dmr;
    std::vector<std::vector<TR>> dmr_save;

    /// @brief whether dmr holds a density matrix calculated from DMK (reset by init_dmr, set by cal_dmr)
    bool _dmr_ready = false;

    /**
     * @brief HContainer for density matrix in real space for grid parallelization
     * same size semantics as dmr
     */
    std::vector<hamilt::HContainer<TR>*> dmr_grid;

    /**
     * @brief density matrix in k space, which is a vector[ik]
     * DMK should be a [spin_mult][_nk][i][j] matrix,
     * whose size is spin_mult * _nk * pv->get_nrow() * pv->get_ncol()
     */
    // std::vector<ModuleBase::ComplexMatrix> dmk;
    std::vector<std::vector<TK>> dmk;

    /**
     * @brief K_Vectors object, which is used to get k-point information
     */
    const std::vector<ModuleBase::Vector3<double>> _kvec_d;

    /**
     * @brief Parallel_Orbitals object, which contain all information of 2D block cyclic distribution
     */
    const Parallel_Orbitals* pv = nullptr;

    /**
     * @brief spin multiplicity used to size the density matrix (1 - none spin and SOC ;
     * 2 - spin polarization). This is NOT the physical nspin (1/2/4); it matches
     * K_Vectors::spin_mult.
     */
    int spin_mult = 1;

    /**
     * @brief the global physical nspin from INPUT (1/2/4).
     * For SOC/noncollinear (nspin==4) the density matrix is stored with spin_mult==1
     * (a single 2x2 spin-block matrix), but cal_dmr/cal_dmr_td must still take the
     * spin-resolved (Pauli) branch. spin_mult cannot distinguish this, so keep the
     * global value here. Equals spin_mult for non-SOC cases.
     */
    int nspin = 1;

    /**
     * @brief real number of k-points
     * _nk is not equal to _kv->get_nks() when spin-polarization is considered
     * _nk = kv->get_nks() / nspin when nspin=2
     */
    int _nk = 0;

    /// temporary pointers for switch DMR, only used with nspin=2
    std::vector<TR> dmr_origin;
    std::vector<TR> dmr_tmp;

    friend void module_dm::cal_dmr<TK, TR>(
        DensityMatrix<TK, TR>& dm,
        std::vector<hamilt::HContainer<TR>*>& dmR_out,
        const int ik_in);
    friend void module_dm::cal_dmr_td<TK, TR>(
        DensityMatrix<TK, TR>& dm,
        std::vector<hamilt::HContainer<TR>*>& dmR_out,
        const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
        const ModuleBase::Vector3<double> At,
        const int ik_in);
    friend void module_dm::cal_dmr_full<TK, TR>(
        const DensityMatrix<TK, TR>& dm,
        hamilt::HContainer<std::complex<double>>* dmR_out,
        const int ik_in);
    friend void module_dm::accumulate_dmr<TK, TR>(
        DensityMatrix<TK, TR>& dm,
        std::vector<hamilt::HContainer<TR>*>& dmR_out,
        const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
        const int ik_in,
        const char* func_name);
    friend void module_dm::add_dmr_real<TK, TR>(
        const DensityMatrix<TK, TR>& dm,
        const DmrBlock& block,
        const int ik_begin,
        const std::vector<std::vector<TK>>& kphase_vec,
        const int ld_hk,
        const int ik_in,
        std::vector<TR*>& dmr_mats);

    friend void module_dm::add_dmr_soc<TK, TR>(
        const DensityMatrix<TK, TR>& dm,
        const DmrBlock& block,
        const int ik_begin,
        const std::vector<std::vector<TK>>& kphase_vec,
        const int ld_hk,
        const int ik_in,
        const int col_stride,
        std::vector<TR*>& dmr_mats);
};

} // namespace module_dm

#endif
