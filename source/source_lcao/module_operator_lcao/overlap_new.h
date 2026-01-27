#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_OPERATOR_LCAO_OVERLAP_NEW_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_OPERATOR_LCAO_OVERLAP_NEW_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include <vector>

namespace hamilt
{

#ifndef __OVERLAPNEWTEMPLATE
#define __OVERLAPNEWTEMPLATE

/// The OverlapNew class template inherits from class T
/// it is used to calculate the overlap of wavefunction basis
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK> or OperatorPW<TK>
/// - TR: data type of real space Hamiltonian, it would be double or std::complex<double>
template <class T>
class OverlapNew : public T
{
};

#endif

/// OverlapNew class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the overlap matrix in real space and fold it to k-space
/// SR = <psi_{mu, 0}|psi_{nu, R}>
/// SK = <psi_{mu, k}|psi_{nu, k}> = \sum_{R} e^{ikR} SR
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class OverlapNew<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    OverlapNew<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                     const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                     hamilt::HContainer<TR>* hR_in,
                                     hamilt::HContainer<TR>* SR_in,
                                     const UnitCell* ucell_in,
                                     const std::vector<double>& orb_cutoff,
                                     const Grid_Driver* GridD_in,
                                     const TwoCenterIntegrator* intor);

    ~OverlapNew();

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

    TK* getSk();

    /**
     * @brief calculate force and stress for overlap operator
     * @param cal_force whether to calculate force
     * @param cal_stress whether to calculate stress
     * @param dmR density matrix in real space
     * @param force output force matrix (nat x 3)
     * @param stress output stress matrix (3 x 3)
     */
    void cal_force_stress(const bool cal_force,
                          const bool cal_stress,
                          const HContainer<double>* dmR,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& stress);

  private:
    const UnitCell* ucell = nullptr;

    std::vector<double> orb_cutoff_;

    hamilt::HContainer<TR>* SR = nullptr;

    const TwoCenterIntegrator* intor_ = nullptr;

    const Grid_Driver* gridD = nullptr;

    bool SR_fixed_done = false;

    /**
     * @brief initialize SR, search the nearest neighbor atoms
     * HContainer is used to store the overlap matrix with specific <I,J,R> atom-pairs
     * the size of SR will be fixed after initialization
     */
    void initialize_SR(const Grid_Driver* GridD_in);

    /**
     * @brief calculate the overlap matrix with specific <I,J,R> atom-pairs
     * nearest neighbor atoms don't need to be calculated again
     * loop the atom-pairs in SR and calculate the overlap matrix
     */
    void calculate_SR();

    /**
     * @brief calculate the SR local matrix of <I,J,R> atom pair
     */
    void cal_SR_IJR(const int& iat1,
                    const int& iat2,
                    const Parallel_Orbitals* paraV,
                    const ModuleBase::Vector3<double>& dtau,
                    TR* data_pointer);

    /**
     * @brief calculate force contribution for atom pair <I,J,R>
     */
    void cal_force_IJR(const int& iat1,
                       const int& iat2,
                       const Parallel_Orbitals* paraV,
                       const std::unordered_map<int, std::vector<double>>& nlm1_all,
                       const std::unordered_map<int, std::vector<double>>& nlm2_all,
                       const hamilt::BaseMatrix<TR>* dmR_pointer,
                       double* force1,
                       double* force2);

    /**
     * @brief calculate stress contribution for atom pair <I,J,R>
     */
    void cal_stress_IJR(const int& iat1,
                        const int& iat2,
                        const Parallel_Orbitals* paraV,
                        const std::unordered_map<int, std::vector<double>>& nlm1_all,
                        const std::unordered_map<int, std::vector<double>>& nlm2_all,
                        const hamilt::BaseMatrix<TR>* dmR_pointer,
                        const ModuleBase::Vector3<double>& dis1,
                        const ModuleBase::Vector3<double>& dis2,
                        double* stress);

    // if k vector is not changed, then do nothing and return
    // default of kvec_d_old is (-10,-10,-10), which is not a valid k vector
    ModuleBase::Vector3<double> kvec_d_old = ModuleBase::Vector3<double>(-10, -10, -10);

  public:
    /**
     * @brief calculate asynchronous overlap matrix for Hefei-NAMD
     * Calculates <phi(t-1)|phi(t)> by shifting atom positions backward
     * @param ucell unit cell with current atomic positions and velocities
     * @param md_dt molecular dynamics time step (in fs)
     * @param paraV parallel orbitals object for matrix distribution
     * @return pointer to the created SR_async container (caller must delete)
     */
    hamilt::HContainer<TR>* calculate_SR_async(const UnitCell& ucell, const double md_dt, const Parallel_Orbitals* paraV);

    /**
     * @brief output asynchronous overlap matrix in CSR format
     * @param istep current ionic step number
     * @param SR_async pointer to the asynchronous overlap matrix container
     * @param precision output precision for floating point numbers
     */
    void output_SR_async_csr(const int istep, hamilt::HContainer<TR>* SR_async, const int precision = 8);
};

} // namespace hamilt
#endif
