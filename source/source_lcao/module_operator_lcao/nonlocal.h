#ifndef NONLOCALNEW_H
#define NONLOCALNEW_H
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

#include <array>
#include <unordered_map>
#include <vector>

namespace hamilt
{

#ifndef __NONLOCALNEWTEMPLATE
#define __NONLOCALNEWTEMPLATE

/// The Nonlocal class template inherits from class T
/// it is used to calculate the non-local pseudopotential of wavefunction basis
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK> or OperatorPW<TK>
/// - TR: data type of real space Hamiltonian, it would be double or std::complex<double>
template <class T>
class Nonlocal : public T
{
};

#endif

/// Nonlocal class template specialization for OperatorLCAO<TK> base class
/// It is used to calculate the non-local pseudopotential matrix in real space and fold it to k-space
/// HR = <psi_{mu, 0}|beta_p1>D_{p1, p2}<beta_p2|psi_{nu, R}>
/// HK = <psi_{mu, k}|beta_p1>D_{p1, p2}<beta_p2|psi_{nu, k}> = \sum_{R} e^{ikR} HR
/// Template parameters:
/// - TK: data type of k-space Hamiltonian
/// - TR: data type of real space Hamiltonian
template <typename TK, typename TR>
class Nonlocal<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    Nonlocal<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                      const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                      hamilt::HContainer<TR>* hR_in,
                                      const UnitCell* ucell_in,
                                      const std::vector<double>& orb_cutoff,
                                      const Grid_Driver* GridD_in,
                                      const TwoCenterIntegrator* intor);
    ~Nonlocal<OperatorLCAO<TK, TR>>();

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|beta_p1>D_{p1, p2}<beta_p2|phi_{\nu, R}>
     */
    virtual void contributeHR() override;

    void cal_force_stress(const bool cal_force,
                          const bool cal_stress,
                          const HContainer<TR>* dmR,
                          ModuleBase::matrix& force,
                          ModuleBase::matrix& stress);

    /**
 * @brief Calculate the derivative of the non-local pseudopotential Hamiltonian
 *        with respect to atomic displacements: dH^{NL}/dτ_I
 *
 * For each atom I and each Cartesian direction d ∈ {0,1,2} (x,y,z),
 * fills an HContainer matrix dhR[d][I] that stores
 * dH^{NL}_{IJ}/dτ_I for all atom-pair blocks (I,J,R).
 *
 * ### Mathematical Background
 *
 * Kleinman-Bylander non-local pseudopotential:
 * \f[
 *   V^{NL}(r, r') = \sum_{I} \sum_{m,m'} \beta_{m}^{I}(r)\, D_{mm'}^{I}\, \beta_{m'}^{I}(r')
 * \f]
 *
 * In LCAO basis:
 * \f[
 *   H_{ij}^{NL} = \sum_{m,m'} \langle \phi_i | \beta_{m}^{I} \rangle
 *                            D_{mm'}^{I}
 *                            \langle \beta_{m'}^{I} | \phi_j \rangle
 * \f]
 *
 * ### Derivative
 * \f[
 * \frac{d H_{ij}^{NL}}{d \tau_I}
 *   = \sum_{m,m'} \Bigg[
 *       \langle \frac{d \phi_i}{d \tau_I} | \beta_{m}^{I} \rangle D_{mm'}^{I}
 *           \langle \beta_{m'}^{I} | \phi_j \rangle
 *     + \langle \phi_i | \frac{d \beta_{m}^{I}}{d \tau_I} \rangle D_{mm'}^{I}
 *           \langle \beta_{m'}^{I} | \phi_j \rangle
 *     + \langle \phi_i | \beta_{m}^{I} \rangle D_{mm'}^{I}
 *           \langle \beta_{m'}^{I} | \frac{d \phi_j}{d \tau_I} \rangle
 *     \Bigg]
 * \f]
 *
 * ### Hellmann-Feynman relation
 * \f[
 *   \langle \phi | \frac{d \beta}{d \tau} \rangle
 *     = - \langle \phi | \nabla \beta \rangle
 *     = \langle \nabla \phi | \beta \rangle
 * \f]
 *
 * ### Computational strategy
 *
 * Intermediate vectors per direction d:
 * \f[
 * tU[d] = \sum_{p_1,p_2}
 *     \langle \nabla \phi_{I1} | \beta \rangle_{p_1, d}
 *     \, D_{p_1, p_2} \,
 *     \langle \beta | \phi_{I2} \rangle_{p_2}
 * \f]
 * \f[
 * tV[d] = \sum_{p_1,p_2}
 *     \langle \phi_{I1} | \beta \rangle_{p_1}
 *     \, D_{p_1, p_2} \,
 *     \langle \beta | \nabla \phi_{I2} \rangle_{p_2, d}
 * \f]
 *
 * Distribution:
 *   - dhR[d][I1] -= tU[d]   (orbital 1 moves)
 *   - dhR[d][I2] -= tV[d]   (orbital 2 moves)
 *   - dhR[d][I0] += tU[d]+tV[d]  (Hellmann-Feynman: nucleus moves)
 *
 * @param[out] dhR  Array of 3 vectors (x,y,z), each nat HContainer pointers.
 */
    void cal_dH(std::array<std::vector<hamilt::HContainer<double>*>, 3>& dhR);

    virtual void set_HR_fixed(void*) override;

  private:
    const UnitCell* ucell = nullptr;

    std::vector<double> orb_cutoff_;

    hamilt::HContainer<TR>* HR_fixed = nullptr;

    // the following variable is introduced temporarily during LCAO refactoring
    const TwoCenterIntegrator* intor_ = nullptr;

    bool allocated = false;

    bool HR_fixed_done = false;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const Grid_Driver* GridD_in);

    /**
     * @brief calculate the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * nearest neighbor atoms don't need to be calculated again
     * loop the atom-pairs in HR and calculate the non-local pseudopotential matrix
     */
    void calculate_HR();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const int& T0,
                    const Parallel_Orbitals* paraV,
                    const std::unordered_map<int, std::vector<double>>& nlm1_all,
                    const std::unordered_map<int, std::vector<double>>& nlm2_all,
                    TR* data_pointer);

    const Grid_Driver* gridD = nullptr;

    /**
     * @brief calculate the atomic Force of <I,J,R> atom pair
     */
    void cal_force_IJR(const int& iat1,
                       const int& iat2,
                       const int& T0,
                       const Parallel_Orbitals* paraV,
                       const std::unordered_map<int, std::vector<double>>& nlm1_all,
                       const std::unordered_map<int, std::vector<double>>& nlm2_all,
                       const hamilt::BaseMatrix<TR>* dmR_pointer,
                       double* force1,
                       double* force2);
    /**
     * @brief calculate the Stress of <I,J,R> atom pair
     */
    void cal_stress_IJR(const int& iat1,
                        const int& iat2,
                        const int& T0,
                        const Parallel_Orbitals* paraV,
                        const std::unordered_map<int, std::vector<double>>& nlm1_all,
                        const std::unordered_map<int, std::vector<double>>& nlm2_all,
                        const hamilt::BaseMatrix<TR>* dmR_pointer,
                        const ModuleBase::Vector3<double>& dis1,
                        const ModuleBase::Vector3<double>& dis2,
                        double* stress);

    std::vector<AdjacentAtomInfo> adjs_all;
};

} // namespace hamilt
#endif
