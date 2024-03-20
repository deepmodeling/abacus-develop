#ifndef DEEPKSLCAO_H
#define DEEPKSLCAO_H
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "operator_lcao.h"
#include "module_elecstate/module_dm/density_matrix.h"

namespace hamilt
{

#ifndef __DEEPKSTEMPLATE
#define __DEEPKSTEMPLATE

/// The DeePKS class template inherits from class T
/// it is used to calculate the Deep Potential Kohn-Sham correction from DeePKS method
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK>
/// - TR: data type of real space Hamiltonian, it would be double or std::complex<double>
template <class T>
class DeePKS : public T
{
};

#endif

template <typename TK, typename TR>
class DeePKS<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    DeePKS<OperatorLCAO<TK, TR>>(Local_Orbital_Charge* loc_in,
                            LCAO_Matrix* LM_in,
                            const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                            HContainer<TR>* hR_in,
                            std::vector<TK>* hK_in,
                            const UnitCell* ucell_in,
                            Grid_Driver* GridD_in,
                            const int& nks_in,
                            elecstate::DensityMatrix<TK,double>* DM_in);
    ~DeePKS();

    /**
     * @brief contribute the DeePKS correction to real space Hamiltonian
     * this function is used for update hR and H_V_delta
    */
    virtual void contributeHR() override;
#ifdef __DEEPKS
    /**
     * @brief contribute the DeePKS correction for each k-point to ld.H_V_delta or ld.H_V_delta_k
     * this function is not used for update hK, but for DeePKS module
     * @param ik: the index of k-point
    */
    virtual void contributeHk(int ik) override;
#endif

  private:
    Local_Orbital_Charge* loc;

    elecstate::DensityMatrix<TK,double>* DM;

    const UnitCell* ucell = nullptr;

    HContainer<TR>* H_V_delta = nullptr;
#ifdef __DEEPKS

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the DeePKS real space Hamiltonian correction with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(Grid_Driver* GridD, const Parallel_Orbitals* paraV);

    /**
     * @brief calculate the DeePKS correction matrix with specific <I,J,R> atom-pairs
     * use the adjs_all to calculate the HR matrix
     */
    void calculate_HR();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const double* hr_in,
        const int& row_size,
        const int& col_size,
        TR* data_pointer);

    void pre_calculate_nlm(const int iat0, std::vector<std::unordered_map<int, std::vector<double>>>& nlm_in);
    std::vector<std::vector<std::unordered_map<int, std::vector<double>>>> nlm_tot;
    /**
     * @brief initialize H_V_delta, search the nearest neighbor atoms
     * used for calculate the DeePKS real space Hamiltonian correction with specific <I,J,R> atom-pairs
    */
    std::vector<AdjacentAtomInfo> adjs_all;
#endif
    const int& nks;
};

} // namespace hamilt
#endif