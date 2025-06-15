#ifndef DEEPKSLCAO_H
#define DEEPKSLCAO_H
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __MLALGOTEMPLATE
#define __MLALGOTEMPLATE

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
    DeePKS<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                 const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                 HContainer<TR>* hR_in,
                                 const UnitCell* ucell_in,
                                 const Grid_Driver* GridD_in,
                                 const TwoCenterIntegrator* intor_orb_alpha,
                                 const LCAO_Orbitals* ptr_orb,
                                 const int& nks_in,
                                 elecstate::DensityMatrix<TK, double>* DM_in
#ifdef __MLALGO
                                 ,
                                 LCAO_Deepks<TK>* ld_in
#endif
    );
    ~DeePKS();

    /**
     * @brief contribute the DeePKS correction to real space Hamiltonian
     * this function is used for update hR and V_delta_R
     */
    virtual void contributeHR() override;
#ifdef __MLALGO
    /**
     * @brief contribute the DeePKS correction for each k-point to V_delta
     * this function is not used for update hK, but for DeePKS module
     * @param ik: the index of k-point
     */
    virtual void contributeHk(int ik) override;

    HContainer<TR>* get_V_delta_R() const
    {
        return this->V_delta_R;
    }
#endif

  private:
    elecstate::DensityMatrix<TK, double>* DM;

    const UnitCell* ucell = nullptr;
    Grid_Driver* gridD = nullptr;

    const Grid_Driver* gd = nullptr;

    HContainer<TR>* V_delta_R = nullptr;

    // the following variable is introduced temporarily during LCAO refactoring
    const TwoCenterIntegrator* intor_orb_alpha_ = nullptr;
    const LCAO_Orbitals* ptr_orb_ = nullptr;

#ifdef __MLALGO

    LCAO_Deepks<TK>* ld = nullptr;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the DeePKS real space Hamiltonian correction with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_HR(const Grid_Driver* GridD);

    /**
     * @brief calculate the DeePKS correction matrix with specific <I,J,R> atom-pairs
     * use the adjs_all to calculate the HR matrix
     */
    void calculate_HR();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_HR_IJR(const double* hr_in, const int& row_size, const int& col_size, TR* data_pointer);

    /**
     * @brief initialize V_delta_R, search the nearest neighbor atoms
     * used for calculate the DeePKS real space Hamiltonian correction with specific <I,J,R> atom-pairs
     */
    std::vector<AdjacentAtomInfo> adjs_all;
#endif
    const int& nks;
};

} // namespace hamilt
#endif
