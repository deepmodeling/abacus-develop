#ifndef TD_CURRENT_H
#define TD_CURRENT_H
#include <unordered_map>
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "td_velocity.h"
#include "module_base/vector3.h"

#ifdef __LCAO
//design to calculate current for length gague
class TD_current
{
  public:
    TD_current(const UnitCell* ucell_in,
               const Grid_Driver* GridD_in,
               const Parallel_Orbitals* paraV,
               const LCAO_Orbitals& orb,
               const TwoCenterIntegrator* intor);
    ~TD_current();

    hamilt::HContainer<std::complex<double>>* get_current_term_pointer(const int& i)const 
    {
        return this->current_term[i];
    }

    void calculate_vcomm_r();
    void calculate_grad_term();

  private:
    const UnitCell* ucell = nullptr;

    const Parallel_Orbitals* paraV = nullptr;

    const LCAO_Orbitals& orb_;

    const Grid_Driver* Grid = nullptr;
    /// @brief Store real space hamiltonian. TD term should include imaginary part, thus it has to be complex type. Only shared between TD operators.
    std::vector<hamilt::HContainer<std::complex<double>>*> current_term = {nullptr, nullptr, nullptr};
    
    const TwoCenterIntegrator* intor_ = nullptr;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_vcomm_r(const Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);
    void initialize_grad_term(const Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_vcomm_r_IJR(const int& iat1,
                         const int& iat2,
                         const int& T0,
                         const Parallel_Orbitals* paraV,
                         const std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& nlm1_all,
                         const std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& nlm2_all,
                         std::complex<double>** current_mat_p);
    void cal_grad_IJR(const int& iat1,
                      const int& iat2,
                      const Parallel_Orbitals* paraV,
                      const ModuleBase::Vector3<double>& dtau,
                      std::complex<double>** current_mat_p);

    /// @brief exact the nearest neighbor atoms from all adjacent atoms
    std::vector<AdjacentAtomInfo> adjs_vcommr;
    std::vector<AdjacentAtomInfo> adjs_grad;

    /// @brief Store the vector potential for td_ekinetic term
    ModuleBase::Vector3<double> cart_At;
};


#endif // __LCAO
#endif // TD_CURRENT_H
