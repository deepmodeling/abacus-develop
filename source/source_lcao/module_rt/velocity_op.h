#ifndef TD_VELOCITY_OP_H
#define TD_VELOCITY_OP_H
#include <unordered_map>
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_base/vector3.h"
#include "source_io/cal_r_overlap_R.h"

//design to calculate velocity operator
template <typename TR>
class Velocity_op
{
  public:
    Velocity_op(const UnitCell* ucell_in,
               const Grid_Driver* GridD_in,
               const Parallel_Orbitals* paraV,
               const LCAO_Orbitals& orb,
               const TwoCenterIntegrator* intor);
    ~Velocity_op();

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

    /// @brief Store real space hamiltonian. TD term should include imaginary part, thus it has to be complex type. Only shared between TD operators.
    std::vector<hamilt::HContainer<std::complex<double>>*> current_term = {nullptr, nullptr, nullptr};
    
    const TwoCenterIntegrator* intor_ = nullptr;
    const TwoCenterIntegrator* intorbeta_ = nullptr;

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
                         const std::vector<std::unordered_map<int, std::vector<double>>>& nlm1_all,
                         const std::vector<std::unordered_map<int, std::vector<double>>>& nlm2_all,
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
    static cal_r_overlap_R r_calculator;
    static bool init_done;
};


#endif // TD_CURRENT_H
