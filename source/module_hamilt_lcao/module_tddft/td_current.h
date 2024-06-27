#ifndef TD_CURRENT_H
#define TD_CURRENT_H
#include <unordered_map>
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_basis/module_ao/ORB_table_phi.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_base/vector3.h"

#ifdef __LCAO

class TD_current
{
  public:
    TD_current(const UnitCell* ucell_in,
               Grid_Driver* GridD_in,
               const Parallel_Orbitals* paraV,
               const ModuleBase::Vector3<double>& cart_At,
               const TwoCenterIntegrator* intor);
    ~TD_current();

    hamilt::HContainer<std::complex<double>>* get_current_term_pointer(const int& i)const 
    {
        return this->current_term[i];
    }

    void initialize_current_term(const hamilt::HContainer<std::complex<double>>* HR, const Parallel_Orbitals* paraV);
    void calculate_vcomm_r();
    void calculate_grad_term();

    /// @brief calculate specific current term or not
    static bool cal_grad;
    static bool cal_vcomm_r;

  private:
    const UnitCell* ucell = nullptr;

    const Parallel_Orbitals* paraV = nullptr;

    /// @brief Store real space hamiltonian. TD term should include imaginary part, thus it has to be complex type. Only shared between TD operators.
    std::vector<hamilt::HContainer<std::complex<double>>*> current_term = {nullptr, nullptr, nullptr};
    Grid_Driver* Grid = nullptr;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_vcomm_r(Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);
    void initialize_grad_term(Grid_Driver* GridD_in, const Parallel_Orbitals* paraV);


    //calculate term [r,Vnl]
    void cal_vcomm_r_IJR(const int& iat1,
                    const int& iat2,
                    const int& T0,
                    const Parallel_Orbitals* paraV,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm1_all,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm1x_all,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm1y_all,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm1z_all,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm2_all,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm2x_all,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm2y_all,
                    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm2z_all,
                    std::complex<double>* data_pointerx,
                    std::complex<double>* data_pointery,
                    std::complex<double>* data_pointerz);
    //calculate term -i⋅∇
    void cal_grad_IJR(const int& iat1,
                      const int& iat2,
                      const Parallel_Orbitals* paraV,
                      const ModuleBase::Vector3<double>& dtau,
                      std::complex<double>* data_pointerx,
                      std::complex<double>* data_pointery,
                      std::complex<double>* data_pointerz);

    /// @brief exact the nearest neighbor atoms from all adjacent atoms
    std::vector<AdjacentAtomInfo> adjs_vcommr;
    std::vector<AdjacentAtomInfo> adjs_grad;
    /// @brief Store the vector potential for td_nonlocal term
    ModuleBase::Vector3<double> cart_At;

    const TwoCenterIntegrator* intor_ = nullptr;

};


#endif // __LCAO
#endif // TD_CURRENT_H