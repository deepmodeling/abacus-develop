#ifndef TD_VELOCITY_OP_H
#define TD_VELOCITY_OP_H

#include <array>
#include <memory>
#include <unordered_map>

#include "source_base/atom_basis_view.h"
#include "source_base/basis_flattened_indexing.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_base/vector3.h"
#include "source_io/module_hs/pos_op_mat.h"

class BlockCyclicDistribution; // Defined in source/source_base/parallel_blkcyc.h

// Cartesian velocity-operator terms for RT-TDDFT currents.
template <typename TR>
class Velocity_op
{
  public:
    Velocity_op(const UnitCell* ucell_in,
               const Grid_Driver* GridD_in,
               const BasisFlattenedIndexing& basis,
               const BlockCyclicDistribution& layout,
               const LCAO_Orbitals& orb,
               const TwoCenterIntegrator* intor);

    hamilt::HContainer<std::complex<double>>* get_current_term_pointer(const int& i) const
    {
        return this->current_term[i].get();
    }
    void calculate_vcomm_r();
    void calculate_grad_term();

  private:
    const UnitCell* ucell = nullptr;

    std::unique_ptr<BasisFlattenedIndexing> basis_;

    const BlockCyclicDistribution& layout_;

    const LCAO_Orbitals& orb_;

    /// @brief Cartesian current terms. TD terms include an imaginary part, so they are complex.
    std::array<std::unique_ptr<hamilt::HContainer<std::complex<double>>>, 3> current_term;

    // Each view borrows TD_info's basis indexing and distribution, which are
    // const references and therefore never reseat, so these stay valid for the
    // whole lifetime of this object.
    std::vector<AtomBasisView> atom_rows_;
    std::vector<AtomBasisView> atom_cols_;

    const TwoCenterIntegrator* intor_ = nullptr;

    /**
     * @brief initialize HR, search the nearest neighbor atoms
     * HContainer is used to store the non-local pseudopotential matrix with specific <I,J,R> atom-pairs
     * the size of HR will be fixed after initialization
     */
    void initialize_vcomm_r(const Grid_Driver* GridD_in);
    void initialize_grad_term(const Grid_Driver* GridD_in);

    void build_atom_basis();

    /**
     * @brief calculate the HR local matrix of <I,J,R> atom pair
     */
    void cal_vcomm_r_IJR(const int& iat1,
                         const int& iat2,
                         const int& T0,
                         const std::vector<std::unordered_map<int, std::vector<double>>>& nlm1_all,
                         const std::vector<std::unordered_map<int, std::vector<double>>>& nlm2_all,
                         std::complex<double>** current_mat_p);
    void cal_grad_IJR(const int& iat1,
                      const int& iat2,
                      const ModuleBase::Vector3<double>& dtau,
                      std::complex<double>** current_mat_p);

    /// @brief exact the nearest neighbor atoms from all adjacent atoms
    std::vector<AdjacentAtomInfo> adjs_vcommr;
    std::vector<AdjacentAtomInfo> adjs_grad;

    /// @brief Store the vector potential for td_ekinetic term
    ModuleBase::Vector3<double> cart_At;
    static Position_op r_calculator;
    static bool init_done;
};


#endif // TD_VELOCITY_OP_H
