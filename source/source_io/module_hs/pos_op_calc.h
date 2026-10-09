#ifndef POS_OP_CALC_H
#define POS_OP_CALC_H

#include "source_base/vector3.h"
#include "source_cell/unitcell.h"
#include "pos_op_basis.h"

#include <vector>

/**
 * @brief Evaluate position-operator matrix elements using the tables
 *        prepared by PosOpBasis.
 *
 * All methods are pure numerical integrations; they do not perform I/O.
 */
class PosOpCalc
{
  public:
    explicit PosOpCalc(const PosOpBasis& basis);

    /**
     * @brief <phi_mu | r_hat | phi_nu> for a pair of numerical atomic orbitals.
     *        Returns the vector (x, y, z) components.
     */
    ModuleBase::Vector3<double> pos_matrix(const ModuleBase::Vector3<double>& R1,
                                           int T1, int L1, int m1, int N1,
                                           const ModuleBase::Vector3<double>& R2,
                                           int T2, int L2, int m2, int N2) const;

    /**
     * @brief Electric-field coupling term <phi_mu | r_hat | grad phi_nu> * Efield.
     */
    ModuleBase::Vector3<double> pos_grad_matrix(const ModuleBase::Vector3<double>& R1,
                                                int T1, int L1, int m1, int N1,
                                                const ModuleBase::Vector3<double>& R2,
                                                int T2, int L2, int m2, int N2,
                                                const ModuleBase::Vector3<double>& Efield,
                                                const ModuleBase::Vector3<double>& dR) const;

    /**
     * @brief <phi_mu | r_hat | beta> for non-local projectors.
     *        nlm[0] = overlap, nlm[1..3] = position matrix elements.
     */
    void pos_beta_matrix(const UnitCell& ucell,
                         std::vector<std::vector<double>>& nlm,
                         const ModuleBase::Vector3<double>& R1,
                         int T1, int L1, int m1, int N1,
                         const ModuleBase::Vector3<double>& R2,
                         int T2) const;

  private:
    const PosOpBasis& basis_;
};

#endif // POS_OP_CALC_H
