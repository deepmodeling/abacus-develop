#ifndef VEFFLCAO_H
#define VEFFLCAO_H
#include "source_base/timer.h"
#include "source_estate/module_pot/potential_new.h"
#include "operator_lcao.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include <vector>

namespace hamilt
{

#ifndef __VEFFTEMPLATE
#define __VEFFTEMPLATE

template <class T>
class Veff : public T
{
};

#endif

/// @brief Effective potential class, used for calculating Hamiltonian with grid integration tools
/// If user want to separate the contribution of V_{eff} into V_{H} and V_{XC} and V_{local pseudopotential} and so on,
/// the user can separate the Potential class into different parts, and construct different Veff class for each part.
/// @tparam TK 
/// @tparam TR 
template <typename TK, typename TR>
class Veff<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    /**
     * @brief Construct a new Veff object
    */
    Veff<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                               const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                               elecstate::Potential* pot_in,
                               hamilt::HContainer<TR>* hR_in,
                               const UnitCell* ucell_in,
                               const std::vector<double>& orb_cutoff,
                               const Grid_Driver* GridD_in,
                               const int& nspin)
        : orb_cutoff_(orb_cutoff), pot(pot_in), ucell(ucell_in),
          gd(GridD_in), OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in)
    {
        this->cal_type = calculation_type::lcao_gint;

        this->initialize_HR(ucell_in, GridD_in);
    }

    ~Veff<OperatorLCAO<TK, TR>>(){};

    /**
     * @brief contributeHR() is used to calculate the HR matrix
     * <phi_{\mu, 0}|V_{eff}|phi_{\nu, R}>
     * the contribution of V_{eff} is calculated by the contribution of V_{H} and V_{XC} and V_{local pseudopotential} and so on.
     * grid integration is used to calculate the contribution Hamiltonian of effective potential
     */
    virtual void contributeHR() override;
  
  const UnitCell* ucell;
  const Grid_Driver* gd;

private:

  std::vector<double> orb_cutoff_;

  // Charge calculating method in LCAO base and contained grid base calculation: DM_R, DM, pvpR_reduced

  elecstate::Potential* pot = nullptr;

  int nspin = 1;

  /**
   * @brief initialize HR, search the nearest neighbor atoms
   * HContainer is used to store the electronic kinetic matrix with specific <I,J,R> atom-pairs
   * the size of HR will be fixed after initialization
   */
  void initialize_HR(const UnitCell* ucell_in, const Grid_Driver* GridD_in);
};

} // namespace hamilt
#endif
