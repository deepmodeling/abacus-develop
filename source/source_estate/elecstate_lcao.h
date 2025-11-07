#ifndef ELECSTATE_LCAO_H
#define ELECSTATE_LCAO_H

#include "elecstate.h"
#include "source_estate/module_dm/density_matrix.h"

#include <vector>

namespace elecstate
{
template <typename TK>
class ElecStateLCAO : public ElecState
{
  public:
    ElecStateLCAO()
    {
    } // will be called by ElecStateLCAO_TDDFT
    ElecStateLCAO(Charge* chr_in,
                  const K_Vectors* klist_in,
                  int nks_in,
                  ModulePW::PW_Basis* rhopw_in,
                  ModulePW::PW_Basis_Big* bigpw_in)
    {
        init_ks(chr_in, klist_in, nks_in, rhopw_in, bigpw_in);
        this->classname = "ElecStateLCAO";
    }

    virtual ~ElecStateLCAO()
    {
    }

    // update charge density for next scf step
    // void getNewRho() override;

    static int out_wfc_lcao;
    static bool need_psi_grid;

    double get_spin_constrain_energy() override;

    // use for pexsi

    /**
     * @brief calculate electronic charge density from pointers of density matrix calculated by pexsi
     * @param pexsi_DM: pointers of density matrix (DMK) calculated by pexsi
     * @param pexsi_EDM: pointers of energy-weighed density matrix (EDMK) calculated by pexsi, needed by MD, will be
     * stored in DensityMatrix::pexsi_EDM
     */
	void dm2rho(std::vector<TK*> pexsi_DM, 
			std::vector<TK*> pexsi_EDM, 
			DensityMatrix<TK, double>* dm);

};

template <typename TK>
int ElecStateLCAO<TK>::out_wfc_lcao = 0;

template <typename TK>
bool ElecStateLCAO<TK>::need_psi_grid = true;

} // namespace elecstate

#endif
