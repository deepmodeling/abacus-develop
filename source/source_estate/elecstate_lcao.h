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
                  ModulePW::PW_Basis_Big* bigpw_in)
    {
        init_ks(chr_in, klist_in, nks_in, bigpw_in);
        this->classname = "ElecStateLCAO";
    }

    virtual ~ElecStateLCAO() = default;

    static int out_wfc_lcao;
    static bool need_psi_grid;

    // use for pexsi

    /**
     * @brief calculate electronic charge density from pointers of density matrix calculated by pexsi
     * @param pexsi_DM: pointers of density matrix (DMK) calculated by pexsi
     * @param edm_pexsi: pointers of energy-weighed density matrix (edmk) calculated by pexsi, needed by MD, will be
     * stored in DensityMatrix::edm_pexsi
     */
	void dm2rho(std::vector<TK*> pexsi_DM,
			std::vector<TK*> edm_pexsi,
			module_dm::DensityMatrix<TK, double>* dm,
			const double omega);

};

template <typename TK>
int ElecStateLCAO<TK>::out_wfc_lcao = 0;

template <typename TK>
bool ElecStateLCAO<TK>::need_psi_grid = true;

} // namespace elecstate

#endif
