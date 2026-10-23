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

    double get_spin_constrain_energy() override;

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

    /**
     * @brief calculate electronic charge density from the density matrix (DMR)
     *
     * Thin wrapper over LCAO_domain::dm2rho so that HSolverLCAO delegates the
     * charge-density calculation through the ElecState interface, mirroring the
     * plane-wave path (ElecStatePW::psiToRho) and the pexsi branch above. This
     * keeps the source_lcao dependency out of source_hsolver.
     *
     * @param omega current unit-cell volume (ucell.omega). Must not be
     *        rhopw->omega, which is stale in variable-cell calculations.
     */
    void dmToRho(std::vector<hamilt::HContainer<double>*>& dmr,
                 int nspin,
                 Charge* chr,
                 const double omega,
                 bool skip_charge = false);

};

template <typename TK>
int ElecStateLCAO<TK>::out_wfc_lcao = 0;

template <typename TK>
bool ElecStateLCAO<TK>::need_psi_grid = true;

} // namespace elecstate

#endif
