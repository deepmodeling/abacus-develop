#ifndef CHG_EXTRA_H
#define CHG_EXTRA_H

#include <vector>

#include "charge.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/stru_fac.h"
#ifdef __MPI
#include "source_base/parallel_grid.h"
#endif

namespace module_charge
{

/// Forward declaration: see chg_atomic.h for the full definition.
struct AtomicRhoCfg;

/**
 * @brief charge extrapolation method
 *
 * This class implements several charge extrapolation methods:
 * pot_order=0 : copy the old potential (nothing is done);
 * pot_order=1 : subtract old atomic charge density and sum the new
 *               if dynamics is done the routine extrapolates also the difference
 *               between the scf charge and the atomic one;
 * pot_order=2 : first order extrapolation:
 *                              \[ \rho(t+dt) = 2\ \rho(t)-\rho(t-dt); \]
 * pot_order=3 : second order extrapolation:
 *                               \[ \rho(t+dt) = \rho(t) + \alpha_0\ (\rho(t) - \rho(t-dt))
 *                               + \beta_0\ (\rho(t-dt)- \rho(t-2 dt)). \]
 *
 * The \(\alpha_0\) and \(\beta_0\) parameters are calculated in find_alpha_and_beta()
 * so that \(|\tau'-\tau(t+dt)|\) is minimum. \(\tau'\) and \(\tau(t+dt)\) are respectively
 *  the atomic positions at time t+dt and the extrapolated one:
 *  \[ \tau(t+dt) = \tau(t) + \alpha_0\ ( \tau(t)    - \tau(t-dt)   )
 *                          + \beta_0\ ( \tau(t-dt) - \tau(t-2 dt) ). \]
 */

class Charge_Extra
{
    public:

    Charge_Extra();

    // The displacement and delta-rho arrays record the history of one
    // molecular-dynamics trajectory; a copy would silently fork that
    // history. Forbid copies until a deep copy is needed.
    Charge_Extra(const Charge_Extra&) = delete;
    Charge_Extra& operator=(const Charge_Extra&) = delete;

    /**
     * @brief Initialization of viriables used in charge extrapolation methods
     *
     * When Esolver is initialized, ucell.natom is not determined
     * As a result, data structures in Charge_Extra cannot be allocated
     * This is a temporary solution by delaying the allocation
     * But after ucell and Esolver are fully decoupled
     * Init_CE will be removed and everything put back in the constructor
     *
     * @param nspin the number of spins
     * @param natom the number of atoms
     * @param nrxx the number of grids
     * @param chg_extrap the charge extrapolation method
     */
    void Init_CE(const int& nspin, const int& natom, const int& nrxx, const std::string chg_extrap);

    /**
     * @brief charge extrapolation method
     *
     * @param Pgrid parallel grids
     * @param ucell the cell information
     * @param chr the charge density
     * @param rhopw plane-wave basis bound to chr
     * @param sf the structure factor
     * @param ofs_running the output stream
     * @param ofs_warning the output stream
     * @param atomic_rho_cfg configuration for atomic_rho (nelec, magnetism, verbosity)
     * @param has_float_data whether the Hamiltonian carries float data; forwarded
     *        to Structure_Factor::setup() when the structure factor is refreshed
     */
    void extrapolate_charge(
        Parallel_Grid* Pgrid,
        UnitCell& ucell,
        Charge* chr,
        const ModulePW::PW_Basis& rhopw,
        Structure_Factor* sf,
        std::ofstream& ofs_running,
        std::ofstream& ofs_warning,
        const AtomicRhoCfg& atomic_rho_cfg,
        const bool has_float_data);

    /**
     * @brief update displacements
     *
     * In the second order extrapolation, the displacements of previous three steps are needed to determine alpha and
     * beta, which are parameters used in this method.
     *
     * @param ucell the cell information
     */
    void update_all_dis(const UnitCell& ucell);

    /**
     * @brief update the difference of charge density
     *
     * @param ucell the cell information
     * @param chr the charge density
     * @param rhopw plane-wave basis bound to chr
     * @param sf the structure factor
     * @param atomic_rho_cfg configuration for atomic_rho (nelec, magnetism, verbosity)
     */
    void update_delta_rho(const UnitCell& ucell,
                          const Charge* chr,
                          const ModulePW::PW_Basis& rhopw,
                          const Structure_Factor* sf,
                          const AtomicRhoCfg& atomic_rho_cfg);

    //==========================================================
    // Accessors over the extrapolation state. The history buffers and the
    // step bookkeeping are driven by Init_CE()/extrapolate_charge(); these
    // expose them for inspection and for seeding a known state, which is
    // what lets callers exercise one extrapolation order at a time.
    //==========================================================

    /// @brief the current step
    int get_istep() const { return istep; }
    void set_istep(const int istep_in) { istep = istep_in; }

    /// @brief the specified charge extrapolation method
    int get_pot_order() const { return pot_order; }
    void set_pot_order(const int pot_order_in) { pot_order = pot_order_in; }

    /// @brief the extrapolation order actually used this step, which
    ///        extrapolate_charge() derives from pot_order and istep
    int get_rho_extr() const { return rho_extr; }

    /// @brief the number of spins
    int get_nspin() const { return nspin; }

    /// @brief displacement histories: dis_old2 = pos_old1 - pos_old2,
    ///        dis_old1 = pos_now - pos_old1, dis_now = pos_next - pos_now
    std::vector<ModuleBase::Vector3<double>>& get_dis_old1() { return dis_old1; }
    std::vector<ModuleBase::Vector3<double>>& get_dis_old2() { return dis_old2; }
    std::vector<ModuleBase::Vector3<double>>& get_dis_now() { return dis_now; }

    /// @brief the last three steps' difference of rho and atomic_rho
    std::vector<std::vector<double>>& get_delta_rho1() { return delta_rho1; }
    std::vector<std::vector<double>>& get_delta_rho2() { return delta_rho2; }
    std::vector<std::vector<double>>& get_delta_rho3() { return delta_rho3; }

    /// @brief parameter used in the second order extrapolation
    double get_alpha() const { return alpha; }

    /// @brief coefficient for the third delta_rho contribution
    double get_beta() const { return beta; }

    /// @brief run find_alpha_and_beta() against the displacement histories
    ///        currently held, so its solution can be checked on its own
    ///        rather than only through a full extrapolate_charge() step
    void find_alpha_and_beta_for_testing(const int& natom,
                                         std::ofstream& ofs_running,
                                         std::ofstream& ofs_warning)
    {
        find_alpha_and_beta(natom, ofs_running, ofs_warning);
    }

  private:
    int istep = 0; ///< the current step
    int pot_order; ///< the specified charge extrapolation method
    int rho_extr;  ///< the actually used method
    int nspin;        ///< the number of spins

    std::vector<ModuleBase::Vector3<double>> dis_old1; ///< dis_old2 = pos_old1 - pos_old2
    std::vector<ModuleBase::Vector3<double>> dis_old2; ///< dis_old1 = pos_now - pos_old1
    std::vector<ModuleBase::Vector3<double>> dis_now;  ///< dis_now = pos_next - pos_now

    std::vector<std::vector<double>> delta_rho1; ///< the last step difference of rho and atomic_rho
    std::vector<std::vector<double>> delta_rho2; ///< the second last step difference of rho and atomic_rho
    std::vector<std::vector<double>> delta_rho3; ///< the third last step difference of rho and atomic_rho

    double alpha; ///< parameter used in the second order extrapolation
    double beta;  ///< coefficient for the third delta_rho contribution

    /**
     * @brief determine alpha and beta
     *
     * @param natom the number of atoms
     * @param ofs_running the output stream
     * @param ofs_warning the output stream
     */
    void find_alpha_and_beta(const int& natom, std::ofstream& ofs_running, std::ofstream& ofs_warning);
};

} // namespace module_charge

#endif
