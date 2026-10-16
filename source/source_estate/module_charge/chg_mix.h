#ifndef CHG_MIX_H
#define CHG_MIX_H
#include "charge.h"
#include "chg_mix_cfg.h"
#include "source_base/module_mixing/mixing.h"
#include "source_base/module_mixing/plain_mixing.h"
#include <functional>
#include <memory>

class Charge_Mixing
{
  /// Charge_Mixing class
  /// This class is used to mix charge density, kinetic energy density and real-space density matrix
  /// This Charge_Mixing class offers the following interfaces:
  /// 1. set_mixing() to set all private mixing parameters
  /// 2. init_mixing() to initialize mixing, including allocating memory for mixing data and reset mixing
  /// 3. mix_rho() to mix charge density
  /// Real-space density matrix mixing is implemented by the stateless
  /// module_charge::init_mixing_dmr/mix_dmr functions in chg_dmr.h; this class
  /// only owns the mixing history buffer, exposed through get_dmr_mdata().
  /// how to use it:
  /// you can (re)start a mixing by calling set_mixing() and init_mixing() before calling mix_rho()

  public:
    Charge_Mixing();
    ~Charge_Mixing();

    /**
     * @brief Set all private mixing parameters from an aggregated config
     * @param cfg mixing parameters and runtime globals (nspin, scf_thr_type, double_grid)
     * @param rhopw_in smooth grid
     * @param rhodpw_in dense grid when double grid is used, otherwise same as rhopw
     * @param omega_in omega for non-linear core correction
     * @param tpiba_in 2*pi/beta for non-linear core correction
     */
    void set_mixing(const MixingConfig& cfg,
                    ModulePW::PW_Basis* rhopw_in,
                    ModulePW::PW_Basis* rhodpw_in,
                    double& omega_in,
                    double& tpiba_in);

    /// Disable Kerker screening for subsequent mix_rho calls.
    /// Used by the non-separate-loop EXX path (exx_lri_interface.hpp)
    /// after EXX convergence: Kerker damping fights the DM update there.
    /// The Kerker kernels read cfg_ (immutable INPUT snapshot), so the
    /// disable flag must live on Charge_Mixing itself rather than mutating cfg_.
    void close_kerker_gg0() { kerker_disabled_ = true; }
    /**
     * @brief initialize mixing, including constructing mixing and allocating memory for mixing data
     * @brief this function should be called at eachiterinit()
     */
    void init_mixing();

    /**
     * @brief charge mixing
     * @param chr pointer of Charge object
     */
    void mix_rho(Charge* chr);

    /**
     * @brief allocate memory of uom_mdata
     * @param uom_size size of DFT+U occupation matrix
     */
    void allocate_mixing_uom(int size_uom);

    /**
     * @brief DFT+U occupation matrix mixing
     * @param uom_in output occupation matrix
     * @param uom_save_in input occupation matrix
     */
    void mix_uom(std::vector<double>& uom_in, std::vector<double>& uom_save_in);

    /**
     * @brief reset mixing, actually we only call init_mixing() to reset mixing instead of this function
     */
    void mix_reset();
    
    // extracting parameters normally these parameters will not be used outside charge mixing
    // while Exx is using them as well as some other places
    const std::string& get_mixing_mode() const {return cfg_.mixing_mode;}
    double get_mixing_beta() const {return cfg_.mixing_beta;}
    int get_mixing_ndim() const {return cfg_.mixing_ndim;}
    Base_Mixing::Mixing* get_mixing() const {return mixing.get();}

    /**
     * @brief mutable access to the real-space density-matrix mixing history
     *
     * The history buffer is owned by Charge_Mixing and driven by the
     * stateless module_charge::init_mixing_dmr/mix_dmr functions in chg_dmr.h.
     */
    Base_Mixing::Mixing_Data& get_dmr_mdata() {return dmr_mdata;}

    /**
     * @brief mutable access to the charge-density mixing history
     *
     * Same ownership story as get_dmr_mdata(): the buffer belongs to
     * Charge_Mixing but is driven through Base_Mixing::Mixing, so callers
     * that reset or inspect the history need a handle on it.
     */
    Base_Mixing::Mixing_Data& get_rho_mdata() {return rho_mdata;}

    /**
     * @brief mutable access to the kinetic-energy-density mixing history
     */
    Base_Mixing::Mixing_Data& get_tau_mdata() {return tau_mdata;}

    /**
     * @brief read-only access to the aggregated mixing config set by set_mixing()
     */
    const MixingConfig& get_mixing_config() const {return cfg_;}

    /**
     * @brief replace the aggregated mixing config without re-running set_mixing()
     *
     * set_mixing() both stores the config and rebuilds the mixing objects. This
     * setter is for callers that need to update the snapshot alone - typically
     * because a derived runtime global (nspin, domag, ...) changed after the
     * mixing objects were already built - and must not disturb the mixing
     * history. Pair it with get_mixing_config() to amend individual fields.
     */
    void set_mixing_config(const MixingConfig& cfg) {cfg_ = cfg;}

    // for mixing restart
    /// which step to restart mixing during SCF
    int mixing_restart_step = 0;
    /// the number of restart mixing during SCF
    int mixing_restart_count = 0;
    /// the label of mixing restart step
    int mixing_restart_last = 0;

    // to calculate the slope of drho curve during SCF, which is used to determine if SCF oscillate
    bool if_scf_oscillate(const int iteration, const double drho,
                           const int iternum_used, const double threshold);

  private:

    // mixing_data
    /// Mixing object for charge and kinetic energy
    std::unique_ptr<Base_Mixing::Mixing> mixing;
    Base_Mixing::Mixing_Data rho_mdata;    ///< Mixing data for charge density
    Base_Mixing::Mixing_Data tau_mdata;    ///< Mixing data for kinetic energy density
    Base_Mixing::Mixing_Data dmr_mdata;    ///< Mixing data for real space density matrix
    Base_Mixing::Mixing_Data uom_mdata;    ///< Mixing data for DFT+U occupation matrix
    std::unique_ptr<Base_Mixing::Plain_Mixing> mixing_highf; ///< The high_frequency part is mixed by plain mixing method.

    //======================================
    // private mixing parameters
    //======================================
    MixingConfig cfg_;                 ///< aggregated mixing config, also holds nspin/scf_thr_type/double_grid
    double* omega = nullptr;                  ///< omega for non-linear core correction
    double* tpiba = nullptr;                  ///< 2*pi/beta for non-linear core correction
    std::vector<double> _drho_history; ///< history of drho used to determine the oscillation, size is scf_nmax

    ModulePW::PW_Basis* rhopw = nullptr;  ///< smooth grid
    ModulePW::PW_Basis* rhodpw = nullptr; ///< dense grid, same as rhopw for ncpp.

    /// Runtime override set by close_kerker_gg0(): short-circuits the
    /// Kerker screening lambdas in mix_rho_recip/mix_rho_real so the
    /// non-separate-loop EXX path can disable Kerker after convergence.
    /// Lives here, not in MixingConfig, because cfg_ is an immutable
    /// INPUT snapshot consumed by the stateless Kerker kernels.
    bool kerker_disabled_ = false;

    /**
     * @brief charge mixing for reciprocal space
     * @param chr pointer of Charge object
     */
    void mix_rho_recip(Charge* chr);

    /**
     * @brief charge mixing for real space
     * @param chr pointer of Charge object
     */
    void mix_rho_real(Charge* chr);
};

#endif
