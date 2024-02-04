
#ifndef CHARGE_MIXING_H
#define CHARGE_MIXING_H
#include "charge.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/module_mixing/mixing.h"
#include "module_base/module_mixing/plain_mixing.h"
#include "module_cell/unitcell.h"
class Charge_Mixing
{
  /// Charge_Mixing class
  /// This class is used to mix charge density, kinetic energy density and real-space density matrix
  /// This Charge_Mixing class offers the following interfaces:
  /// 1. set_mixing() to set all private mixing parameters
  /// 2. init_mixing() to initialize mixing, including allocating memory for mixing data and reset mixing
  /// 3. mix_rho() to mix charge density
  /// 4. mix_dmr() to mix real-space density matrix
  /// how to use it:
  /// you can (re)start a mixing by calling set_mixing() and init_mixing() before calling mix_rho() or mix_dmr()

  public:
    Charge_Mixing();
    ~Charge_Mixing();

    /**
     * @brief reset mixing
     */
    void mix_reset();

    /**
     * @brief charge mixing
     * @param chr pointer of Charge object
     */
    void mix_rho(Charge* chr);

    /**
     * @brief density matrix mixing, only for LCAO
     * @param DM pointer of DensityMatrix object
     */
    void mix_dmr(elecstate::DensityMatrix<double, double>* DM);
    void mix_dmr(elecstate::DensityMatrix<std::complex<double>, double>* DM);

    /**
     * @brief Set all private mixing paramters
     * @param mixing_mode_in mixing mode: "plain", "broyden", "pulay"
     * @param mixing_beta_in mixing beta
     * @param mixing_ndim_in mixing ndim
     * @param mixing_gg0_in mixing gg0 for Kerker screen
     * @param mixing_tau_in whether to use tau mixing
     * @param mixing_beta_mag_in mixing beta for magnetism
     * @param mixing_gg0_mag_in mixing gg0 for Kerker screen for magnetism
     * @param mixing_gg0_min_in minimum kerker coefficient
     * @param mixing_angle_in mixing angle for nspin=4
     * @param mixing_dmr_in whether to mixing real space density matrix
     */
    void set_mixing(const std::string& mixing_mode_in,
                    const double& mixing_beta_in,
                    const int& mixing_ndim_in,
                    const double& mixing_gg0_in,
                    const bool& mixing_tau_in,
                    const double& mixing_beta_mag_in,
                    const double& mixing_gg0_mag_in,
                    const double& mixing_gg0_min_in,
                    const double& mixing_angle_in,
                    const bool& mixing_dmr_in);

    /**
     * @brief allocate memory of dmr_mdata
     * @param nnr size of real-space density matrix
     */
    void allocate_mixing_dmr(int nnr);

    /**
     * @brief Get the drho
     *
     */
    double get_drho(Charge* chr, const double nelec);
    
    /**
     * @brief Set the smooth and dense grids
     * @param rhopw_in smooth grid
     * @param rhodpw_in dense grid when double grid is used, otherwise same as rhopw
     */
    void set_rhopw(ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis* rhodpw_in);

    // extracting parameters normally these parameters will not be used outside charge mixing
    // while Exx is using them as well as some other places
    const std::string& get_mixing_mode() const {return mixing_mode;}
    double get_mixing_beta() const {return mixing_beta;}
    int get_mixing_ndim() const {return mixing_ndim;}
    double get_mixing_gg0() const {return mixing_gg0;}
    Base_Mixing::Mixing* get_mixing() const {return mixing;}

  private:
  
    // mixing_data
    Base_Mixing::Mixing* mixing = nullptr; ///< Mixing object to mix charge density, kinetic energy density and compensation density
    Base_Mixing::Mixing_Data rho_mdata;    ///< Mixing data for charge density
    Base_Mixing::Mixing_Data tau_mdata;    ///< Mixing data for kinetic energy density
    Base_Mixing::Mixing_Data nhat_mdata;   ///< Mixing data for compensation density
    Base_Mixing::Mixing_Data dmr_mdata;    ///< Mixing data for real space density matrix
    Base_Mixing::Plain_Mixing* mixing_highf = nullptr; ///< The high_frequency part is mixed by plain mixing method.

    //======================================
    // private mixing parameters
    //======================================
    std::string mixing_mode = "broyden"; ///< mixing mode: "plain", "broyden", "pulay"
    double mixing_beta = 0.8;            ///< mixing beta for density
    double mixing_beta_mag = 1.6;        ///< mixing beta for magnetism
    int mixing_ndim = 8;                 ///< mixing ndim for broyden and pulay
    double mixing_gg0 = 0.0;             ///< mixing gg0 for Kerker screen
    bool mixing_tau = false;             ///< whether to use tau mixing
    double mixing_gg0_mag = 0.0;         ///< mixing gg0 for Kerker screen for magnetism
    double mixing_gg0_min = 0.1;         ///< minimum kerker coefficient
    double mixing_angle = 0.0;           ///< mixing angle for nspin=4
    bool mixing_dmr = false;             ///< whether to mixing real space density matrix

    bool new_e_iteration = true;

    ModulePW::PW_Basis* rhopw = nullptr;  ///< smooth grid
    ModulePW::PW_Basis* rhodpw = nullptr; ///< dense grid, same as rhopw for ncpp.

    /**
     * @brief charge mixing for reciprocal space
     * @param chr pointer of Charge object
     */
    void mix_rho_recip_new(Charge* chr);

    /**
     * @brief charge mixing for real space
     * @param chr pointer of Charge object
     */
    void mix_rho_real(Charge* chr);

    /**
     * @brief Kerker screen method for reciprocal space
     * @param rhog charge density in reciprocal space
     */
    void Kerker_screen_recip(std::complex<double>* rhog);
    void Kerker_screen_recip_new(std::complex<double>* rhog);

    /**
     * @brief Kerker screen method for real space
     * @param rho charge density in real space
     */
    void Kerker_screen_real(double* rho);

    /**
     * @brief Inner product of two complex vectors
     * 
     */
    double inner_product_recip(std::complex<double>* rho1, std::complex<double>* rho2);
    double inner_product_recip_new1(std::complex<double>* rho1, std::complex<double>* rho2);
    double inner_product_recip_new2(std::complex<double>* rho1, std::complex<double>* rho2);

    /**
     * @brief Inner product of two double vectors
     *
     */
    double inner_product_real(double* rho1, double* rho2);

    double rhog_dot_product(const std::complex<double>* const* const rhog1,
                            const std::complex<double>* const* const rhog2) const;

    /**
     * @brief divide rho/tau to smooth and high frequency parts
     * @param data_d dense data
     * @param data_s smooth data
     * @param data_hf high frequency data = dense data - smooth data
     *
     */
    void divide_data(std::complex<double>* data_d, std::complex<double>*& data_s, std::complex<double>*& data_hf);
    /**
     * @brief gather smooth and high frequency parts to rho/tau
     * @param data_d dense data
     * @param data_s smooth data
     * @param data_hf high frequency data = dense data - smooth data
     *  
     */
    void combine_data(std::complex<double>* data_d, std::complex<double>*& data_s, std::complex<double>*& data_hf);
    /**
     * @brief clean smooth and high frequency parts
     * @param data_d dense data
     * @param data_s smooth data
     * @param data_hf high frequency data = dense data - smooth data
     *
     */
    void clean_data(std::complex<double>*& data_s, std::complex<double>*& data_hf);
};

#endif
