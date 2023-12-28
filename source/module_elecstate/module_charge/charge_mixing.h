
#ifndef CHARGE_MIXING_H
#define CHARGE_MIXING_H
#include "charge.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/module_mixing/mixing.h"
#include "module_base/module_mixing/plain_mixing.h"
#include "module_cell/unitcell.h"
class Charge_Mixing
{
  public:
    Charge_Mixing();
    ~Charge_Mixing();
    Base_Mixing::Mixing* mixing = nullptr; ///< Mixing object to mix charge density, kinetic energy density and compensation density
    Base_Mixing::Mixing_Data rho_mdata;    ///< Mixing data for charge density
    Base_Mixing::Mixing_Data tau_mdata;    ///< Mixing data for kinetic energy density
    Base_Mixing::Mixing_Data nhat_mdata;   ///< Mixing data for compensation density

    Base_Mixing::Plain_Mixing* mixing_highf = nullptr; ///< The high_frequency part is mixed by plain mixing method.

    /**
     * @brief reset mixing
     *
     */
    void mix_reset();

    /**
     * @brief charge mixing
     *
     */
    void mix_rho(Charge* chr);

    /**
     * @brief charge mixing for reciprocal space
     *
     */
    void mix_rho_recip(Charge* chr);
    void mix_rho_recip_new(Charge* chr);

    /**
     * @brief charge mixing for real space
     *
     */
    void mix_rho_real(Charge* chr);

    /**
     * @brief Kerker screen method for reciprocal space
     *
     */
    void Kerker_screen_recip(std::complex<double>* rhog);
    void Kerker_screen_recip_new(std::complex<double>* rhog);

    /**
     * @brief Kerker screen method for real space
     *
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

    /**
     * @brief Set the mixing object
     *
     * @param mixing_mode_in mixing mode: "plain", "broyden", "pulay"
     * @param mixing_beta_in mixing beta
     * @param mixing_ndim_in mixing ndim
     * @param mixing_gg0_in mixing gg0 for Kerker screen
     * @param mixing_tau_in whether to use tau mixing
     * @param mixing_beta_mag_in mixing beta for magnetism
     * @param mixing_gg0_mag_in mixing gg0 for magnetism
     * @param mixing_angle_in mixing angle for non-colinear angle mixing
     */
    void set_mixing(const std::string& mixing_mode_in,
                    const double& mixing_beta_in,
                    const int& mixing_ndim_in,
                    const double& mixing_gg0_in,
                    const bool& mixing_tau_in,
                    const double& mixing_beta_mag_in,
                    const double& mixing_gg0_mag_in,
                    const double& mixing_angle_in);

    /**
     * @brief Get the drho
     *
     */
    double get_drho(Charge* chr, const double nelec);

    // init pwrho and rhodpw
    
    /**
     * @brief Set the smooth and dense grids
     * 
     * @param rhopw_in smooth grid
     * @param rhodpw_in dense grid when double grid is used, otherwise same as rhopw
     */
    void set_rhopw(ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis* rhodpw_in);

    // extracting parameters
    // normally these parameters will not be used
    // outside charge mixing, but Exx is using them
    // as well as some other places
    const std::string& get_mixing_mode() const
    {
        return mixing_mode;
    }
    double get_mixing_beta() const
    {
        return mixing_beta;
    }
    int get_mixing_ndim() const
    {
        return mixing_ndim;
    }
    double get_mixing_gg0() const
    {
        return mixing_gg0;
    }

  private:
    //======================================
    // General parameters
    //======================================
    std::string mixing_mode = "broyden"; ///< mixing mode: "plain", "broyden", "pulay"
    double mixing_beta = 0.8;            ///< mixing beta for density
    double mixing_beta_mag = 1.6;        ///< mixing beta for magnetism
    int mixing_ndim = 8;                 ///< mixing ndim for broyden and pulay
    double mixing_gg0 = 0.0;             ///< Kerker screen gg0 for density
    double mixing_gg0_mag = 0.0;         ///< Kerker screen gg0 for magnetism
    double mixing_angle = -10.0;         ///< mixing beta for non-colinear angle mixing
    bool mixing_tau = false;             ///< whether to use tau mixing
    int _nspin = 1;                      ///< number of spin      
    int _mixing_rho_type_num = 1;        ///< how many types of rho stored in mixing_data. only charge density: 1; charge density and kinetic energy density: 2 
    int _mixing_rho_unit_num = 1;        ///< number of data_unit for each rho_block in mixing_data, depending on nspin and whether mixing_tau is used

    bool new_e_iteration = true;

    ModulePW::PW_Basis* rhopw = nullptr;  ///< smooth grid
    ModulePW::PW_Basis* rhodpw = nullptr; ///< dense grid, same as rhopw for ncpp.

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
