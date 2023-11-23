
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
    Base_Mixing::Mixing* mixing = nullptr;
    Base_Mixing::Mixing_Data rho_mdata;
    Base_Mixing::Mixing_Data tau_mdata;
    Base_Mixing::Mixing_Data nhat_mdata;

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
    void Kerker_screen_real_test(double* rho);

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
     */
    void set_mixing(const std::string& mixing_mode_in,
                    const double& mixing_beta_in,
                    const int& mixing_ndim_in,
                    const double& mixing_gg0_in,
                    const bool& mixing_tau_in,
                    const double& mixing_beta_mag_in = 1.6); // mohan add mixing_gg0_in 2014-09-27

    // /**
    //  * @brief use auto set
    //  *
    //  */
    // void need_auto_set();

    // /**
    //  * @brief auto set mixing gg0 and mixing_beta
    //  *
    //  */
    // void auto_set(const double& bandgap_in, const UnitCell& ucell_);

    /**
     * @brief Get the drho
     *
     */
    double get_drho(Charge* chr, const double nelec);

    // init pwrho and rhodpw
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
    std::string mixing_mode = "broyden";
    double mixing_beta = 0.8;
    double mixing_beta_mag = 1.6;
    int mixing_ndim = 8;
    double mixing_gg0 = 0.0; // mohan add 2014-09-27
    bool mixing_tau = false;

    bool new_e_iteration = true;

    ModulePW::PW_Basis* rhopw = nullptr;  ///< smooth grid
    ModulePW::PW_Basis* rhodpw = nullptr; ///< dense grid, same as rhopw for ncpp.
    bool autoset = false;

  private:
    double rhog_dot_product(const std::complex<double>* const* const rhog1,
                            const std::complex<double>* const* const rhog2) const;

    /**
     * @brief divide rho/tau to smooth and high frequency parts
     *
     */
    void divide_data(std::complex<double>* data_d, std::complex<double>*& data_s, std::complex<double>*& data_hf);
    /**
     * @brief gather smooth and high frequency parts to rho/tau
     *
     */
    void combine_data(std::complex<double>* data_d, std::complex<double>*& data_s, std::complex<double>*& data_hf);
    /**
     * @brief clean smooth and high frequency parts
     *
     */
    void clean_data(std::complex<double>*& data_s, std::complex<double>*& data_hf);
};

#endif
