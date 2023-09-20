
#ifndef CHARGE_MIXING_H
#define CHARGE_MIXING_H
#include "charge.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/module_mixing/mixing.h"
#include "module_cell/unitcell.h"
class Charge_Mixing
{
  public:
    Charge_Mixing();
    ~Charge_Mixing();
    Base_Mixing::Mixing* mixing = nullptr;
    Base_Mixing::Mixing_Data rho_mdata;
    Base_Mixing::Mixing_Data tau_mdata;

    void mix_rho(const int& iter, Charge* chr);
    void mix_rho_recip(const int& iter, Charge* chr);
    void mix_rho_real(const int& iter, Charge* chr);

    void Kerker_screen_recip(std::complex<double>* rhog);
    void Kerker_screen_real(double* rho);

    double inner_dot_recip(std::complex<double>* rho1, std::complex<double>* rho2);
    double inner_dot_real(double* rho1, double* rho2);

    //======================================
    // General interfaces, in charge_mixing.cpp
    //======================================
    void set_mixing(const std::string& mixing_mode_in,
                    const double& mixing_beta_in,
                    const int& mixing_ndim_in,
                    const double& mixing_gg0_in,
                    const bool& mixing_tau_in); // mohan add mixing_gg0_in 2014-09-27

    void need_auto_set();
    void auto_set(const double& bandgap_in, const UnitCell& ucell_);

    double get_drho(Charge* chr, const double nelec);

    // init pwrho, sunliang add 2023-05-08
    void set_rhopw(ModulePW::PW_Basis* rhopw_in);

    // extracting parameters
	// normally these parameters will not be used
	// outside charge mixing, but Exx is using them
    // as well as some other places
    const std::string &get_mixing_mode() const {return mixing_mode;}
    double get_mixing_beta() const {return mixing_beta;}
    int get_mixing_ndim() const {return mixing_ndim;}
    double get_mixing_gg0() const {return mixing_gg0;}

  private:
    //======================================
    // General parameters
    //======================================
    std::string mixing_mode = "broyden";
    double mixing_beta = 0.7;
    int mixing_ndim = 8;
    double mixing_gg0 = 0.0; // mohan add 2014-09-27
    bool mixing_tau = false;

    bool new_e_iteration = true;

    ModulePW::PW_Basis* rhopw = nullptr;
    bool autoset = false;

  private:
    double rhog_dot_product(const std::complex<double>* const* const rhog1,
                            const std::complex<double>* const* const rhog2) const;
};

#endif
