#include "charge_mixing.h"

#include "module_parameter/parameter.h"
#include "module_base/module_mixing/broyden_mixing.h"
#include "module_base/module_mixing/pulay_mixing.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

Charge_Mixing::Charge_Mixing()
{
}

Charge_Mixing::~Charge_Mixing()
{
    delete this->mixing;
    delete this->mixing_highf;
}

void Charge_Mixing::set_mixing(const std::string& mixing_mode_in,
                               const double& mixing_beta_in,
                               const int& mixing_ndim_in,
                               const double& mixing_gg0_in,
                               const bool& mixing_tau_in,
                               const double& mixing_beta_mag_in,
                               const double& mixing_gg0_mag_in,
                               const double& mixing_gg0_min_in,
                               const double& mixing_angle_in,
                               const bool& mixing_dmr_in)
{
    // get private mixing parameters
    this->mixing_mode = mixing_mode_in;
    this->mixing_beta = mixing_beta_in;
    this->mixing_beta_mag = mixing_beta_mag_in;
    this->mixing_ndim = mixing_ndim_in;
    this->mixing_gg0 = mixing_gg0_in;
    this->mixing_tau = mixing_tau_in;
    this->mixing_gg0_mag = mixing_gg0_mag_in;
    this->mixing_gg0_min = mixing_gg0_min_in;
    this->mixing_angle = mixing_angle_in;
    this->mixing_dmr = mixing_dmr_in;

    // check the paramters
    if (this->mixing_beta > 1.0 || this->mixing_beta < 0.0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "You'd better set mixing_beta to [0.0, 1.0]!");
    }
    if (PARAM.inp.nspin >= 2 && this->mixing_beta_mag < 0.0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "You'd better set mixing_beta_mag >= 0.0!");
    }

    if (!(this->mixing_mode == "plain" || this->mixing_mode == "broyden" || this->mixing_mode == "pulay"))
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing mode is not implemended yet,coming soon.");
    }

    // print into running.log
    GlobalV::ofs_running<<"\n----------- Double Check Mixing Parameters Begin ------------"<<std::endl;
    GlobalV::ofs_running<<"mixing_type: "<< this->mixing_mode <<std::endl;
    GlobalV::ofs_running<<"mixing_beta: "<< this->mixing_beta <<std::endl;
    GlobalV::ofs_running<<"mixing_gg0: "<< this->mixing_gg0 <<std::endl;
    GlobalV::ofs_running<<"mixing_gg0_min: "<< PARAM.inp.mixing_gg0_min <<std::endl;
    if (PARAM.inp.nspin==2 || PARAM.inp.nspin==4)
    {
        GlobalV::ofs_running<<"mixing_beta_mag: "<< this->mixing_beta_mag <<std::endl;
        GlobalV::ofs_running<<"mixing_gg0_mag: "<< PARAM.inp.mixing_gg0_mag <<std::endl;
    }
    if (PARAM.inp.mixing_angle > 0)
    {
        GlobalV::ofs_running<<"mixing_angle: "<< PARAM.inp.mixing_angle <<std::endl;
    }
    GlobalV::ofs_running<<"mixing_ndim: "<< this->mixing_ndim <<std::endl;
    GlobalV::ofs_running<<"----------- Double Check Mixing Parameters End ------------"<<std::endl;

    return;
}

void Charge_Mixing::init_mixing()
{
    // this init should be called at the 1-st iteration of each scf loop

    ModuleBase::TITLE("Charge_Mixing", "init_mixing");
    ModuleBase::timer::tick("Charge_Mixing", "init_mixing");

    // (re)construct mixing object
    if (this->mixing_mode == "broyden")
    {
        delete this->mixing;
        this->mixing = new Base_Mixing::Broyden_Mixing(this->mixing_ndim, this->mixing_beta);
    }
    else if (this->mixing_mode == "plain")
    {
        delete this->mixing;
        this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
    }
    else if (this->mixing_mode == "pulay")
    {
        delete this->mixing;
        this->mixing = new Base_Mixing::Pulay_Mixing(this->mixing_ndim, this->mixing_beta);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing mode is not implemended yet,coming soon.");
    }

    if ( PARAM.globalv.double_grid)
    {
        // ONLY smooth part of charge density is mixed by specific mixing method
        // The high_frequency part is mixed by plain mixing method.
        delete this->mixing_highf;
        this->mixing_highf = new Base_Mixing::Plain_Mixing(this->mixing_beta);
    }

    // allocate memory for mixing data, if exists, free it first and then allocate new memory
    // initailize rho_mdata
    if (PARAM.inp.scf_thr_type == 1)
    {  
        if (PARAM.inp.nspin == 4 && PARAM.inp.mixing_angle > 0 )
        {
            this->mixing->init_mixing_data(this->rho_mdata,
                                        this->rhopw->npw * 2,
                                        sizeof(std::complex<double>));
        }
        else
        {
            this->mixing->init_mixing_data(this->rho_mdata,
                                        this->rhopw->npw * PARAM.inp.nspin,
                                        sizeof(std::complex<double>));
        }
    }
    else
    {
        if (PARAM.inp.nspin == 4 && PARAM.inp.mixing_angle > 0 )
        {
            this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * 2, sizeof(double));
        }
        else
        {
            this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * PARAM.inp.nspin, sizeof(double));
        }
    }
    
    // initailize tau_mdata
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        if (PARAM.inp.scf_thr_type == 1)
        {
            this->mixing->init_mixing_data(this->tau_mdata,
                                           this->rhopw->npw * PARAM.inp.nspin,
                                           sizeof(std::complex<double>));
        }
        else
        {
            this->mixing->init_mixing_data(this->tau_mdata, this->rhopw->nrxx * PARAM.inp.nspin, sizeof(double));
        }
    }

    // initailize nhat_mdata
#ifdef USE_PAW
    if(PARAM.inp.use_paw) { this->mixing->init_mixing_data(this->nhat_mdata, this->rhopw->nrxx * PARAM.inp.nspin, sizeof(double));
}
#endif

    ModuleBase::timer::tick("Charge_Mixing", "init_mixing");

    return;
}

void Charge_Mixing::set_rhopw(ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis* rhodpw_in)
{
    this->rhopw = rhopw_in;
    this->rhodpw = rhodpw_in;
}

void Charge_Mixing::mix_reset()
{
    this->mixing->reset();
    this->rho_mdata.reset();
    // initailize tau_mdata
    if ((XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) && mixing_tau)
    {
        this->tau_mdata.reset();
    }
    // reset for paw
#ifdef USE_PAW
    this->nhat_mdata.reset();
#endif
}

bool Charge_Mixing::if_scf_oscillate(const int iteration, const double drho, const int iternum_used, const double threshold)
{
    ModuleBase::TITLE("Charge_Mixing", "if_scf_oscillate");
    ModuleBase::timer::tick("Charge_Mixing", "if_scf_oscillate");

    if(this->_drho_history.size() == 0)
    {
        this->_drho_history.resize(PARAM.inp.scf_nmax);
    }

    // add drho into history
    this->_drho_history[iteration - 1] = drho;

    if(threshold >= 0) // close the function
    {
        return false;
    }

    // check if the history is long enough
    if(iteration < iternum_used + this->mixing_restart_last)
    {
        return false;
    }

    // calculate the slope of the last iternum_used iterations' drho
    double slope = 0.0;

    // Least Squares Method
    // this part is too short, so I do not design it as a free function in principle
    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
    for (int i = iteration - iternum_used; i < iteration; i++)
    {
        sumX += i;
        sumY += std::log10(this->_drho_history[i]);
        sumXY += i * std::log10(this->_drho_history[i]);
        sumXX += i * i;
    }
    double numerator = iternum_used * sumXY - sumX * sumY;
    double denominator = iternum_used * sumXX - sumX * sumX;
    if (denominator == 0) {
        return false;
    }
    slope =  numerator / denominator;

    // if the slope is less than the threshold, return true
    if(slope > threshold)
    {
        return true;
    }

    return false;

    ModuleBase::timer::tick("Charge_Mixing", "if_scf_oscillate");
  
}