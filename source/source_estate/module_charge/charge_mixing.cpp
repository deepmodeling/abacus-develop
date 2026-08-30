#include "charge_mixing.h"

#include <functional>
#include <limits>

#include "source_io/module_parameter/parameter.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_base/module_mixing/pulay_mixing.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"

Charge_Mixing::Charge_Mixing()
{
    this->mixing = nullptr;
    this->mixing_uom = nullptr;
    this->mixing_highf = nullptr;
}

Charge_Mixing::~Charge_Mixing()
{
    if(this->mixing != nullptr)
	{
		delete this->mixing;
		this->mixing = nullptr;
	}

	if(this->mixing_highf != nullptr)
	{
		delete this->mixing_highf;
        this->mixing_highf = nullptr;
	}

    if (this->mixing_uom != nullptr)
    {
        delete this->mixing_uom;
        this->mixing_uom = nullptr;
    }
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
                               const bool& mixing_dmr_in,
                               double& omega_in,
                               double& tpiba_in)
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
    this->omega = &omega_in;
    this->tpiba = &tpiba_in;
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
    //GlobalV::ofs_running << "\n\n";
    GlobalV::ofs_running << "\n";
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                 "
           "   |" << std::endl;
    GlobalV::ofs_running << " | Setup charge mixing parameters                                  "
           "   |" << std::endl;
    GlobalV::ofs_running << " |                                                                 "
           "   |" << std::endl;
    GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<" << std::endl;
    GlobalV::ofs_running << "\n";


    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_type", this->mixing_mode);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_beta", this->mixing_beta);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0", this->mixing_gg0);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0_min", PARAM.inp.mixing_gg0_min);

    if (PARAM.inp.nspin==2 || PARAM.inp.nspin==4)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_beta_mag", this->mixing_beta_mag);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0_mag", PARAM.inp.mixing_gg0_mag);
    }
    if (PARAM.inp.mixing_angle > 0)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_angle", PARAM.inp.mixing_angle);
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_ndim", this->mixing_ndim);

    return;
}

void Charge_Mixing::init_mixing()
{
    // this init should be called at the 1-st iteration of each scf loop

    ModuleBase::TITLE("Charge_Mixing", "init_mixing");
    ModuleBase::timer::start("Charge_Mixing", "init_mixing");

    // (re)construct the charge mixer.  The independent UOM mixer is
    // initialized explicitly by the DFT+U call path.

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

    delete this->mixing_uom;
    this->mixing_uom = nullptr;

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
    if ((XC_Functional::get_ked_flag()) && mixing_tau)
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

    ModuleBase::timer::end("Charge_Mixing", "init_mixing");

    return;
}

void Charge_Mixing::init_mixing_uom()
{
    ModuleBase::TITLE("Charge_Mixing", "init_mixing_uom");
    ModuleBase::timer::start("Charge_Mixing", "init_mixing_uom");

    delete this->mixing_uom;
    // The occupation matrix is a small, strongly coupled nonlinear variable.
    // A conservative plain step avoids binding its history to the much
    // larger charge-density mixer.
    this->mixing_uom = new Base_Mixing::Plain_Mixing(0.5 * this->mixing_beta);

    ModuleBase::timer::end("Charge_Mixing", "init_mixing_uom");
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
    if ((XC_Functional::get_ked_flag()) && mixing_tau)
    {
        this->tau_mdata.reset();
    }
    if (this->mixing_uom != nullptr)
    {
        this->mixing_uom->reset();
        this->uom_mdata.reset();
    }
}

bool Charge_Mixing::if_scf_oscillate(const int iteration, const double drho, const int iternum_used, const double threshold)
{
    ModuleBase::TITLE("Charge_Mixing", "if_scf_oscillate");

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
}

void Charge_Mixing::allocate_mixing_uom(int uom_size)
{
    ModuleBase::TITLE("Charge_Mixing", "allocate_mixing_uom");
    ModuleBase::timer::start("Charge_Mixing", "allocate_mixing_uom");
    // For nspin=2, uom_size already includes both spin channels
    // (pot_uterm_pw.size() = pot_index * 2 for nspin=2)
    // So uom_fold should always be 1.  UOM has an independent history.
    if (this->mixing_uom == nullptr)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "UOM mixing object is not initialized.");
    }
    this->mixing_uom->init_mixing_data(this->uom_mdata, uom_size, sizeof(double));
    this->uom_mdata.reset();
    ModuleBase::timer::end("Charge_Mixing", "allocate_mixing_uom");
    return;
}

void Charge_Mixing::mix_uom(std::vector<double>& uom_in, std::vector<double>& uom_save_in)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_uom");
    ModuleBase::timer::start("Charge_Mixing", "mix_uom");

    if (uom_in.empty())
    {
        ModuleBase::timer::end("Charge_Mixing", "mix_uom");
        return;
    }
    if (this->mixing_uom == nullptr || this->uom_mdata.length != uom_in.size()
        || uom_save_in.size() != uom_in.size())
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "UOM mixing data is not initialized consistently.");
    }

    double* uom_value_out = uom_in.data();
    double* uom_value_in = uom_save_in.data();
    const std::size_t length = this->uom_mdata.length;

    // The UOM vector has its own history and its own contractive step.  It
    // must not reuse the charge-density mixer: the two vectors have different
    // lengths and represent different nonlinear variables.
    this->mixing_uom->push_data(this->uom_mdata,
                                uom_value_in,
                                uom_value_out,
                                nullptr,
                                true);

    auto inner_product_uom = [length](double* first, double* second) {
        double value = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : value)
#endif
        for (std::size_t i = 0; i < length; ++i)
        {
            value += first[i] * second[i];
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(value);
#endif
        return value;
    };

    this->mixing_uom->cal_coef(this->uom_mdata, inner_product_uom);
    this->mixing_uom->mix_data(this->uom_mdata, uom_value_out);

    ModuleBase::timer::end("Charge_Mixing", "mix_uom");
#ifdef __MPI
    // Synchronize mixed uom across all ranks to prevent divergence
    // after multiple mixing steps.
    if (uom_in.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing::mix_uom", "UOM buffer is too large for MPI broadcast count");
    }
    Parallel_Common::bcast_double(uom_in.data(), static_cast<int>(uom_in.size()));
#endif
    return;
}
