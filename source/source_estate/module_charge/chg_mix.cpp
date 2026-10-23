#include "chg_mix.h"
#include "chg_drho.h"
#include "chg_precond.h"
#include "chg_rho_detail.h"
#include "chg_tau.h"
#include "chg_uspp.h"

#include <functional>
#include <memory>

#include "source_base/module_mixing/broyden_mixing.h"
#include "source_base/module_mixing/pulay_mixing.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

Charge_Mixing::Charge_Mixing()
{
    // unique_ptr members default-construct to nullptr
}

Charge_Mixing::~Charge_Mixing()
{
    // unique_ptr members (mixing, mixing_highf) are released automatically
}

void Charge_Mixing::set_mixing(const MixingConfig& cfg,
                               ModulePW::PW_Basis* rhopw_in,
                               ModulePW::PW_Basis* rhodpw_in,
                               double& omega_in,
                               double& tpiba_in)
{
    // store the aggregated config; init_mixing/mix_rho and the stateless
    // Kerker kernels all read nspin, scf_thr_type, double_grid, mixing_gg0,
    // mixing_gg0_mag, mixing_gg0_min, mixing_angle, mixing_dmr from cfg_
    // instead of PARAM/GlobalV. cfg_ is treated as an immutable INPUT
    // snapshot; runtime overrides (e.g. close_kerker_gg0) live as flags on
    // Charge_Mixing itself, never by mutating cfg_.
    this->cfg_ = cfg;
    // store the smooth and dense grids
    this->rhopw = rhopw_in;
    this->rhodpw = rhodpw_in;
    // omega and tpiba are pointers to external runtime state (cell volume
    // and lattice constant) that changes across SCF iterations; they are
    // not INPUT parameters and therefore stay out of MixingConfig.
    this->omega = &omega_in;
    this->tpiba = &tpiba_in;
    // check the paramters
    if (this->cfg_.mixing_beta > 1.0 || this->cfg_.mixing_beta < 0.0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "You'd better set mixing_beta to [0.0, 1.0]!");
    }
    if (this->cfg_.nspin >= 2 && this->cfg_.mixing_beta_mag < 0.0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "You'd better set mixing_beta_mag >= 0.0!");
    }

    if (!(this->cfg_.mixing_mode == "plain" || this->cfg_.mixing_mode == "broyden" || this->cfg_.mixing_mode == "pulay"))
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


    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_type", this->cfg_.mixing_mode);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_beta", this->cfg_.mixing_beta);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0", cfg_.mixing_gg0);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0_min", cfg_.mixing_gg0_min);

    if (this->cfg_.nspin==2 || this->cfg_.nspin==4)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_beta_mag", this->cfg_.mixing_beta_mag);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_gg0_mag", cfg_.mixing_gg0_mag);
    }
    if (cfg_.mixing_angle > 0)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_angle", cfg_.mixing_angle);
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "mixing_ndim", this->cfg_.mixing_ndim);

    return;
}

void Charge_Mixing::init_mixing()
{
    // this init should be called at the 1-st iteration of each scf loop

    ModuleBase::TITLE("Charge_Mixing", "init_mixing");
    ModuleBase::timer::start("Charge_Mixing", "init_mixing");

    /// Fail fast when set_mixing was skipped: the grid sizes below would
    /// otherwise dereference a null pointer.
    if (this->rhopw == nullptr)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing",
                                 "set_mixing must be called before init_mixing");
    }

    // (re)construct mixing object
    if (this->cfg_.mixing_mode == "broyden")
    {
        this->mixing = std::unique_ptr<Base_Mixing::Broyden_Mixing>(
            new Base_Mixing::Broyden_Mixing(this->cfg_.mixing_ndim, this->cfg_.mixing_beta));
    }
    else if (this->cfg_.mixing_mode == "plain")
    {
        this->mixing = std::unique_ptr<Base_Mixing::Plain_Mixing>(
            new Base_Mixing::Plain_Mixing(this->cfg_.mixing_beta));
    }
    else if (this->cfg_.mixing_mode == "pulay")
    {
        this->mixing = std::unique_ptr<Base_Mixing::Pulay_Mixing>(
            new Base_Mixing::Pulay_Mixing(this->cfg_.mixing_ndim, this->cfg_.mixing_beta));
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing", "This Mixing mode is not implemended yet,coming soon.");
    }

    if ( this->cfg_.double_grid)
    {
        // ONLY smooth part of charge density is mixed by specific mixing method
        // The high_frequency part is mixed by plain mixing method.
        this->mixing_highf = std::unique_ptr<Base_Mixing::Plain_Mixing>(
            new Base_Mixing::Plain_Mixing(this->cfg_.mixing_beta));
    }

    // allocate memory for mixing data, if exists, free it first and then allocate new memory
    // initailize rho_mdata
    if (this->cfg_.scf_thr_type == 1)
    {
        if (this->cfg_.nspin == 4 && this->cfg_.mixing_angle > 0 )
        {
            this->mixing->init_mixing_data(this->rho_mdata,
                                        this->rhopw->npw * 2,
                                        sizeof(std::complex<double>));
        }
        else
        {
            this->mixing->init_mixing_data(this->rho_mdata,
                                        this->rhopw->npw * this->cfg_.nspin,
                                        sizeof(std::complex<double>));
        }
    }
    else
    {
        if (this->cfg_.nspin == 4 && this->cfg_.mixing_angle > 0 )
        {
            this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * 2, sizeof(double));
        }
        else
        {
            this->mixing->init_mixing_data(this->rho_mdata, this->rhopw->nrxx * this->cfg_.nspin, sizeof(double));
        }
    }

    // initailize tau_mdata
    if (cfg_.mixing_tau)
    {
        if (this->cfg_.scf_thr_type == 1)
        {
            this->mixing->init_mixing_data(this->tau_mdata,
                                           this->rhopw->npw * this->cfg_.nspin,
                                           sizeof(std::complex<double>));
        }
        else
        {
            this->mixing->init_mixing_data(this->tau_mdata, this->rhopw->nrxx * this->cfg_.nspin, sizeof(double));
        }
    }

    ModuleBase::timer::end("Charge_Mixing", "init_mixing");

    return;
}

void Charge_Mixing::mix_reset()
{
    this->mixing->reset();
    this->rho_mdata.reset();
    // initailize tau_mdata
    if (cfg_.mixing_tau)
    {
        this->tau_mdata.reset();
    }
}

bool Charge_Mixing::if_scf_oscillate(const int iteration, const double drho,
                                     const int iternum_used, const double threshold)
{
    ModuleBase::TITLE("Charge_Mixing", "if_scf_oscillate");

    /// Fail fast when set_mixing was skipped: cfg_.scf_nmax is otherwise an
    /// indeterminate value, which would resize the history to zero and then
    /// write out of bounds (heap corruption instead of a clear error).
    if (this->cfg_.scf_nmax <= 0)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing",
                                 "set_mixing must be called before if_scf_oscillate");
    }
    if (iteration < 1 || iteration > this->cfg_.scf_nmax)
    {
        ModuleBase::WARNING_QUIT("Charge_Mixing",
                                 "iteration must be within [1, scf_nmax]");
    }

    if(this->_drho_history.size() == 0)
    {
        this->_drho_history.resize(this->cfg_.scf_nmax);
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
    // (uterm_mat.size() = pot_index * 2 for nspin=2)
    // So uom_fold should always be 1
    this->mixing->init_mixing_data(this->uom_mdata, uom_size, sizeof(double));
    this->uom_mdata.reset();
    ModuleBase::timer::end("Charge_Mixing", "allocate_mixing_uom");
    return;
}

void Charge_Mixing::mix_uom(std::vector<double>& uom_in, std::vector<double>& uom_save_in)
{
    ModuleBase::TITLE("Charge_Mixing", "mix_uom");
    ModuleBase::timer::start("Charge_Mixing", "mix_uom");
    double* uom_value_out = uom_in.data();
    double* uom_value_in = uom_save_in.data();
    // For all nspin cases, uom_array layout is already fully sized
    // and mixing operates on the entire array
    this->mixing->push_data(this->uom_mdata, uom_value_in, uom_value_out, nullptr, false);
    this->mixing->mix_data(this->uom_mdata, uom_value_out);
    ModuleBase::timer::end("Charge_Mixing", "mix_uom");
#ifdef __MPI
    // Synchronize mixed uom across all ranks to prevent divergence
    // after multiple Pulay steps (same pattern as mix_dmr)
    Parallel_Common::bcast_double(uom_in.data(), uom_in.size());
#endif
    return;
}
