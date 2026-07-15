#include "esolver.h"

#include "esolver_ks_pw.h"
#include "esolver_sdft_pw.h"
#include "source_base/module_device/device.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __LCAO
#include "esolver_dm2rho.h"
#include "esolver_double_xc.h"
#include "esolver_gets.h"
#include "esolver_ks_lcao.h"
#include "esolver_ks_lcao_tddft.h"
#include "esolver_ks_lcaopw.h"
#include "source_lcao/module_lr/esolver_lrtd_lcao.h"
#include "source_base/module_external/blacs_connector.h"
#endif
#include "esolver_dp.h"
#include "esolver_nep.h"
#include "esolver_lj.h"
#include "esolver_of.h"
#include "esolver_of_tddft.h"

#include <stdexcept>

namespace ModuleESolver
{

#ifdef __LCAO
template <typename TK, typename TR>
class ESolver_KS_LR_LCAO : public ESolver
{
  public:
    ESolver_KS_LR_LCAO(const Input_para& input, UnitCell& ucell)
        : input_(input), ucell_(ucell), ks_(new ESolver_KS_LCAO<TK, TR>())
    {
        this->classname = "ESolver_KS_LR_LCAO";
    }

    void before_all_runners(UnitCell& ucell, const Input_para& input) override
    {
        ks_->before_all_runners(ucell, input);
    }

    void runner(UnitCell& ucell,
                const int istep,
                const ModuleContext::SimulationContext& context) override
    {
        if (!lr_)
        {
            ks_->runner(ucell, istep, context);
            std::cout << " PREPARING FOR EXCITED STATES." << std::endl;
            lr_.reset(new LR::ESolver_LR<TK, TR>(std::move(*ks_), input_, ucell_));
            ks_.reset();
            lr_->before_all_runners(ucell, input_);
        }
        lr_->runner(ucell, istep, context);
        this->conv_esolver = lr_->conv_esolver;
    }

    void after_all_runners(UnitCell& ucell,
                           const ModuleContext::SimulationContext& context) override
    {
        if (lr_)
        {
            lr_->after_all_runners(ucell, context);
        }
    }

    double cal_energy() override
    {
        return lr_ ? lr_->cal_energy() : ks_->cal_energy();
    }

    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override
    {
        if (lr_)
        {
            lr_->cal_force(ucell, force);
        }
        else
        {
            ks_->cal_force(ucell, force);
        }
    }

    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override
    {
        if (lr_)
        {
            lr_->cal_stress(ucell, stress);
        }
        else
        {
            ks_->cal_stress(ucell, stress);
        }
    }

  private:
    Input_para input_;
    UnitCell& ucell_;
    std::unique_ptr<ESolver_KS_LCAO<TK, TR>> ks_;
    std::unique_ptr<LR::ESolver_LR<TK, TR>> lr_;
};
#endif

std::string determine_type()
{
    std::string esolver_type = "none";
    if (PARAM.inp.basis_type == "pw")
    {
        if (PARAM.inp.esolver_type == "sdft")
        {
            esolver_type = "sdft_pw";
        }
        else if (PARAM.inp.esolver_type == "ofdft")
        {
            esolver_type = "ofdft";
        }
        else if (PARAM.inp.esolver_type == "tdofdft")
        {
            esolver_type = "tdofdft";
        }
        else if (PARAM.inp.esolver_type == "ksdft")
        {
            esolver_type = "ksdft_pw";
        }
    }
    else if (PARAM.inp.basis_type == "lcao_in_pw")
    {
#ifdef __LCAO
        if (PARAM.inp.esolver_type == "sdft")
        {
            esolver_type = "sdft_pw";
        }
        else if (PARAM.inp.esolver_type == "ksdft")
        {
            esolver_type = "ksdft_lip";
        }
#else
        ModuleBase::WARNING_QUIT("ESolver", "Calculation involving numerical orbitals must be compiled with __LCAO");
#endif
    }
    else if (PARAM.inp.basis_type == "lcao")
    {
#ifdef __LCAO
        if (PARAM.inp.esolver_type == "tddft")
        {
            esolver_type = "ksdft_lcao_tddft";
        }
        else if (PARAM.inp.esolver_type == "ksdft")
        {
            esolver_type = "ksdft_lcao";
        }
        else if (PARAM.inp.esolver_type == "ks-lr")
        {
            esolver_type = "ksdft_lr_lcao";
        }
        else if (PARAM.inp.esolver_type == "lr")
        {
            esolver_type = "lr_lcao";
        }
#else
        ModuleBase::WARNING_QUIT("ESolver", "Calculation involving numerical orbitals must be compiled with __LCAO");
#endif
    }

    if (PARAM.inp.esolver_type == "lj")
    {
        esolver_type = "lj_pot";
    }
    else if (PARAM.inp.esolver_type == "dp")
    {
        esolver_type = "dp_pot";
    }
    else if (PARAM.inp.esolver_type == "nep")
    {
        esolver_type = "nep_pot";
    }
    else if (esolver_type == "none")
    {
        ModuleBase::WARNING_QUIT("ESolver", "No such esolver_type combined with basis_type");
    }

    GlobalV::ofs_running << "\n #ENERGY SOLVER# " << esolver_type << std::endl;

    auto device_info = PARAM.inp.device;

    for (char& c: device_info)
    {
        if (std::islower(c))
        {
            c = std::toupper(c);
        }
    }
    base_device::information::output_device_info(std::cout, PARAM.inp.device);
    base_device::information::output_device_info(GlobalV::ofs_running, PARAM.inp.device);
    /***auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "hipGetDeviceInfo took " << duration.count() << " seconds" << std::endl;***/
    return esolver_type;
}

// Some API to operate E_Solver
ESolver* init_esolver(const Input_para& inp, UnitCell& ucell)
{
    // determine type of esolver based on INPUT information
    const std::string esolver_type = determine_type();

    // initialize the corresponding Esolver child class
    if (esolver_type == "ksdft_pw")
    {
#if ((defined __CUDA) || (defined __ROCM))
        if (PARAM.inp.device == "gpu")
        {
            if (PARAM.inp.precision == "single")
            {
                return new ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>();
            }
            else
            {
                return new ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>();
            }
        }
#endif
        if (PARAM.inp.precision == "single")
        {
            return new ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>();
        }
        else
        {
            return new ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>();
        }
    }
    else if (esolver_type == "sdft_pw")
    {
#if ((defined __CUDA) || (defined __ROCM))
        if (PARAM.inp.device == "gpu")
        {
            // if (PARAM.inp.precision == "single")
            // {
            //     return new ESolver_SDFT_PW<std::complex<float>, base_device::DEVICE_GPU>();
            // }
            // else
            // {
            return new ESolver_SDFT_PW<std::complex<double>, base_device::DEVICE_GPU>();
            // }
        }
#endif
        // if (PARAM.inp.precision == "single")
        // {
        // 	return new ESolver_SDFT_PW<std::complex<float>, base_device::DEVICE_CPU>();
        // }
        // else
        // {
        return new ESolver_SDFT_PW<std::complex<double>, base_device::DEVICE_CPU>();
        // }
    }
#ifdef __LCAO
    else if (esolver_type == "ksdft_lip")
    {
        if (PARAM.inp.precision == "single")
        {
            return new ESolver_KS_LIP<std::complex<float>>();
        }
        else
        {
            return new ESolver_KS_LIP<std::complex<double>>();
        }
    }
    else if (esolver_type == "ksdft_lcao")
    {
        if (PARAM.inp.calculation == "get_s")
        {
            if (PARAM.globalv.gamma_only_local)
            {
                ModuleBase::WARNING_QUIT("ESolver", "get_s is not implemented for gamma_only");
            }
            else
            {
                return new ESolver_GetS();
            }
        }
        else if (PARAM.inp.deepks_out_base != "none")
        {
            if (PARAM.globalv.gamma_only_local)
            {
                return new ESolver_DoubleXC<double, double>();
            }
            else if (PARAM.inp.nspin < 4)
            {
                return new ESolver_DoubleXC<std::complex<double>, double>();
            }
            else
            {
                return new ESolver_DoubleXC<std::complex<double>, std::complex<double>>();
            }
        }
        else if (PARAM.inp.dm_to_rho)
        {
            if (PARAM.globalv.gamma_only_local)
            {
                ModuleBase::WARNING_QUIT("ESolver", "dm_to_rho is not implemented for gamma_only");
            }
            else if (PARAM.inp.nspin < 4)
            {
                return new ESolver_DM2rho<std::complex<double>, double>();
            }
            else
            {
                return new ESolver_DM2rho<std::complex<double>, std::complex<double>>();
            }
        }
        else
        {
            if (PARAM.globalv.gamma_only_local)
            {
                return new ESolver_KS_LCAO<double, double>();
            }
            else if (PARAM.inp.nspin < 4)
            {
                return new ESolver_KS_LCAO<std::complex<double>, double>();
            }
            else
            {
                return new ESolver_KS_LCAO<std::complex<double>, std::complex<double>>();
            }
        }
    }
    else if (esolver_type == "ksdft_lcao_tddft")
    {
        if (PARAM.inp.nspin < 4)
        {
#if ((defined __CUDA) /* || (defined __ROCM) */)
            if (PARAM.inp.device == "gpu")
            {
                return new ESolver_KS_LCAO_TDDFT<double, base_device::DEVICE_GPU>();
            }
#endif
            return new ESolver_KS_LCAO_TDDFT<double, base_device::DEVICE_CPU>();
        }
        else
        {
#if ((defined __CUDA) /* || (defined __ROCM) */)
            if (PARAM.inp.device == "gpu")
            {
                return new ESolver_KS_LCAO_TDDFT<std::complex<double>, base_device::DEVICE_GPU>();
            }
#endif
            return new ESolver_KS_LCAO_TDDFT<std::complex<double>, base_device::DEVICE_CPU>();
        }
    }
    else if (esolver_type == "lr_lcao")
    {
        // use constructor rather than Init function to initialize reference (instead of pointers) to ucell
        if (PARAM.globalv.gamma_only_local)
        {
            return new LR::ESolver_LR<double, double>(inp, ucell);
        }
        else
        {
            return new LR::ESolver_LR<std::complex<double>, double>(inp, ucell);
        }
    }
    else if (esolver_type == "ksdft_lr_lcao")
    {
        if (PARAM.globalv.gamma_only_local)
        {
            return new ESolver_KS_LR_LCAO<double, double>(inp, ucell);
        }
        else if (PARAM.inp.nspin < 4)
        {
            return new ESolver_KS_LR_LCAO<std::complex<double>, double>(inp, ucell);
        }
        else
        {
            ModuleBase::WARNING_QUIT("init_esolver", "ksdft_lr_lcao supports only nspin=1 or nspin=2");
            return nullptr;
        }
    }
#endif
    else if (esolver_type == "ofdft")
    {
        return new ESolver_OF();
    }
    else if (esolver_type == "tdofdft")
    {
        return new ESolver_OF_TDDFT();
    }
    else if (esolver_type == "lj_pot")
    {
        return new ESolver_LJ();
    }
    else if (esolver_type == "dp_pot")
    {
        return new ESolver_DP(PARAM.mdp.pot_file);
    }
    else if (esolver_type == "nep_pot")
    {
        return new ESolver_NEP(PARAM.mdp.pot_file);
    }
    throw std::invalid_argument("esolver_type = " + std::string(esolver_type) + ". Wrong in " + std::string(__FILE__)
                                + " line " + std::to_string(__LINE__));
}


} // namespace ModuleESolver
