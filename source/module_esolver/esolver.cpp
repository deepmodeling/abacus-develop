#include "esolver.h"

#include "esolver_ks_pw.h"
#include "esolver_sdft_pw.h"
#include "module_base/module_device/device.h"
#include "module_parameter/parameter.h"
#ifdef __LCAO
#include "esolver_gets.h"
#include "esolver_ks_lcao.h"
#include "esolver_ks_lcao_tddft.h"
#include "esolver_ks_lcaopw.h"
#include "module_lr/esolver_lrtd_lcao.h"
extern "C"
{
#include "module_base/blacs_connector.h"
}
#endif
#include "esolver_dp.h"
#include "esolver_lj.h"
#include "esolver_of.h"
#include "module_parameter/md_parameter.h"
#include <stdexcept>

namespace ModuleESolver
{

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
        else if (PARAM.inp.esolver_type == "ksdft")
        {
            esolver_type = "ksdft_pw";
        }
    }
    else if (PARAM.inp.basis_type == "lcao_in_pw")
    {
#ifdef __LCAO
		if(PARAM.inp.esolver_type == "sdft")
		{
			esolver_type = "sdft_pw";
		}
		else if(PARAM.inp.esolver_type == "ksdft")
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
    else if (esolver_type == "none")
    {
        ModuleBase::WARNING_QUIT("ESolver", "No such esolver_type combined with basis_type");
    }

    GlobalV::ofs_running << " The esolver type has been set to : " << esolver_type << std::endl;

    auto device_info = PARAM.inp.device;

	for (char &c : device_info)
	{
		if (std::islower(c))
		{
			c = std::toupper(c);
		}
	}
	if (GlobalV::MY_RANK == 0)
	{
		std::cout << " RUNNING WITH DEVICE  : " << device_info << " / "
			<< base_device::information::get_device_info(PARAM.inp.device) << std::endl;
	}

    GlobalV::ofs_running << "\n RUNNING WITH DEVICE  : " << device_info << " / "
                         << base_device::information::get_device_info(PARAM.inp.device) << std::endl;

    return esolver_type;
}


//Some API to operate E_Solver
ESolver* init_esolver(const Input_para& inp, UnitCell& ucell)
{
	//determine type of esolver based on INPUT information
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
		if (PARAM.globalv.gamma_only_local)
		{
            if (PARAM.inp.calculation == "get_S")
            {
                return new ESolver_GetS<double, double>();
            }
            else
            {
                return new ESolver_KS_LCAO<double, double>();
            }
        }
        else if (PARAM.inp.nspin < 4)
        {
            if (PARAM.inp.calculation == "get_S")
            {
                return new ESolver_GetS<std::complex<double>, double>();
            }
            else
            {
                return new ESolver_KS_LCAO<std::complex<double>, double>();
            }
        }
        else
        {
            if (PARAM.inp.calculation == "get_S")
            {
                return new ESolver_GetS<std::complex<double>, std::complex<double>>();
            }
            else
            {
                return new ESolver_KS_LCAO<std::complex<double>, std::complex<double>>();
            }
        }
    }
    else if (esolver_type == "ksdft_lcao_tddft")
	{
		return new ESolver_KS_LCAO_TDDFT();
    }
    else if (esolver_type == "lr_lcao")
    {
        // use constructor rather than Init function to initialize reference (instead of pointers) to ucell
        if (PARAM.globalv.gamma_only_local) { return new LR::ESolver_LR<double, double>(inp, ucell); }
        else { return new LR::ESolver_LR<std::complex<double>, double>(inp, ucell); }
    }
    else if (esolver_type == "ksdft_lr_lcao")
    {
        // initialize the 1st ESolver_KS
        ModuleESolver::ESolver* p_esolver = nullptr;
        if (PARAM.globalv.gamma_only_local)
        {
            p_esolver = new ESolver_KS_LCAO<double, double>();
        }
        else if (PARAM.inp.nspin < 4)
        {
            p_esolver = new ESolver_KS_LCAO<std::complex<double>, double>();
        }
        else
        {
            p_esolver = new ESolver_KS_LCAO<std::complex<double>, std::complex<double>>();
        }
        p_esolver->before_all_runners(inp, ucell);
        p_esolver->runner(0, ucell); // scf-only
        // force and stress is not needed currently,
        // they will be supported after the analytical gradient
        // of LR-TDDFT is implemented.
        // after_all_runners() is for output, it is not needed here.
        std::cout << "Setting up the esolver for excited states..." << std::endl;
        // initialize the 2nd ESolver_LR at the temporary pointer
        ModuleESolver::ESolver* p_esolver_lr = nullptr;
        if (PARAM.globalv.gamma_only_local)
            p_esolver_lr = new LR::ESolver_LR<double, double>(
                std::move(*dynamic_cast<ModuleESolver::ESolver_KS_LCAO<double, double>*>(p_esolver)),
                inp,
                ucell);
        else
            p_esolver_lr = new LR::ESolver_LR<std::complex<double>, double>(
                std::move(*dynamic_cast<ModuleESolver::ESolver_KS_LCAO<std::complex<double>, double>*>(p_esolver)),
                inp,
                ucell);
        // clean the 1st ESolver_KS and swap the pointer
        ModuleESolver::clean_esolver(p_esolver, false); // do not call Cblacs_exit, remain it for the 2nd ESolver
        return p_esolver_lr;
    }
#endif
	else if(esolver_type == "ofdft")
	{
		return new ESolver_OF();
	}
	else if (esolver_type == "lj_pot")
	{
		return new ESolver_LJ();
	}
	else if (esolver_type == "dp_pot")
	{
		return new ESolver_DP(PARAM.mdp.pot_file);
	}
	throw std::invalid_argument("esolver_type = "+std::string(esolver_type)+". Wrong in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
}

void clean_esolver(ESolver * &pesolver, const bool lcao_cblacs_exit)
{
// Zhang Xiaoyang modified in 2024/7/6:
// Note: because of the init method of serial lcao hsolver
// it needs no release step for it, or this [delete] will cause Segmentation Fault
// Probably it will be modified later.
#ifdef __MPI
	delete pesolver;
#ifdef __LCAO
    if (lcao_cblacs_exit) {
        Cblacs_exit(1);
}
#endif
#endif
}
}
