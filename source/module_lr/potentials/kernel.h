#pragma once
#include "module_basis/module_pw/pw_basis.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_cell/unitcell.h"
// #include <ATen/tensor.h>

namespace LR
{
    class KernelXC
    {
    public:
        KernelXC() {};
        ~KernelXC() {};

        void cal_kernel(const Charge* chg_gs/* ground state*/, const UnitCell* ucell, int& nspin);

        const std::vector<double>& get_kernel(const std::string& name) { return kernel_set_[name]; }
        const std::vector<ModuleBase::Vector3<double>>& get_grad_kernel(const std::string& name) { return grad_kernel_set_[name]; }


    protected:
        // xc kernel for LR-TDDFT
        void f_xc_libxc(const int& nspin, const double& omega, const double& tpiba, const Charge* chg_gs);

        const ModulePW::PW_Basis* rho_basis_ = nullptr;
        std::map<std::string, std::vector<double>> kernel_set_; // [kernel_type][nrxx][nspin]
        std::map<std::string, std::vector<ModuleBase::Vector3<double>>> grad_kernel_set_;// [kernel_type][nrxx][nspin],  intermediate terms for GGA
    };
}

