#pragma once
#include "module_hsolver/hsolver.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_psi/psi.h"
namespace LR
{
    template<typename T, typename Device = base_device::DEVICE_CPU>
    class HSolverLR : public hsolver::HSolver<T, Device>
    {
        using Real = typename GetTypeReal<T>::type;
        const int& nk;
        const int& npairs;
        const int& ispin_solve;
        const bool out_wfc_lr = false;
    public:
        HSolverLR(const int& nk_in, const int& npairs_in, const int& ispin_solve_in = 0, const bool& out_wfc_lr_in = false)
            :nk(nk_in), npairs(npairs_in), out_wfc_lr(out_wfc_lr_in), ispin_solve(ispin_solve_in) {};
        virtual Real set_diagethr(const int istep, const int iter, const Real ethr) override
        {
            this->diag_ethr = ethr;
            return ethr;
        }
        virtual void solve(hamilt::Hamilt<T, Device>* pHamilt,
            psi::Psi<T, Device>& psi,
            elecstate::ElecState* pes,
            const std::string method_in,
            const bool skip_charge = false) override;
    };
};