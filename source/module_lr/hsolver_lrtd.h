#pragma once
#include "module_hsolver/hsolver.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_psi/psi.h"
namespace LR
{
    template<typename T, typename Device = base_device::DEVICE_CPU>
    class HSolverLR : public hsolver::HSolver<T, Device>
    {
    private:
        using Real = typename GetTypeReal<T>::type;
        const int nks = 0;
        const int npairs = 0;
        const bool out_wfc_lr = false;
    public:
        HSolverLR(const int nks_in, const int npairs_in, const bool out_wfc_lr_in = false)
            :nks(nks_in), npairs(npairs_in), out_wfc_lr(out_wfc_lr_in) {};
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