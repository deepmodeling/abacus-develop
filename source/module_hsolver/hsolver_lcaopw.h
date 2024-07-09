#ifndef HSOLVERLIP_H
#define HSOLVERLIP_H

#include "hsolver.h"
#include "module_base/macros.h"
#include "module_base/module_device/types.h"
namespace hsolver {

    // LCAO-in-PW does not support GPU now.
    template <typename T>
    class HSolverLIP : public HSolver<T, base_device::DEVICE_CPU>
    {
    private:
        // Note GetTypeReal<T>::type will 
        // return T if T is real type(float, double), 
        // otherwise return the real type of T(complex<float>, complex<double>)
        using Real = typename GetTypeReal<T>::type;

    public:

        HSolverLIP(ModulePW::PW_Basis_K* wfc_basis_in);

        /// @brief solve function for lcao_in_pw
        /// @param pHamilt interface to hamilt
        /// @param psi reference to psi
        /// @param pes interface to elecstate
        /// @param transform transformation matrix between lcao and pw
        /// @param skip_charge 
        void solve(hamilt::Hamilt<T>* pHamilt,
            psi::Psi<T>& psi,
            elecstate::ElecState* pes,
            psi::Psi<T>& transform,
            const bool skip_charge) override;

    protected:

        ModulePW::PW_Basis_K* wfc_basis = nullptr;

        std::vector<Real> eigenvalues;
        using castmem_2d_2h_op
            = base_device::memory::cast_memory_op<double, Real, base_device::DEVICE_CPU, base_device::DEVICE_CPU>;
    };

} // namespace hsolver

#endif