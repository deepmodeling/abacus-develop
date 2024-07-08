#ifndef ESOLVER_KS_LIP_H
#define ESOLVER_KS_LIP_H
#include "module_esolver/esolver_ks_pw.h"
#include "module_hsolver/hsolver_lcaopw.h"
namespace ModuleESolver
{

    template <typename T>
    class ESolver_KS_LIP : public ESolver_KS_PW<T, base_device::DEVICE_CPU>
    {
    private:
        using Real = typename GetTypeReal<T>::type;

    public:
        ESolver_KS_LIP();

        ~ESolver_KS_LIP() = default;

        /// All the other interfaces except this one are the same as ESolver_KS_PW.
        virtual void hamilt2density(const int istep, const int iter, const double ethr) override;

    protected:

        virtual void allocate_hsolver() override;
        virtual void deallocate_hsolver() override;

    };
} // namespace ModuleESolver
#endif
