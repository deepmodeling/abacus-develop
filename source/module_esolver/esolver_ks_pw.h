#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
// #include "Basis_PW.h"
// #include "Estate_PW.h"
// #include "Hamilton_PW.h"
// #include "H2E_pw.h"


namespace ModuleESolver
{

    template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
    class ESolver_KS_PW : public ESolver_KS<FPTYPE, Device>
    {
    public:
        ESolver_KS_PW();
        ~ESolver_KS_PW();
        void Init(Input& inp, UnitCell& cell) override;
        void cal_Energy(FPTYPE& etot) override;
        void cal_Force(ModuleBase::matrix& force) override;
        void cal_Stress(ModuleBase::matrix& stress) override;
        virtual void hamilt2density(const int istep, const int iter, const FPTYPE ethr) override;
        virtual void hamilt2estates(const FPTYPE ethr) override;
        virtual void nscf() override;
        void postprocess() override;
        //calculate conductivities with Kubo-Greenwood formula
        void KG(const int nche_KG, const FPTYPE fwhmin, const FPTYPE wcut, 
             const FPTYPE dw_in, const int times, ModuleBase::matrix& wg);

    protected:
        virtual void beforescf(const int istep) override;
        virtual void eachiterinit(const int istep, const int iter) override;
        virtual void updatepot(const int istep, const int iter) override;
        virtual void eachiterfinish(const int iter) override;
        virtual void afterscf(const int istep) override;
        virtual void othercalculation(const int istep)override;

        //temporary, this will be removed in the future;
        //Init Global class
        void Init_GlobalC(Input& inp, UnitCell& cell);
        //calculate conductivities from j-j correlation function
        void calcondw(const int nt,const FPTYPE dt, const FPTYPE fwhmin, const FPTYPE wcut, const FPTYPE dw_in, FPTYPE *ct11, FPTYPE *ct12, FPTYPE *ct22);


    private:
        // It copies the function in Threshold_Elec class.
        // After all ESolver, HSolver are constructed, Class Electrons and Threshold_Elec should be deleted.
        void print_eigenvalue(std::ofstream& ofs);
        
        Device * ctx = {};
        psi::AbacusDevice_t device = {};
        psi::Psi<std::complex<FPTYPE>, Device>* kspw_psi = nullptr;
        using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU, Device>;
    };
}  // namespace ModuleESolver
#endif
