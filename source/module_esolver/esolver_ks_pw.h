#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/velocity_pw.h"
#include "module_psi/psi_initializer.h"
#include <module_base/macros.h>

// #include "Basis_PW.h"
// #include "Estate_PW.h"
// #include "Hamilton_PW.h"
// #include "H2E_pw.h"


namespace ModuleESolver
{

    template<typename T, typename Device = psi::DEVICE_CPU>
    class ESolver_KS_PW : public ESolver_KS<T, Device>
    {
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        ESolver_KS_PW();
        ~ESolver_KS_PW();
        void Init(Input& inp, UnitCell& cell) override;
        void init_after_vc(Input& inp, UnitCell& cell) override;
        double cal_Energy() override;
        void cal_Force(ModuleBase::matrix& force) override;
        void cal_Stress(ModuleBase::matrix& stress) override;
        virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
        virtual void hamilt2estates(const double ethr) override;
        virtual void nscf() override;
        void postprocess() override;

        /**
         * @brief calculate Onsager coefficients Lmn(\omega) and conductivities with Kubo-Greenwood formula
         * 
         * @param fwhmin FWHM for delta function
         * @param wcut cutoff \omega for Lmn(\omega)
         * @param dw_in \omega step
         * @param dt_in time step
         * @param wg wg(ik,ib) occupation for the ib-th band in the ik-th kpoint
         */
        void KG(const double fwhmin,
                const double wcut,
                const double dw_in,
                const double dt_in,
                ModuleBase::matrix& wg);
        
        /**
         * @brief calculate the response function Cmn(t) for currents
         * 
         * @param ik k point
         * @param nt number of steps of time
         * @param dt time step
         * @param decut ignore dE which is larger than decut
         * @param wg wg(ik,ib) occupation for the ib-th band in the ik-th kpoint
         * @param velop velocity operator
         * @param ct11 C11(t)
         * @param ct12 C12(t)
         * @param ct22 C22(t)
         */
        void jjcorr_ks(const int ik,
                       const int nt,
                       const double dt,
                       const double decut,
                       ModuleBase::matrix& wg,
                       hamilt::Velocity& velop,
                       double* ct11,
                       double* ct12,
                       double* ct22);

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
        void calcondw(const int nt,const double dt, const double fwhmin, const double wcut, const double dw_in, double *ct11, double *ct12, double *ct22);
        /// @brief allocate psi_init the new psi_initializer
        void allocate_psi_init();
        /// @brief initialize psi
        void initialize_psi();
    protected:
        psi::Psi<std::complex<double>, Device>* psi = nullptr;   //hide the psi in ESolver_KS for tmp use
    private:
        psi_initializer* psi_init = nullptr;
        Device * ctx = {};
        psi::AbacusDevice_t device = {};
        psi::Psi<T, Device>* kspw_psi = nullptr;
        psi::Psi<std::complex<double>, Device>* __kspw_psi = nullptr;
        using castmem_2d_d2h_op = psi::memory::cast_memory_op<std::complex<double>, T, psi::DEVICE_CPU, Device>;
    };
}  // namespace ModuleESolver
#endif
