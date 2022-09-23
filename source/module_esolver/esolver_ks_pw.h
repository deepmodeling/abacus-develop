#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
// #include "Basis_PW.h"
// #include "Estate_PW.h"
// #include "Hamilton_PW.h"
// #include "H2E_pw.h"


namespace ModuleESolver
{

    class ESolver_KS_PW : public ESolver_KS
    {
    public:
        ESolver_KS_PW();
        ~ESolver_KS_PW();
        void Init(Input& inp, UnitCell_pseudo& cell) override;
        void cal_Energy(energy& en) override;
        void cal_Force(ModuleBase::matrix& force) override;
        void cal_Stress(ModuleBase::matrix& stress) override;
        virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
        virtual void hamilt2estates(const double ethr) override;
        virtual void nscf() override;
        void postprocess() override;
        //calculate conductivities with Kubo-Greenwood formula
        void KG(const int nche_KG, const double fwhmin, const double wcut, 
             const double dw_in, const int times);

    protected:
        virtual void beforescf(const int istep) override;
        virtual void eachiterinit(const int istep, const int iter) override;
        virtual void updatepot(const int istep, const int iter) override;
        virtual void eachiterfinish(const int iter) override;
        virtual void afterscf(const int istep) override;
        virtual void othercalculation(const int istep)override;

        //temporary, this will be removed in the future;
        //Init Global class
        void Init_GlobalC(Input& inp, UnitCell_pseudo& cell);
        //calculate conductivities from j-j correlation function
        void calcondw(const int nt,const double dt,const double fwhmin,const double wcut,const double dw_in,double*ct11,double*ct12,double *ct22);


    private:
        // It copies the function in Threshold_Elec class.
        // After all ESolver, HSolver are constructed, Class Electrons and Threshold_Elec should be deleted.
        void print_eigenvalue(std::ofstream& ofs);

    };
}
#endif
