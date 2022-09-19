#ifndef ESOLVER_KS_LCAO_H
#define ESOLVER_KS_LCAO_H
#include "./esolver_ks.h"

#include "src_lcao/record_adj.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/LCAO_hamilt.h"
#include "module_orbital/ORB_control.h"

namespace ModuleESolver
{

    class ESolver_KS_LCAO : public ESolver_KS
    {
    public:
        ESolver_KS_LCAO();
        ~ESolver_KS_LCAO();

        void Init(Input& inp, UnitCell_pseudo& cell) override;

        void cal_Energy(energy& en) override;
        void cal_Force(ModuleBase::matrix& force) override;
        void cal_Stress(ModuleBase::matrix& stress) override;
        void postprocess() override;
        void nscf() override;
        void get_S();

    protected:
        virtual void beforescf(const int istep) override;
        virtual void eachiterinit(const int istep, const int iter) override;
        virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
        virtual void updatepot(const int istep, const int iter) override;
        virtual void eachiterfinish(const int iter) override;
        virtual void afterscf(const int istep) override;
        virtual bool do_after_converge(int& iter) override;
        int two_level_step = 0;

        virtual void othercalculation(const int istep)override;
        ORB_control orb_con;    //Basis_LCAO
        Record_adj RA;
        Local_Orbital_wfc LOWF;
        Local_Orbital_Charge LOC;
        LCAO_Hamilt UHM;
        LCAO_Matrix LM;


        // Temporarily store the stress to unify the interface with PW,
        // because it's hard to seperate force and stress calculation in LCAO.
        // The copy costs memory and time ! 
        // Are there any better way to solve this problem ?
        ModuleBase::matrix scs;
        bool have_force = false;

        void Init_Basis_lcao(ORB_control& orb_con, Input& inp, UnitCell_pseudo& ucell);

        // output subfuncs: implemented in src_io/write_HS_R.cpp
        void output_HS_R(
            const int &istep,
            const std::string& SR_filename = "data-SR-sparse_SPIN0.csr",
            const std::string& HR_filename_up = "data-HR-sparse_SPIN0.csr",
            const std::string HR_filename_down = "data-HR-sparse_SPIN1.csr",
            const bool& binary = false,
            const double& sparse_threshold = 1e-10
        ); //LiuXh add 2019-07-15, modify in 2021-12-3
        void output_SR(const std::string& SR_filename, const bool& binary = false, const double& sparse_threshold = 1e-10);

        //--------------common for all calculation, not only scf-------------
        // set matrix and grid integral
        void set_matrix_grid(Record_adj& ra);
        void beforesolver(const int istep);
        //----------------------------------------------------------------------
    };



}
#endif
