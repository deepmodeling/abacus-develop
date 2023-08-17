#ifndef ESOLVER_KS_LCAO_H
#define ESOLVER_KS_LCAO_H
#include "esolver_ks.h"

#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_basis/module_ao/ORB_control.h"
#ifdef __EXX
#include "module_ri/Mix_DMk_2D.h"
#include "module_ri/Exx_LRI_interface.h"
#endif
#include "module_io/output_dm.h"
#include "module_io/output_dm1.h"
#include "module_io/output_mat_sparse.h"
#include "module_basis/module_nao/two_center_bundle.h"
#include <memory>
namespace ModuleESolver
{

    class ESolver_KS_LCAO : public ESolver_KS<double>
    {
    public:
        ESolver_KS_LCAO();
        ~ESolver_KS_LCAO();

        void Init(Input& inp, UnitCell& cell) override;
        void init_after_vc(Input& inp, UnitCell& cell) override;

        double cal_Energy() override;
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

        virtual void othercalculation(const int istep)override;
        ORB_control orb_con;    //Basis_LCAO
        Record_adj RA;
        Local_Orbital_wfc LOWF;
        Local_Orbital_Charge LOC;
        LCAO_Hamilt UHM;
        LCAO_Matrix LM;
        Grid_Technique GridT;

        std::unique_ptr<TwoCenterBundle> two_center_bundle;

        // Temporarily store the stress to unify the interface with PW,
        // because it's hard to seperate force and stress calculation in LCAO.
        // The copy costs memory and time ! 
        // Are there any better way to solve this problem ?
        ModuleBase::matrix scs;
        bool have_force = false;

        void Init_Basis_lcao(ORB_control& orb_con, Input& inp, UnitCell& ucell);

        //--------------common for all calculation, not only scf-------------
        // set matrix and grid integral
        void set_matrix_grid(Record_adj& ra);
        void beforesolver(const int istep);
        //----------------------------------------------------------------------

        /// @brief create ModuleIO::Output_DM object to output density matrix
        ModuleIO::Output_DM create_Output_DM(int is, int iter);

        /// @brief create ModuleIO::Output_DM1 object to output sparse density matrix
        ModuleIO::Output_DM1 create_Output_DM1(int istep);

        /// @brief create ModuleIO::Output_Mat_Sparse object to output sparse density matrix of H, S, T, r
        ModuleIO::Output_Mat_Sparse create_Output_Mat_Sparse(int istep);

        /// @brief check if skip the corresponding output in md calculation
        bool md_skip_out(std::string calculation, int istep, int interval);

#ifdef __EXX
        std::shared_ptr<Exx_LRI_Interface<double>> exd = nullptr;
        std::shared_ptr<Exx_LRI_Interface<std::complex<double>>> exc = nullptr;
        std::shared_ptr<Exx_LRI<double>> exx_lri_double = nullptr;
        std::shared_ptr<Exx_LRI<std::complex<double>>> exx_lri_complex = nullptr;
#endif
    };



}
#endif
