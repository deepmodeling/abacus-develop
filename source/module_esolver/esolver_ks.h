#ifndef ESOLVER_KS_H
#define ESOLVER_KS_H
#include "./esolver_fp.h"
#include "string.h"
#include "fstream"
// #include "estates.h"
// #include "h2e.h"
namespace ModuleESolver
{

class ESolver_KS: public ESolver_FP
{
    public:
        ESolver_KS();
        // HSolver* phsol;
        double diag_ethr; // diag threshold
        double scf_thr;   // scf threshold
        double drho;      // the difference between rho_in (before HSolver) and rho_out (After HSolver)
        int maxniter;     // maximum iter steps for scf
        int niter;        // iter steps actually used in scf
        int out_freq_elec;// frequency for output
        
        virtual void Run(const int istep, UnitCell_pseudo& cell) override;

        // calculate electron density from a specific Hamiltonian
        virtual void hamilt2density(const int istep, const int iter, const double ethr);
        // get
        virtual int getniter() override;

    protected:
        // Something to do before SCF iterations.
        virtual void beforescf(int istep){}; 
        // Something to do before hamilt2density function in each iter loop.
        virtual void eachiterinit(const int istep, const int iter){}; 
        // Something to do after hamilt2density function in each iter loop.
        virtual void eachiterfinish(const int iter, const bool conv){}; 
        // Something to do after SCF iterations when SCF is converged or comes to the max iter step.
        virtual void afterscf(const int iter, const bool conv){};
        // <Temporary> It should be replaced by a function in Hamilt Class
        virtual void updatepot(const int istep, const int iter, const bool conv){};
        

    //TOOLS:
    protected:
        // Set ethr for hsolver
        void set_ethr(const int istep, const int iter);
        // Print the headline on the screen:
        // ITER   ETOT(eV)       EDIFF(eV)      DRHO    TIME(s) 
        void printhead();
        // Print inforamtion in each iter
        // G1    -3.435545e+03  0.000000e+00   3.607e-01  2.862e-01
        void printiter(const bool conv, const int iter, const double drho, const double duration, const double ethr);
        // Write the headline in the running_log file
        // "PW/LCAO" ALGORITHM --------------- ION=   1  ELEC=   1--------------------------------
        void writehead(std::ofstream &ofs_running, const int istep, const int iter);
        void reset_diagethr(std::ofstream &ofs_running, const double hsover_error);



    protected:
        std::string basisname; //PW or LCAO

};
}
#endif