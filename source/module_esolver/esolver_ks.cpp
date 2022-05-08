#include "esolver_ks.h"
#include <iostream>
#include <algorithm>
#include "time.h"

//--------------Temporary----------------
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
#include "../src_pw/charge_broyden.h"
#include "../module_base/timer.h"
//---------------------------------------

namespace ModuleESolver
{

ESolver_KS::ESolver_KS()
{
    classname = "ESolver_KS";
    basisname = "PLEASE ADD BASISNAME FOR CURRENT ESOLVER.";
    diag_ethr = GlobalV::PW_DIAG_THR; 
    scf_thr = GlobalV::SCF_THR; 
    drho = 0.0;
    maxniter = GlobalV::SCF_NMAX;
    niter = maxniter;
    out_freq_elec = GlobalV::OUT_FREQ_ELEC;
}

void ESolver_KS:: hamilt2density(const int istep, const int iter, const double ethr)
{
    ModuleBase::timer:: tick(this->classname,"hamilt2density");
    //Temporarily, before HSolver is constructed, it should be overrided by
    //LCAO, PW, SDFT and TDDFT.
    //After HSolver is constructed, LCAO, PW, SDFT should delete their own
    //hamilt2density() and use:
    //this->phsol->solve(this->phamilt, this->pes, this->wf, ETHR);
    ModuleBase::timer:: tick(this->classname,"hamilt2density");
}


void ESolver_KS:: Run(const int istep, UnitCell_pseudo& cell)
{
    if (!(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md"
        || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
#ifdef __MPI
        || Exx_Global::Hybrid_Type::No!=GlobalC::exx_global.info.hybrid_type 
#endif
        )
    {
        this->othercalculation(istep);
    }
    else
    {
        ModuleBase::timer::tick(this->classname, "Run");

        this->printhead(); //print the headline on the screen.
        this->beforescf(istep); //Something else to do before the iter loop

        bool firstscf = true;
        bool conv_elec = false;
        this->niter = this->maxniter;
        int iter = 1;
        for (iter = 1; iter <= this->maxniter; ++iter)
        {
            writehead(GlobalV::ofs_running, istep, iter);
            clock_t iterstart, iterend;
            iterstart = std::clock();
            set_ethr(istep, iter);
            eachiterinit(istep, iter);

            this->hamilt2density(istep, iter, this->diag_ethr);

            // double drho = this->estate.caldr2(); 
            // EState should be used after it is constructed.
            drho = GlobalC::CHR.get_drho();
            double hsolver_error = 0.0;
            if (firstscf)
            {
                firstscf = false;
                hsolver_error = this->diag_ethr * std::max(1.0, GlobalC::CHR.nelec);
                // The error of HSolver is larger than drho, so a more precise HSolver should be excuconv_elected.
                if (hsolver_error > drho)
                {
                    reset_diagethr(GlobalV::ofs_running, hsolver_error);
                    this->hamilt2density(istep, iter, this->diag_ethr);
                    drho = GlobalC::CHR.get_drho();
                    hsolver_error = this->diag_ethr * std::max(1.0, GlobalC::CHR.nelec);
                }
            }

            conv_elec = (drho < this->scf_thr);

            // If drho < hsolver_error in the first iter or drho < scf_thr, we do not change rho.
            if (drho < hsolver_error || conv_elec)
            {
                if (drho < hsolver_error)    GlobalV::ofs_warning << " drho < hsolver_error, keep charge density unchanged." << std::endl;
            }
            else
            {
                //charge mixing
                //conv_elec = this->estate.mix_rho();
                GlobalC::CHR.mix_rho(iter);
            }

            // Hamilt should be used after it is constructed.
            // this->phamilt->update(conv_elec);
            updatepot(istep, iter, conv_elec);
            eachiterfinish(iter, conv_elec);
            iterend = std::clock();
            double duration = double(iterend - iterstart) / CLOCKS_PER_SEC;
            printiter(conv_elec, iter, drho, duration, diag_ethr);
            if (conv_elec)
            {
                this->niter = iter;
                break;
            }
        }
        afterscf(iter, conv_elec);

        ModuleBase::timer::tick(this->classname, "Run");
    }

    return;
};

//<Temporary> It should be a function of Diag_H class in the future.
void ESolver_KS:: set_ethr(const int istep, const int iter)
{
//It is too complex now and should be modified.
    if(iter == 1)
    {
        if (abs(this->diag_ethr - 1.0e-2) < 1.0e-10)
        {
            if (GlobalC::pot.init_chg == "file")
            {
                //======================================================
                // if you think that the starting potential is good
                // do not spoil it with a louly first diagonalization:
                // set a strict this->diag_ethr in the input file ()diago_the_init
                //======================================================
                this->diag_ethr = 1.0e-5;
            }
            else
            {
                //=======================================================
                // starting atomic potential is probably far from scf
                // don't waste iterations in the first diagonalization
                //=======================================================
                this->diag_ethr = 1.0e-2;
            }
        }
        if(GlobalV::FINAL_SCF) this->diag_ethr = 1.0e-2;
        
        if(GlobalV::CALCULATION=="md"||GlobalV::CALCULATION=="relax"||GlobalV::CALCULATION=="cell-relax")
        {
            this->diag_ethr = std::max(this->diag_ethr, INPUT.pw_diag_thr);
        }
        
    }
    else
    {
        if (iter == 2)
        {
            this->diag_ethr = 1.e-2;
        }
		this->diag_ethr = std::min( this->diag_ethr, 0.1*this->drho/ std::max(1.0, GlobalC::CHR.nelec));

    }
    if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") 
    {
        this->diag_ethr = 0.0;
    }
}

void ESolver_KS:: printhead()
{
    std::cout << " " << std::setw(7) << "ITER";
    if(GlobalV::NSPIN==2)
    {
        std::cout << std::setw(10) << "TMAG";
        std::cout << std::setw(10) << "AMAG";
    }
    std::cout << std::setw(15) << "ETOT(eV)";
    std::cout << std::setw(15) << "EDIFF(eV)";
    std::cout << std::setw(11) << "DRHO";
	std::cout << std::setw(11) << "TIME(s)" << std::endl;
}

void ESolver_KS::printiter(const bool conv_elec, const int iter, const double drho, const double duration, const double ethr)
{
    GlobalC::en.print_etot(conv_elec, iter, drho, duration, ethr);
}

void ESolver_KS:: writehead(std::ofstream &ofs_running, const int istep, const int iter)
{
    ofs_running
        << "\n "
        << this->basisname
        << " ALGORITHM --------------- ION=" << std::setw(4) << istep + 1
        << "  ELEC=" << std::setw(4) << iter
        << "--------------------------------\n";
}

void ESolver_KS:: reset_diagethr(std::ofstream &ofs_running, const double hsover_error)
{
    ofs_running << " Notice: Threshold on eigenvalues was too large.\n";
    ModuleBase::WARNING("scf","Threshold on eigenvalues was too large.");
    ofs_running << " hsover_error=" << hsover_error << " > DRHO=" << drho << std::endl;
    ofs_running << " Origin diag_ethr = " << this->diag_ethr << std::endl;
    this->diag_ethr =0.1 *  drho / GlobalC::CHR.nelec; 
    ofs_running << " New    diag_ethr = " << this->diag_ethr << std::endl;
}

int ESolver_KS:: getniter()
{
    return this->niter;
}


}
