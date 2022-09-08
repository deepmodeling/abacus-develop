#include "esolver.h"
#include "esolver_ks_pw.h"
#include "esolver_sdft_pw.h"
#include "esolver_ks_lcao.h"
#include "esolver_ks_lcao_tddft.h"
#include "esolver_of.h"
#include "esolver_lj.h"
#include "esolver_dp.h"
#include "module_md/MD_parameters.h"

namespace ModuleESolver
{
    void ESolver::printname()
    {
        std::cout << classname << std::endl;
    }

    string determine_type()
    {
        string esolver_type;
        if (GlobalV::BASIS_TYPE == "pw")
        {
            if (GlobalV::CALCULATION.substr(0, 3) == "sto")
                esolver_type = "sdft_pw";
            else
                esolver_type = "ksdft_pw";
        }
#ifdef __LCAO
        else if (GlobalV::BASIS_TYPE == "lcao_in_pw")
        {
            if (GlobalV::CALCULATION.substr(0, 3) == "sto")
                esolver_type = "sdft_pw";
            else
                esolver_type = "ksdft_pw";            
        }
        else if (GlobalV::BASIS_TYPE == "lcao")
        {
            esolver_type = "ksdft_lcao";
            if (INPUT.tddft == 1)
                esolver_type = "ksdft_lcao_tddft";
        }
#endif
        if(INPUT.mdp.md_ensolver == "LJ")
        {
            esolver_type = "lj_pot";
        }
        if(INPUT.mdp.md_ensolver == "DP")
        {
            esolver_type = "dp_pot";
        }

        GlobalV::ofs_running << " The esolver type has been set to : " << esolver_type << std::endl;
        return esolver_type;
    }

    //Some API to operate E_Solver
    void init_esolver(ESolver*& p_esolver)
    {
        //determine type of esolver based on INPUT information
        string esolver_type = determine_type();

        //initialize the corresponding Esolver child class
        if (esolver_type == "ksdft_pw")
        {
            p_esolver = new ESolver_KS_PW();
        }
        else if (esolver_type == "ksdft_lcao")
        {
            p_esolver = new ESolver_KS_LCAO();
        }
        else if (esolver_type == "ksdft_lcao_tddft")
        {
            p_esolver = new ESolver_KS_LCAO_TDDFT();
        }
        else if (esolver_type == "sdft_pw")
        {
            p_esolver = new ESolver_SDFT_PW();
        }
        //  else if(esolver_type == "ofdft")
        //  {
        //      p_esolver = new OFDFT();
        //  }
        else if (esolver_type == "lj_pot")
        {
            p_esolver = new ESolver_LJ();
        }
        else if (esolver_type == "dp_pot")
        {
            p_esolver = new ESolver_DP();
        }
    }

    void clean_esolver(ESolver*& pesolver)
    {
        delete pesolver;
    }

}
