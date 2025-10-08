#include "source_lcao/module_ri/setup_exx.h"

void Exx_NAO::init0()
{
#ifdef __EXX
    // 1. currently this initialization must be put in constructor rather than `before_all_runners()`
    //  because the latter is not reused by ESolver_LCAO_TDDFT,
    //  which cause the failure of the subsequent procedure reused by ESolver_LCAO_TDDFT
    // 2. always construct but only initialize when if(cal_exx) is true
    //  because some members like two_level_step are used outside if(cal_exx)
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exd = std::make_shared<Exx_LRI_Interface<TK, double>>(GlobalC::exx_info.info_ri);
    }
    else
    {
        this->exc = std::make_shared<Exx_LRI_Interface<TK, std::complex<double>>>(GlobalC::exx_info.info_ri);
    }
#endif
}

void Exx_NAO::before_runner()
{
#ifdef __EXX
    if (inp.calculation == "scf" || inp.calculation == "relax" || inp.calculation == "cell-relax"
        || inp.calculation == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            if (inp.init_wfc != "file")
            { // if init_wfc==file, directly enter the EXX loop
                XC_Functional::set_xc_first_loop(ucell);
            }

            // initialize 2-center radial tables for EXX-LRI
            if (GlobalC::exx_info.info_ri.real_number)
            {
                this->exd->init(MPI_COMM_WORLD, ucell, this->kv, orb_);
                this->exd->exx_before_all_runners(this->kv, ucell, this->pv);
            }
            else
            {
                this->exc->init(MPI_COMM_WORLD, ucell, this->kv, orb_);
                this->exc->exx_before_all_runners(this->kv, ucell, this->pv);
            }
        }
    }
#endif
}
