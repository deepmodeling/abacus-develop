#include "source_lcao/setup_exx.h"

template <typename TK>
Exx_NAO<TK>::Exx_NAO(){}

template <typename TK>
Exx_NAO<TK>::~Exx_NAO(){}


template <typename TK>
void Exx_NAO<TK>::init()
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

template <typename TK>
void Exx_NAO<TK>::before_runner(
		UnitCell& ucell, // unitcell
		K_Vectors &kv, // k points
        const LCAO_Orbitals &orb, // orbital info 
        const Parallel_Orbitals &pv, // parallel orbitals
		const Input_para& inp)
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
                this->exd->init(MPI_COMM_WORLD, ucell, kv, orb);
                this->exd->exx_before_all_runners(kv, ucell, pv);
            }
            else
            {
                this->exc->init(MPI_COMM_WORLD, ucell, kv, orb);
                this->exc->exx_before_all_runners(kv, ucell, pv);
            }
        }
    }
#endif
}

template class Exx_NAO<double>;
template class Exx_NAO<std::complex<double>>;
