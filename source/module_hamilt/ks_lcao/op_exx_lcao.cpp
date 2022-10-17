#include "op_exx_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_xc/xc_functional.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"

namespace hamilt
{

template class OperatorEXX<OperatorLCAO<double>>;

template class OperatorEXX<OperatorLCAO<std::complex<double>>>;

template<typename T>
void OperatorEXX<OperatorLCAO<T>>::contributeHR()
{

}

template<>
void OperatorEXX<OperatorLCAO<double>>::contributeHk(int ik)
{
#ifdef __MPI //liyuanbo 2022/2/23
    // Peize Lin add 2016-12-03
    auto &exx_lcao = GlobalC::exx_lcao;
    auto &exx_global = GlobalC::exx_global;
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
        if( Exx_Global::Hybrid_Type::HF == exx_lcao.info.hybrid_type ) //HF
        {
            exx_lcao.add_Hexx(ik, 1, *this->LM);
        }
        else if( Exx_Global::Hybrid_Type::PBE0 == exx_lcao.info.hybrid_type )			// PBE0
        {
            exx_lcao.add_Hexx(ik, exx_global.info.hybrid_alpha, *this->LM);
        }
        else if( Exx_Global::Hybrid_Type::SCAN0 == exx_lcao.info.hybrid_type )			// SCAN0
        {
            exx_lcao.add_Hexx(ik, exx_global.info.hybrid_alpha, *this->LM);
        }
        else if( Exx_Global::Hybrid_Type::HSE  == exx_lcao.info.hybrid_type )			// HSE
        {
            exx_lcao.add_Hexx(ik, exx_global.info.hybrid_alpha, *this->LM);
        }
    }
#endif
}

template<>
void OperatorEXX<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
#ifdef __MPI //liyuanbo 2022/2/23
    // Peize Lin add 2016-12-03
    auto &exx_lcao = GlobalC::exx_lcao;
    auto &exx_global = GlobalC::exx_global;
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
        if( Exx_Global::Hybrid_Type::HF  == exx_lcao.info.hybrid_type )				// HF
        {
            exx_lcao.add_Hexx(ik,1, *this->LM);
        }
        else if( Exx_Global::Hybrid_Type::PBE0  == exx_lcao.info.hybrid_type )			// PBE0
        {
            exx_lcao.add_Hexx(ik,exx_global.info.hybrid_alpha, *this->LM);
        }
        else if( Exx_Global::Hybrid_Type::SCAN0  == exx_lcao.info.hybrid_type )			// SCAN0
        {
            exx_lcao.add_Hexx(ik,exx_global.info.hybrid_alpha, *this->LM);
        }
        else if( Exx_Global::Hybrid_Type::HSE  == exx_lcao.info.hybrid_type )			// HSE
        {
            exx_lcao.add_Hexx(ik,exx_global.info.hybrid_alpha, *this->LM);
        }
    }
#endif
}

}