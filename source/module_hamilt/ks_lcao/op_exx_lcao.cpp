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
#ifdef __EXX
    // Peize Lin add 2016-12-03
    auto &exx_lri = GlobalC::exx_lri_double;
    auto &exx_info = GlobalC::exx_info;
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
        //if(Exx_Info::Hybrid_Type::HF  == GlobalC::exx_info.info_global.hybrid_type)	// Peize Lin delete 2022.11.13
        //{
        //    exx_info.info_global.hybrid_alpha = 1.0;
        //}
        RI_2D_Comm::add_Hexx(
            ik,
            exx_info.info_global.hybrid_alpha,
            exx_lri.Hexxs,
            *this->LM->ParaV,
            *this->LM);
    }
#endif
}

template<>
void OperatorEXX<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
#ifdef __EXX
    // Peize Lin add 2016-12-03
    auto &exx_lri = GlobalC::exx_lri_complex;
    auto &exx_info = GlobalC::exx_info;
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
        //if(Exx_Info::Hybrid_Type::HF  == GlobalC::exx_info.info_global.hybrid_type)	// Peize Lin delete 2022.11.13
        //{
        //    exx_info.info_global.hybrid_alpha = 1.0;
        //}
        RI_2D_Comm::add_Hexx(
            ik,
            exx_info.info_global.hybrid_alpha,
            exx_lri.Hexxs,
            *this->LM->ParaV,
            *this->LM);
    }
#endif
}
} // namespace hamilt