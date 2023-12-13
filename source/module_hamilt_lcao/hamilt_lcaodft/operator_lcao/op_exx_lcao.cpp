#ifdef __EXX
#include "op_exx_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_ri/RI_2D_Comm.h"

namespace hamilt
{

template class OperatorEXX<OperatorLCAO<double, double>>;

template class OperatorEXX<OperatorLCAO<std::complex<double>, double>>;

template class OperatorEXX<OperatorLCAO<std::complex<double>, std::complex<double>>>;

template<typename TK, typename TR>
void OperatorEXX<OperatorLCAO<TK, TR>>::contributeHR()
{

}

// double and complex<double> are the same temperarily
template<>
void OperatorEXX<OperatorLCAO<double, double>>::contributeHk(int ik)
{
    // Peize Lin add 2016-12-03
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
		if(GlobalC::exx_info.info_ri.real_number)
			RI_2D_Comm::add_Hexx(
				kv,
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
                *this->LM->Hexxd,
                *this->LM->ParaV,
                *this->hK);
		else
			RI_2D_Comm::add_Hexx(
				kv,
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
                *this->LM->Hexxc,
                *this->LM->ParaV,
                *this->hK);
    }
}

template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    // Peize Lin add 2016-12-03
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
		if(GlobalC::exx_info.info_ri.real_number)
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				*this->LM->Hexxd,
                *this->LM->ParaV,
                *this->hK);
		else
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				*this->LM->Hexxc,
                *this->LM->ParaV,
                *this->hK);
    }
}

template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHk(int ik)
{
    // Peize Lin add 2016-12-03
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
		if(GlobalC::exx_info.info_ri.real_number)
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				*this->LM->Hexxd,
                *this->LM->ParaV,
                *this->hK);
		else
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				*this->LM->Hexxc,
                *this->LM->ParaV,
                *this->hK);
    }
}

} // namespace hamilt
#endif