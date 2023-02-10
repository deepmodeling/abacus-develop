#include "op_exx_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/hamilt_lcaodft/global_fp.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace hamilt
{

template class OperatorEXX<OperatorLCAO<double>>;

template class OperatorEXX<OperatorLCAO<std::complex<double>>>;

template<typename T>
void OperatorEXX<OperatorLCAO<T>>::contributeHR()
{

}

// double and complex<double> are the same temperarily
template<>
void OperatorEXX<OperatorLCAO<double>>::contributeHk(int ik)
{
#ifdef __EXX
    // Peize Lin add 2016-12-03
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
		if(GlobalC::exx_info.info_ri.real_number)
			RI_2D_Comm::add_Hexx(
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				GlobalC::exx_lri_double.Hexxs,
				*this->LM->ParaV,
				*this->LM);
		else
			RI_2D_Comm::add_Hexx(
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				GlobalC::exx_lri_complex.Hexxs,
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
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
		if(GlobalC::exx_info.info_ri.real_number)
			RI_2D_Comm::add_Hexx(
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				GlobalC::exx_lri_double.Hexxs,
				*this->LM->ParaV,
				*this->LM);
		else
			RI_2D_Comm::add_Hexx(
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				GlobalC::exx_lri_complex.Hexxs,
				*this->LM->ParaV,
				*this->LM);
    }
#endif
}
} // namespace hamilt