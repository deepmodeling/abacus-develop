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
                this->Hexxd == nullptr ? *this->LM->Hexxd : *this->Hexxd,
                *this->LM->ParaV,
                *this->hK);
		else
			RI_2D_Comm::add_Hexx(
				kv,
				ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
                this->Hexxc == nullptr ? *this->LM->Hexxc : *this->Hexxc,
                *this->LM->ParaV,
                *this->hK);
    }
}

template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    // std::cout << "\n\n\n******\n begin contributeHk() double c d \n******\n\n\n";
    // Peize Lin add 2016-12-03
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
		if(GlobalC::exx_info.info_ri.real_number)
        {
            // std::cout << "\n\n\n******\n begin contributeHk() double \n******\n\n\n";
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				this->Hexxd == nullptr ? *this->LM->Hexxd : *this->Hexxd,
                *this->LM->ParaV,
                *this->hK);
            // std::cout << "\n\n\n******\n successfully run add_Hexx() double \n******\n\n\n";
        }
		else
        {
            if(this->symbol == 1) std::cout << "\n\n\n******\n begin contributeHk() complex \n******\n\n\n";    ///// run here 
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				this->Hexxc == nullptr ? *this->LM->Hexxc : *this->Hexxc,
                *this->LM->ParaV,
                *this->hK,
                this->symbol);
            if(this->symbol == 1) std::cout << "\n\n\n******\n successfully run add_Hexx() complex \n******\n\n\n";    ////// and don't arrive here
        }
    }
}

template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHk(int ik)
{
    std::cout << "\n\n\n******\n begin contributeHk() \n******\n\n\n";
    // Peize Lin add 2016-12-03
    if(XC_Functional::get_func_type()==4 || XC_Functional::get_func_type()==5)
    {
        std::cout << "\n\n\n******\n begin contributeHk() double \n******\n\n\n";
		if(GlobalC::exx_info.info_ri.real_number)
        {
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				this->Hexxd == nullptr ? *this->LM->Hexxd : *this->Hexxd,
                *this->LM->ParaV,
                *this->hK);
        std::cout << "\n\n\n******\n successfully run add_Hexx() double \n******\n\n\n";
        }
		// else
        //     RI_2D_Comm::add_Hexx(
        //         kv,
        //         ik,
		// 		GlobalC::exx_info.info_global.hybrid_alpha,
		// 		this->Hexxc == nullptr ? *this->LM->Hexxc : *this->Hexxc,
        //         *this->LM->ParaV,
        //         *this->hK);
		else
        {
            std::cout << "\n\n\n******\n begin contributeHk() complex \n******\n\n\n";
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
				GlobalC::exx_info.info_global.hybrid_alpha,
				this->Hexxc == nullptr ? *this->LM->Hexxc : *this->Hexxc,
                *this->LM->ParaV,
                *this->hK);
            std::cout << "\n\n\n******\n successfully run add_Hexx() complex \n******\n\n\n";
        }
    }
}

} // namespace hamilt
#endif