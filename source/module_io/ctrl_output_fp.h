#ifndef CTRL_OUTPUT_FP_H 
#define CTRL_OUTPUT_FP_H 

namespace ModuleIO
{
	template <typename TK, typename TR>
		void ctrl_output_fp(UnitCell& ucell, 
				elecstate::ElecStateLCAO<TK>* pelec, 
				const int istep);
}
#endif
