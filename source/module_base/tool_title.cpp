#include "tool_title.h"

#ifdef __NORMAL
#else
#include "global_variable.h"
#endif

namespace ModuleBase
{
//==========================================================
// GLOBAL FUNCTION :
// NAME : TITLE( title for each function )
//==========================================================
void TITLE(const std::string &class_function_name)
{
    return;//no output

#ifdef __NORMAL
    std::cout<<" ==> "<<class_function_name<<std::endl;
#else
	if(GlobalV::ofs_running) // mohan add 2009-08-25 in case the function called before allocate GlobalV::ofs_running
	{
   		GlobalV::ofs_running<<" ==> "<<class_function_name<<std::endl;
	}
#endif
}

void TITLE(const std::string &class_name,const std::string &function_name)
{
    return;//no output
#ifdef __NORMAL
    std::cout<<" ==> "<<class_name<<"::"<<function_name<<std::endl;
#else
	if(GlobalV::ofs_running) // mohan add 2009-08-25 in case the function called before allocate GlobalV::ofs_running
	{
   		GlobalV::ofs_running<<" ==> "<<class_name<<"::"<<function_name<<std::endl;
	}
#endif
    return;
}

void TITLE(std::ofstream &ofs,const std::string &class_name,const std::string &function_name)
{
    return;// no output
#ifdef __NORMAL
    std::cout<<"\n\n ==> "<<class_name<<"::"<<function_name<<std::endl;
#else
	if(GlobalV::ofs_running)
	{
    	ofs<<" ==> "<<class_name<<"::"<<function_name<<std::endl;
	}
#endif
    return;
}

}