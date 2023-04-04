#include "tool_title.h"

#ifdef __NORMAL
#else
#include "global_variable.h"
#endif

#include "global_function.h"
#include "timer.h"

namespace ModuleBase
{
//==========================================================
// GLOBAL FUNCTION :
// NAME : TITLE( title for each function )
//==========================================================

void TITLE(const std::string &class_name,const std::string &function_name,const bool disable)
{
    if (disable)
    {
        return;//no output
    }
#ifdef __NORMAL
    std::cout<<" ==> "<<class_name<<"::"<<function_name<<"\t"
		   <<ModuleBase::GlobalFunc::MemAvailable()/1024.0/1024<<" GB\t"
		   <<ModuleBase::timer::print_until_now()<<" s"<<std::endl;
#else
	if(GlobalV::ofs_running) // mohan add 2009-08-25 in case the function called before allocate GlobalV::ofs_running
	{
   		GlobalV::ofs_running<<" ==> "<<class_name<<"::"<<function_name<<"\t"
		   <<ModuleBase::GlobalFunc::MemAvailable()/1024.0/1024<<" GB\t"
		   <<ModuleBase::timer::print_until_now()<<" s"<<std::endl;
	}
#endif
}

void TITLE(std::ofstream &ofs,const std::string &class_name,const std::string &function_name,const bool disable)
{
    if (disable)
    {
        return;//no output
    }
#ifdef __NORMAL
    std::cout<<"\n\n ==> "<<class_name<<"::"<<function_name<<"\t"
		   <<ModuleBase::GlobalFunc::MemAvailable()/1024.0/1024<<" GB\t"
		   <<ModuleBase::timer::print_until_now()<<" s"<<std::endl;
#else
	if(GlobalV::ofs_running)
	{
    	ofs<<" ==> "<<class_name<<"::"<<function_name<<"\t"
		   <<ModuleBase::GlobalFunc::MemAvailable()/1024.0/1024<<" GB\t"
		   <<ModuleBase::timer::print_until_now()<<" s"<<std::endl;
	}
#endif
}

}
