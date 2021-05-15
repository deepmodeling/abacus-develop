#include "tool_title.h"

#ifdef __NORMAL
#else
#include "global_variable.h"
#endif

//==========================================================
// GLOBAL FUNCTION :
// NAME : TITLE( title for each function )
//==========================================================
void TITLE(const string &class_function_name)
{
	return;//no output

#ifdef __NORMAL
    cout<<" ==> "<<class_function_name<<endl;
#else
	if(ofs_running) // mohan add 2009-08-25 in case the function called before allocate ofs_running
	{
   		ofs_running<<" ==> "<<class_function_name<<endl;
	}
#endif
}

void TITLE(const string &class_name,const string &function_name)
{
	return;//no output
#ifdef __NORMAL
    cout<<" ==> "<<class_name<<"::"<<function_name<<endl;
#else
	if(ofs_running) // mohan add 2009-08-25 in case the function called before allocate ofs_running
	{
   		ofs_running<<" ==> "<<class_name<<"::"<<function_name<<endl;
	}
#endif
    return;
}

void TITLE(ofstream &ofs,const string &class_name,const string &function_name)
{
	return;// no output
#ifdef __NORMAL
    cout<<"\n\n ==> "<<class_name<<"::"<<function_name<<endl;
#else
	if(ofs_running)
	{
    	ofs<<" ==> "<<class_name<<"::"<<function_name<<endl;
	}
#endif
    return;
}
