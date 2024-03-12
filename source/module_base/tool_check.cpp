#include "tool_check.h"

#include "tool_quit.h"

namespace ModuleBase
{

void CHECK_NAME(std::ifstream &ifs,const std::string &name_in,bool quit)
{
    std::string name;
    ifs >> name;
    if ( name != name_in)
    {
		if(quit)
		{
			//GlobalV::ofs_warning << "\n name = " <<name;
			//GlobalV::ofs_warning << "\n should be = " << name_in;
			std::cout << "\n name = " <<name;
			std::cout << "\n should be = " << name_in;
        	WARNING_QUIT("CHECK_NAME","Some parameter name is wrong!");
		}
		else
		{
        	std::cout <<"\n Can not match : "<<name<<"(readin)  "<<name_in<<std::endl;
		}
    }
    return;
}

void CHECK_INT(std::ifstream &ifs,const int &v,bool quit)
{
	int v_in;
	ifs >> v_in;
	if( v!= v_in)
	{
		if(quit)
		{
			std::cout << "\n value = " << v_in;
			std::cout << "\n should be = " << v;
			WARNING_QUIT("CHECK_INT","Some parameter name is wrong!");
		}
		else
		{
			std::cout <<"\n Can not match well: "<<v_in<<"(readin)  "<<v<<std::endl;
		}
	}
	return;
}

void CHECK_DOUBLE(std::ifstream &ifs,const double &v,bool quit)
{
	const double tiny = 1.0e-5;
	double v_in;
	ifs >> v_in;
	if( fabs(v - v_in) > tiny )
	{
		if(quit)
		{
			std::cout << " read in value = " << v_in << std::endl;
			std::cout << " the value should be = " << v << std::endl;
			WARNING_QUIT("CHECK_DOUBLE","the name of parameter wrong!");
		}
		else
		{
			std::cout <<" can not match well (1.0e-5): "<< v_in <<"(readin)  "<<v<<std::endl;
		}
	}
	return;
}

void CHECK_STRING(std::ifstream &ifs,const std::string &v,bool quit)
{
	std::string v_in;
	ifs >> v_in;
	if( v_in != v )
	{
		if(quit)
		{
			std::cout << " read in value = " << v_in << std::endl;
			std::cout << " the value should be = " << v << std::endl;
			WARNING_QUIT("CHECK_STRING","the name of parameter wrong!");
		}
		else
		{
			std::cout <<" can not match well : "<<v_in<<"(readin)  "<<v<<std::endl;
		}
	}
	return;
}

}
