#include "tool_quit.h"
#include <cstdlib>
#ifdef __MPI
#include "mpi.h"
#endif

#ifdef __NORMAL

#else
#include "global_variable.h"
#include "global_file.h"
#include "timer.h"
#include "memory.h"
#endif

namespace ModuleBase
{
namespace
{
std::string g_quit_out_dir;
std::string g_quit_calculation;

[[noreturn]] void quit_impl(const int ret, const bool print_timer)
{
#ifdef __NORMAL
    (void)print_timer;
#else
    if (print_timer)
    {
        ModuleBase::timer::finish(GlobalV::ofs_running, GlobalV::MY_RANK == 0, false);
        std::cout << " See output information in : " << g_quit_out_dir << std::endl;
    }
    ModuleBase::Global_File::close_all_log(GlobalV::MY_RANK);
#endif

    std::exit(ret);
}
}

void set_quit_out_dir(const std::string& dir)
{
    g_quit_out_dir = dir;
}

const std::string& get_global_out_dir()
{
    return g_quit_out_dir;
}

void set_quit_calculation(const std::string& calculation)
{
    g_quit_calculation = calculation;
}

//==========================================================
// GLOBAL FUNCTION :
// NAME : WARNING( write information into GlobalV::ofs_warning)
// NAME : QUIT( exit the running program)
// NAME : WARNING_QUIT( write information into
// 		  GlobalV::ofs_warning , and then quit)
//==========================================================
void WARNING(const std::string &file,const std::string &description)
{
#ifdef __NORMAL

#else

    if (GlobalV::MY_RANK==0)
    {
//		std::cout << "\n "<<file<<"  Warning : "<<description<<std::endl;
        GlobalV::ofs_warning << " " << file <<"  warning : "<< description<<std::endl;
    }
#endif
    return;
}

void QUIT()
{
	QUIT(0);
}

void QUIT(int ret)
{
    quit_impl(ret, true);
}


void WARNING_QUIT(const std::string &file,const std::string &description)
{
    WARNING_QUIT(file, description, 1);
}

void WARNING_QUIT(const std::string &file,const std::string &description,int ret)
{
#ifdef __NORMAL

	std::cout << " ---------------------------------------------------------" << std::endl;
	std::cout << "                         !ERROR!                          " << std::endl;
	std::cout << " ---------------------------------------------------------" << std::endl;
    std::cout << " For detailed manual of ABACUS, please see the website" << std::endl;
    std::cout << " https://abacus.deepmodeling.com" << std::endl;
	std::cout << " For any questions, propose issues on the website" << std::endl;
    std::cout << " https://github.com/deepmodeling/abacus-develop/issues" << std::endl;

#else
	std::cout << " " << std::endl;
	std::cout << " ---------------------------------------------------------" << std::endl;
	std::cout << "                         !ERROR!                          " << std::endl;
	std::cout << " ---------------------------------------------------------" << std::endl;
	std::cout << " " << std::endl;
	std::cout << " " << description << std::endl;
	std::cout << " " << std::endl;
	std::cout << " CHECK IN FILE : " << g_quit_out_dir << "warning.log" << std::endl;
	std::cout << " " << std::endl;
	std::cout << " For detailed manual of ABACUS, please see the website" << std::endl;
	std::cout << " https://abacus.deepmodeling.com" << std::endl;
	std::cout << " For any questions, propose issues on the website" << std::endl;
	std::cout << " https://github.com/deepmodeling/abacus-develop/issues" << std::endl;
	std::cout << " ---------------------------------------------------------" << std::endl;
	std::cout << "                         !ERROR!                          " << std::endl;
	std::cout << " ---------------------------------------------------------" << std::endl;


	GlobalV::ofs_running << " ---------------------------------------------------------" << std::endl;
	GlobalV::ofs_running << "                         !ERROR!                          " << std::endl;
	GlobalV::ofs_running << " ---------------------------------------------------------" << std::endl;
	GlobalV::ofs_running << std::endl;
	GlobalV::ofs_running << " " << description << std::endl;
	GlobalV::ofs_running << std::endl;
	GlobalV::ofs_running << " CHECK IN FILE : " << g_quit_out_dir << "warning.log" << std::endl;
	GlobalV::ofs_running << std::endl;
	GlobalV::ofs_running << " For detailed manual of ABACUS, please see the website" << std::endl;
	GlobalV::ofs_running << " https://abacus.deepmodeling.com" << std::endl;
	GlobalV::ofs_running << " For any questions, propose issues on the website" << std::endl;
	GlobalV::ofs_running << " https://github.com/deepmodeling/abacus-develop/issues" << std::endl;
	GlobalV::ofs_running << " ---------------------------------------------------------" << std::endl;
	GlobalV::ofs_running << "                         !ERROR!                          " << std::endl;
	GlobalV::ofs_running << " ---------------------------------------------------------" << std::endl;

	WARNING(file,description);
    GlobalV::ofs_running<<" Check in file : "<<g_quit_out_dir<<"warning.log"<<std::endl;

#endif

    quit_impl(ret, false);
}

//Check and print warning information for all cores.
//Maybe in the future warning.log should be replaced by error.log.
void CHECK_WARNING_QUIT(const bool error_in, const std::string &file, const std::string &description)
{
#ifdef __NORMAL
    if(error_in) std::cout << description << std::endl;
#else
	if(error_in)
	{
		//All cores will print inforamtion
		std::cout.clear();
		if(!GlobalV::ofs_running.is_open())
		{
			std::string logfile = g_quit_out_dir + "running_" + g_quit_calculation + ".log";
			GlobalV::ofs_running.open( logfile.c_str(), std::ios::app );
		}
		if(!GlobalV::ofs_warning.is_open())
		{
			std::string warningfile = g_quit_out_dir + "warning.log";
			GlobalV::ofs_warning.open( warningfile.c_str(), std::ios::app );
		}

		//print error information
		std::cout << " ---------------------------------------------------------" << std::endl;
		std::cout << " ERROR! " << description << std::endl;
		std::cout << " CHECK IN FILE : " << g_quit_out_dir << "warning.log" << std::endl;
		std::cout << " ---------------------------------------------------------" << std::endl;
		GlobalV::ofs_running << " ERROR! CHECK IN FILE : " << g_quit_out_dir << "warning.log" << std::endl;
		GlobalV::ofs_warning << std::endl;
		GlobalV::ofs_warning << " ERROR! " << file << ", core " << GlobalV::MY_RANK+1 << ": " << description << std::endl;
		GlobalV::ofs_warning << std::endl;
		exit(1);
	}
#endif
	return;
}

}
