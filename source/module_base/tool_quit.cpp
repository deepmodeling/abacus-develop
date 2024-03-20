#include "tool_quit.h"

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

void QUIT(void)
{
	QUIT(0);
}

void QUIT(int ret)
{

#ifdef __NORMAL

#else
    ModuleBase::timer::finish(GlobalV::ofs_running , !GlobalV::MY_RANK);

    ModuleBase::Global_File::close_all_log(GlobalV::MY_RANK);

    //if (GlobalV::MY_RANK==0)
    //{
    //    ModuleBase::Memory::print_all( GlobalV::ofs_running ) ;
    //}
    std::cout<<" See output information in : "<<GlobalV::global_out_dir<<std::endl;
#endif

    exit(ret);
}


void WARNING_QUIT(const std::string &file,const std::string &description)
{
	WARNING_QUIT(file, description, 0);
}

void WARNING_QUIT(const std::string &file,const std::string &description,int ret)
{
#ifdef __NORMAL

		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << "                         NOTICE                           " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

#else

    //std::cout<<" ----------- SOMETHING TO WARN YOU ! -------" << std::endl;
	// 41: red background
	// 42: green background
	// 31: red
	// 32: greem
	// 33: yellow
	// 34: blue
	// 35: zi
	// 36: qing
	// 37: white
	if(GlobalV::COLOUR)
	{
		//printf( "\e[32m%s\e[0m\n", " -------------- SOMETHING TO WARN YOU ! ------------");
		//printf( " \e[32m%s\e[0m\n", description.c_str());
		//printf( " \e[32m%s\e[0m", "CHECK IN FILE : ");
		//printf( "\e[32m%s\e[0m", GlobalV::global_out_dir.c_str());
		//printf( "\e[32m%s\e[0m\n", "warning.log");
		//printf( "\e[32m%s\e[0m\n", " ---------------------------------------------------");
		printf( "[32m%s[0m\n", " -------------- SOMETHING TO WARN YOU ! ------------");
		printf( " [32m%s[0m\n", description.c_str());
		printf( " [32m%s[0m", "CHECK IN FILE : ");
		printf( "[32m%s[0m", GlobalV::global_out_dir.c_str());
		printf( "[32m%s[0m\n", "warning.log");
		printf( "[32m%s[0m\n", " ---------------------------------------------------");
	}
	else
	{
		std::cout << " " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << "                         NOTICE                           " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << " " << std::endl;
		std::cout << " " << description << std::endl;
		std::cout << " CHECK IN FILE : " << GlobalV::global_out_dir << "warning.log" << std::endl;
		std::cout << " " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << "                         NOTICE                           " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;


		GlobalV::ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		GlobalV::ofs_running << "                         NOTICE                           " << std::endl;
		GlobalV::ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		GlobalV::ofs_running << std::endl;
		GlobalV::ofs_running << " " << description << std::endl;
		GlobalV::ofs_running << " CHECK IN FILE : " << GlobalV::global_out_dir << "warning.log" << std::endl;
		GlobalV::ofs_running << std::endl;
		GlobalV::ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		GlobalV::ofs_running << "                         NOTICE                           " << std::endl;
		GlobalV::ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	}

	WARNING(file,description);
    GlobalV::ofs_running<<" Check in file : "<<GlobalV::global_out_dir<<"warning.log"<<std::endl;

#endif

    QUIT(ret);
}

//Check and print warning information for all cores.
//Maybe in the future warning.log should be replaced by error.log.
void CHECK_WARNING_QUIT(const bool error_in, const std::string &file,const std::string &description)
{
#ifdef __NORMAL
// only for UT, do nothing here
#else
	if(error_in)
	{
		//All cores will print inforamtion
		std::cout.clear();
		if(!GlobalV::ofs_running.is_open())
		{
			std::string logfile = GlobalV::global_out_dir + "running_" + GlobalV::CALCULATION + ".log";
			GlobalV::ofs_running.open( logfile.c_str(), std::ios::app );
		}
		if(!GlobalV::ofs_warning.is_open())
		{
			std::string warningfile = GlobalV::global_out_dir + "warning.log";
			GlobalV::ofs_warning.open( warningfile.c_str(), std::ios::app );
		}

		//print error information
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << " ERROR! " << description << std::endl;
		std::cout << " CHECK IN FILE : " << GlobalV::global_out_dir << "warning.log" << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		GlobalV::ofs_running << " ERROR! CHECK IN FILE : " << GlobalV::global_out_dir << "warning.log" << std::endl;
		GlobalV::ofs_warning << std::endl;
		GlobalV::ofs_warning << " ERROR! " << file << ", core " << GlobalV::MY_RANK+1 << ": " << description << std::endl;
		GlobalV::ofs_warning << std::endl;
		exit(1);
	}
#endif
	return;
}

}
