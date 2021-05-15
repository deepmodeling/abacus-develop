#include "tool_quit.h"

#ifdef __NORMAL

#else
#include "global_variable.h"
#include "global_file.h"
#include "../src_parallel/parallel_common.h"
#include "timer.h"
#include "memory.h"
#endif

//==========================================================
// GLOBAL FUNCTION :
// NAME : WARNING( write information into ofs_warning)
// NAME : QUIT( exit the running program)
// NAME : WARNING_QUIT( write information into
// 		  ofs_warning , and then quit)
//==========================================================
void WARNING(const string &file,const string &description)
{
#ifdef __NORMAL

#else

    if (MY_RANK==0)
    {
//		cout << "\n "<<file<<"  Warning : "<<description<<endl;
        ofs_warning << " " << file <<"  warning : "<< description<<endl;
    }
#endif
    return;
}

void QUIT(void)
{

#ifdef __NORMAL

#else
    timer::finish(ofs_running , !MY_RANK);

    Global_File::close_all_log(MY_RANK);

    if (MY_RANK==0)
    {
        Memory::print_all( ofs_running ) ;
    }
    cout<<" See output information in : "<<global_out_dir<<endl;
#endif

#ifdef __MPI
    MPI_Finalize();
#endif
    exit(0);
}


void WARNING_QUIT(const string &file,const string &description)
{
#ifdef __NORMAL

		cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "                         NOTICE                           " << endl;
		cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

#else

    //cout<<" ----------- SOMETHING TO WARN YOU ! -------" << endl;
	// 41: red background
	// 42: green background
	// 31: red
	// 32: greem
	// 33: yellow
	// 34: blue
	// 35: zi
	// 36: qing
	// 37: white
	if(COLOUR)
	{
		//printf( "\e[32m%s\e[0m\n", " -------------- SOMETHING TO WARN YOU ! ------------");
		//printf( " \e[32m%s\e[0m\n", description.c_str());
		//printf( " \e[32m%s\e[0m", "CHECK IN FILE : ");
		//printf( "\e[32m%s\e[0m", global_out_dir.c_str());
		//printf( "\e[32m%s\e[0m\n", "warning.log"); 
		//printf( "\e[32m%s\e[0m\n", " ---------------------------------------------------"); 
		printf( "[32m%s[0m\n", " -------------- SOMETHING TO WARN YOU ! ------------");
		printf( " [32m%s[0m\n", description.c_str());
		printf( " [32m%s[0m", "CHECK IN FILE : ");
		printf( "[32m%s[0m", global_out_dir.c_str());
		printf( "[32m%s[0m\n", "warning.log"); 
		printf( "[32m%s[0m\n", " ---------------------------------------------------"); 
	}
	else
	{
		cout << " " << endl;
		cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "                         NOTICE                           " << endl;
		cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << " " << endl;
		cout << " " << description << endl;
		cout << " CHECK IN FILE : " << global_out_dir << "warning.log" << endl;
		cout << " " << endl;
		cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "                         NOTICE                           " << endl;
		cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;


		ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		ofs_running << "                         NOTICE                           " << endl;
		ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		ofs_running << endl;
		ofs_running << " " << description << endl;
		ofs_running << " CHECK IN FILE : " << global_out_dir << "warning.log" << endl;
		ofs_running << endl;
		ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		ofs_running << "                         NOTICE                           " << endl;
		ofs_running << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	}

	WARNING(file,description);
    ofs_running<<" Check in file : "<<global_out_dir<<"warning.log"<<endl;

#endif

    QUIT();
}
