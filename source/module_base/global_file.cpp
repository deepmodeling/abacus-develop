//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================
#include "global_file.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include "global_function.h"
#include "global_variable.h"
#include "module_base/parallel_common.h"
#include "module_base/parallel_reduce.h"

//----------------------------------------------------------
// EXPLAIN : Be Called in input.cpp
//----------------------------------------------------------
namespace ModuleBase
{
void ModuleBase::Global_File::make_dir_out(
    const std::string &suffix,
	const std::string &calculation,
    const bool &out_dir,
    const int rank,
    const bool &restart,
    const bool out_alllog)
{
//----------------------------------------------------------
// USE STL FUNCTION
// NAME : system
//----------------------------------------------------------

    std::string prefix ;

#ifdef __EPM
#ifdef __MPI
    prefix = "OUT_EPM_MPI.";
#else
    prefix = "OUT_EPM.";
#endif
#else
    prefix = "OUT.";
#endif

    GlobalV::global_out_dir = prefix + suffix + "/";
    GlobalV::global_stru_dir = GlobalV::global_out_dir + "STRU/";
    GlobalV::global_matrix_dir = GlobalV::global_out_dir + "matrix/";

#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    int make_dir = 0;
	// mohan update 2011-05-03
    std::string command0 =  "test -d " + GlobalV::global_out_dir + " || mkdir " + GlobalV::global_out_dir;

	int times = 0;
	while(times<GlobalV::NPROC)
	{
		if(rank==times)
		{
			if ( system( command0.c_str() ) == 0 )
			{
				std::cout << " MAKE THE DIR         : " << GlobalV::global_out_dir << std::endl;
				make_dir = 1;
			}
			else
			{
				std::cout << " PROC " << rank << " CAN NOT MAKE THE DIR !!! " << std::endl;
				make_dir = 0;
			}
        }
#ifdef __MPI
		Parallel_Reduce::reduce_int_all(make_dir);
#endif
		if(make_dir>0)break;
		++times;
	}

#ifdef __MPI
	if(make_dir==0)
	{
		std::cout << " CAN NOT MAKE THE OUT DIR......." << std::endl;
		ModuleBase::QUIT();
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif

    if(calculation == "md")
    {
        int make_dir_stru = 0;
        std::string command1 =  "test -d " + GlobalV::global_stru_dir + " || mkdir " + GlobalV::global_stru_dir;

        times = 0;
        while(times<GlobalV::NPROC)
        {
            if(rank==times)
            {
                if ( system( command1.c_str() ) == 0 )
                {
                    std::cout << " MAKE THE STRU DIR    : " << GlobalV::global_stru_dir << std::endl;
                    make_dir_stru = 1;
                }
                else
                {
                    std::cout << " PROC " << rank << " CAN NOT MAKE THE STRU DIR !!! " << std::endl;
                    make_dir_stru = 0;
                }
            }
#ifdef __MPI
            Parallel_Reduce::reduce_int_all(make_dir_stru);
#endif
            if(make_dir_stru>0) break;
            ++times;
        }

#ifdef __MPI
        if(make_dir_stru==0)
        {
            std::cout << " CAN NOT MAKE THE STRU DIR......." << std::endl;
            ModuleBase::QUIT();
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    // make dir for HS matrix output in md calculation
    if((out_dir) && calculation == "md")
    {
        int make_dir_matrix = 0;
        std::string command1 =  "test -d " + GlobalV::global_matrix_dir + " || mkdir " + GlobalV::global_matrix_dir;

        times = 0;
        while(times<GlobalV::NPROC)
        {
            if(rank==times)
            {
                if ( system( command1.c_str() ) == 0 )
                {
                    std::cout << " MAKE THE MATRIX DIR    : " << GlobalV::global_stru_dir << std::endl;
                    make_dir_matrix = 1;
                }
                else
                {
                    std::cout << " PROC " << rank << " CAN NOT MAKE THE MATRIX DIR !!! " << std::endl;
                    make_dir_matrix = 0;
                }
            }
#ifdef __MPI
            Parallel_Reduce::reduce_int_all(make_dir_matrix);
#endif
            if(make_dir_matrix>0) break;
            ++times;
        }

#ifdef __MPI
        if(make_dir_matrix==0)
        {
            std::cout << " CAN NOT MAKE THE MATRIX DIR......." << std::endl;
            ModuleBase::QUIT();
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    std::stringstream ss,ss1;

    // mohan add 2010-09-12
    if(out_alllog)
    {
	    ss << "running_" << calculation << "_" << rank + 1;
	    open_log(GlobalV::ofs_running, ss.str(), calculation, restart);
        #if defined(__CUDA) || defined(__ROCM)
        open_log(GlobalV::ofs_device, "device" + std::to_string(rank), calculation, restart);
        #endif
    }
    else
    {
	    if(rank==0)
	    {
		    ss << "running_" << calculation;
		    open_log(GlobalV::ofs_running, ss.str(), calculation, restart);
            #if defined(__CUDA) || defined(__ROCM)
            open_log(GlobalV::ofs_device, "device", calculation, restart);
            #endif
	    }
    }

    if(rank==0)
    {
        open_log(GlobalV::ofs_warning, "warning", calculation, restart);
    }

#ifdef GATHER_INFO
    open_log(GlobalV::ofs_info, "math_info_" + std::to_string(rank), calculation, restart);
#endif

    return;
}

void ModuleBase::Global_File::make_dir_atom(const std::string &label)
{
//----------------------------------------------------------
// EXPLAIN : generate atom dir for each type of atom
//----------------------------------------------------------
    std::stringstream ss;
    ss << GlobalV::global_out_dir << label << "/";
    ModuleBase::GlobalFunc::MAKE_DIR(ss.str());
    return;
}

void ModuleBase::Global_File::open_log(std::ofstream &ofs, const std::string &fn, const std::string &calculation, const bool &restart)
{
//----------------------------------------------------------
// USE GLOBAL VARIABLE :
// GlobalV::global_out_dir : (default dir to store "*.log" file)
//----------------------------------------------------------
    std::stringstream ss;
    ss << GlobalV::global_out_dir << fn << ".log";

    if(calculation == "md" && restart)
    {
        ofs.open(ss.str(), std::ios::app);
    }
    else
    {
        ofs.open( ss.str() );
    }
//	ofs << " WELCOME TO MESIA PROGRAM." << std::endl;
//	ofs << " OPEN "<<fn<<".log"<<" DONE."<<std::endl;
    return;
}

void ModuleBase::Global_File::close_log( std::ofstream &ofs,const std::string &fn)
{
	if(ofs)
	{
    	ofs.close();
	}
    //ofs << "CLOSE "<<fn<<".log"<<" DONE."<<std::endl;
    return;
}

void ModuleBase::Global_File::close_all_log(const int rank, const bool out_alllog)
{
//----------------------------------------------------------
// USE GLOBAL VARIABLES :
// NAME : GlobalV::ofs_running
// NAME : GlobalV::ofs_warning
// NAME : ofs_recon
// NAME : ofs_sph_proj
// NAME : ofs_build
//----------------------------------------------------------

	// mohan update 2011-01-13
    std::stringstream ss;
	if(out_alllog)
	{
    	ss << "running_" << GlobalV::CALCULATION << "_cpu" << rank << ".log";
    	close_log(GlobalV::ofs_running,ss.str());
        #if defined(__CUDA) || defined(__ROCM)
        close_log(GlobalV::ofs_device, "device" + std::to_string(rank));
        #endif
	}
	else
	{
		if(rank==0)
		{
    		ss << "running_" << GlobalV::CALCULATION << ".log";
    		close_log(GlobalV::ofs_running,ss.str());
            #if defined(__CUDA) || defined(__ROCM)
            close_log(GlobalV::ofs_device, "device");
            #endif
		}
	}

    if (rank==0)
    {
        close_log(GlobalV::ofs_warning, "warning.log");
    }

#ifdef GATHER_INFO
    close_log(GlobalV::ofs_info, "math_info");
#endif
    return;
}
}
