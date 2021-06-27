//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================
#include "driver.h"
#include "src_parallel/parallel_global.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"

#include <ctime>

void calculate();

int main(int argc, char **argv)
{
    Parallel_Global::read_mpi_parameters(argc,argv);

    calculate();

    return 0;
}


void calculate()
{

	time_t time_start = std::time(NULL);

	timer::start();

	//----------------------------------------------------------
	// main program for doing electronic structure calculations
	//----------------------------------------------------------
	Driver DD;
	DD.init();

	time_t	time_finish= std::time(NULL);

	// print out information before ABACUS ends
	cout << "\n START  Time  : " << ctime(&time_start);
	cout << " FINISH Time  : " << ctime(&time_finish);
	cout << " TOTAL  Time  : " << difftime(time_finish, time_start) << endl;
	cout << " SEE INFORMATION IN : "<<global_out_dir<<endl;

	ofs_running << "\n Start  Time  : " << ctime(&time_start);
	ofs_running << " Finish Time  : " << ctime(&time_finish);

	double total_time = difftime(time_finish, time_start);
	int hour = total_time / 3600;
	int mins = ( total_time - 3600 * hour ) / 60;
	int secs = total_time - 3600 * hour - 60 * mins ;
	ofs_running << " Total  Time  : " << hour << " h "
	            << mins << " mins "
	            << secs << " secs "<< endl;

#ifdef __MPI
    MPI_Finalize();
#endif

    return;
}
