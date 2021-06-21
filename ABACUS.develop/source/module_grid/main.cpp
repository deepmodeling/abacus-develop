//#include "timer.h"
#include <ctime>

void calculate();

int main(int argc, char **argv)
{

    calculate();

    return 0;
}


void calculate()
{
/*
	time_t time_start = std::time(NULL);

//	timer::start();

	//----------------------------------------------------------
	// main program for doing electronic structure calculations
	//----------------------------------------------------------
//	Driver DD;
//	DD.init();

	time_t	time_finish= std::time(NULL);

	// print out information before ABACUS ends
	cout << "\n START  Time  : " << ctime(&time_start);
	cout << " FINISH Time  : " << ctime(&time_finish);
	cout << " TOTAL  Time  : " << difftime(time_finish, time_start) << endl;

	double total_time = difftime(time_finish, time_start);
	int hour = total_time / 3600;
	int mins = ( total_time - 3600 * hour ) / 60;
	int secs = total_time - 3600 * hour - 60 * mins ;
	cout << " Total  Time  : " << hour << " h "
	            << mins << " mins "
	            << secs << " secs "<< endl;
*/

    return;
}
