//#include "timer.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

void calculate();

int main(int argc, char **argv)
{

	cout << "Hello, this is the ADJ module of ABACUS." << endl;

	cout << "The module searchs for the adjacent atoms for a given atomic position" << endl;

	cout << "Right now, the module is still empty, soon we will have more tests." << endl;

    calculate();

    return 0;
}


void calculate()
{
	ofstream ofs("log.txt");

//	ooo.set_orb_tables();

	ofs.close();

	cout << "--------------------" << endl;
	cout << " Have a great day! " << endl;
	cout << "--------------------" << endl;

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
