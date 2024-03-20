#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>


void calculate();

int main(int argc, char **argv)
{

	std::cout << "Hello, this is the 'base' module of ABACUS." << std::endl;

	std::cout << "The module searchs for the neighboring atoms for a given atomic position" << std::endl;

	std::cout << "Right now, the module is still empty, soon we will have more tests." << std::endl;

    calculate();

    return 0;
}


void calculate()
{
	std::ofstream ofs("log.txt");

//	ooo.set_orb_tables();

	ofs.close();

	std::cout << "--------------------" << std::endl;
	std::cout << " Have a great day! " << std::endl;
	std::cout << "--------------------" << std::endl;

/*
	time_t time_start = std::time(NULL);

//	ModuleBase::timer::start();

	//----------------------------------------------------------
	// main program for doing electronic structure calculations
	//----------------------------------------------------------
//	Driver DD;
//	DD.init();

	time_t	time_finish= std::time(NULL);

	// print out information before ABACUS ends
	std::cout << "\n START  Time  : " << ctime(&time_start);
	std::cout << " FINISH Time  : " << ctime(&time_finish);
	std::cout << " TOTAL  Time  : " << difftime(time_finish, time_start) << std::endl;

	double total_time = difftime(time_finish, time_start);
	int hour = total_time / 3600;
	int mins = ( total_time - 3600 * hour ) / 60;
	int secs = total_time - 3600 * hour - 60 * mins ;
	std::cout << " Total  Time  : " << hour << " h "
	            << mins << " mins "
	            << secs << " secs "<< std::endl;
*/

    return;
}
