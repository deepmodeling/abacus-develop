#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "symmetry.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h"

void calculate();

int main(int argc, char **argv)
{

	std::cout << "Hello, this is the 'symmetry' module of ABACUS." << std::endl;

	std::cout << "The module does symmetry analysis for an input geometry." << std::endl;
	Parallel_Global::read_mpi_parameters(argc,argv);
	//std::cout << "Right now, the module is still empty, soon we will have more tests." << std::endl;

    calculate();

    return 0;
}


void calculate()
{
	//std::ofstream ofs("log.txt");
	std::ofstream ofs_running("log.txt");
	//std::ofstream ofs("useless.txt");
	std::ofstream ofs_warning("warning.txt");
	std::ifstream ifs("INPUT");
	UnitCell ucell;
	ModuleSymmetry::Symmetry symm;
	ifs >> ucell.ntype;
	ucell.latName = "none";
	ifs.close();
	output out;
	ucell.setup_cell(
	"STRU", 
	ofs_running);
	std::cout << "set up cell classic done." << std::endl;
	symm.analy_sys(ucell, out, ofs_running);
	ofs_running.close();
//	ooo.set_orb_tables();

	//ofs.close();

	std::cout << "--------------------" << std::endl;
	std::cout << " Have a great day! " << std::endl;
	std::cout << "--------------------" << std::endl;


    return;
}
