#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

int main(int argc, char **argv)
{
	string input_file = argv[1];

	ifstream inp(input_file.c_str(), ios::in);

	if (!inp) 
	{
		cout << "Can't find " << input_file << " !" << endl;
		return 1;
	}
    else
    {
//        cout << "Find the " << input_file << endl;
    }

	int nx = 0;
	int ny = 0;
	int nz = 0;
	int nr = 0;
	double mx = 0.0;
	double my = 0.0;
	double mz = 0.0;
	double tmp = 0.0;
	std::string tmpstring;
	getline(inp, tmpstring);
	getline(inp, tmpstring);
	getline(inp, tmpstring);
	inp >> nx >> mx >> tmp >> tmp;
	inp >> ny >> tmp >> my >> tmp;
	inp >> nz >> tmp >> tmp >> mz;
	getline(inp, tmpstring);
	getline(inp, tmpstring);

	nr = nx * ny * nz;

	double sum = 0.0;
	double env = 0.0;
	for (int i=0; i< nr; i++)
	{
		inp >> env;
		sum += env;
	}

	double ne = 0.0;

	ne = sum * mx * my * mz;
	std::cout<<setprecision(10)<<ne<<std::endl;
	
	return 0;
}
