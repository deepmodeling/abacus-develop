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

	int nx = 0;
	int ny = 0;
	int nz = 0;
	int nr = 0;
	double a0 = 0.0;
	std::string tmpstring;
	getline(inp, tmpstring);
	inp >> a0;
	//std::cout<< "a0 "<< a0<<std::endl;
        for (int i = 0; i < 9; i++) getline(inp, tmpstring);
	inp >> nx >> ny >> nz;

	nr = nx*ny*nz;

	double sum = 0.0;
	double env = 0.0;
	for (int i=0; i< nr; i++)
	{
		inp >> env;
		sum += env*env;
	}

	double ne = 0.0;

	a0 /= 0.529177;
	ne = sum * a0 * a0 * a0 / nr;
	//std::cout<<nx<<"--"<<ny<<"--"<<nz<<std::endl;
	//std::cout<<" sum = "<<sum<< " ne = "<<ne<<std::endl;
	std::cout<<ne<<std::endl;

	
	return 0;
}
