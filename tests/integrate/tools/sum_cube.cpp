#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
    string input_file = argv[1];

    ifstream inp(input_file.c_str(), ios::in);

    if (!inp)
    {
        cout << "Can't find " << input_file << " !" << endl;
        return 1;
    }

    // skip the first two lines
    string tmpstring;
    for (int i = 0; i < 2; i++)
    {
        getline(inp, tmpstring);
    }

    // read the 3rd line: number of atoms + origin coordinates
    int natom;
    double origin_x, origin_y, origin_z;
    inp >> natom >> origin_x >> origin_y >> origin_z;
    getline(inp, tmpstring);

    // read the grid vectors (support non-orthogonal)
    double v1[3], v2[3], v3[3];
    int nx, ny, nz;
    inp >> nx >> v1[0] >> v1[1] >> v1[2];
    inp >> ny >> v2[0] >> v2[1] >> v2[2];
    inp >> nz >> v3[0] >> v3[1] >> v3[2];

    // calculate the volume element |v1 · (v2 × v3)|
    double volume = fabs(v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) - v1[1] * (v2[0] * v3[2] - v2[2] * v3[0])
                         + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]));

    getline(inp, tmpstring);

    // skip the atom coordinates
    for (int i = 0; i < natom; ++i)
    {
        getline(inp, tmpstring);
    }

    int nr = nx * ny * nz;

    double sum = 0.0;
    double env = 0.0;
    for (int i = 0; i < nr; i++)
    {
        inp >> env;
        sum += env;
    }

    double ne = 0.0;

    ne = sum * volume;
    std::cout << setprecision(10) << ne << std::endl;

    return 0;
}
