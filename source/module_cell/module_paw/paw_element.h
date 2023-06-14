// The Paw_Element class reads information from PAW xml input files
// and stores them in corresponding data structures
// Since most of the operations will be carried out using libpaw from ABINIT
// only part of information is necessary

#ifndef PAW_ELEMENT
#define PAW_ELEMENT

#include <vector>
#include <string>
#include <fstream>

class Paw_Element
{
    public:

    Paw_Element(){};
    ~Paw_Element(){};

    void read_paw_xml(std::string filename); //read info from paw file

    private:

    std::string scan_file(std::ifstream &ifs, std::string pattern); //scan for line containing certain pattern from file

    double Zat, core, val; //atomic number, core electron & valence electron
    char   symbol[2]; //element symbol

    double rcut; //radius of augmentation sphere
    int    nr, nrcut; //size of radial grid; radial grid point corresponding to rcut
    int    nl, nstates; //number of different l quantum numbers; number of channels

    std::vector<int> lstate; //l quantum number of each state

    int    mstates; //#. m states (for each l state, there will be 2l+1 m states)
    std::vector<int> mstate; //m quantum number of each mstate
    std::vector<int> im_to_il; //map from mstate to l channel (namely nstates)

    std::vector<double> rr, dr; //radial grid and increments
    std::vector<std::vector<double>> ptilde_r; //projector functions in real and reciprocal space

};

#endif