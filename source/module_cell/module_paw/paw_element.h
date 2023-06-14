// The Paw_Element class reads information from PAW xml input files
// and stores them in corresponding data structures
// Since most of the operations will be carried out using libpaw from ABINIT
// only part of information is necessary

#ifndef PAW_ELEMENT
#define PAW_ELEMENT

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Paw_Element
{
    public:

    Paw_Element(){};
    ~Paw_Element(){};

    void read_paw_xml(std::string filename); //read info from paw file
    void transform_ptilde(); //projectors from real to reciprocal space

    private:

    double Zat, core, val; //atomic number, core electron & valence electron
    std::string symbol; //element symbol

    double rcut; //radius of augmentation sphere
    int    nr, nrcut; //size of radial grid; radial grid point corresponding to rcut
    //note : nrcut will be the upper bound of later radial integrations
    
    int    nstates; //number of channels (quantum numbers n,l)
    std::vector<int> lstate; //l quantum number of each channel

    int    mstates; //#. m states (for each (n,l) channel, there will be 2l+1 m states)
    std::vector<int> mstate; //m quantum number of each mstate
    std::vector<int> im_to_istate; //map from mstate to (n,l) channel (namely nstates)

    std::vector<double> rr, dr; //radial grid and increments
    std::vector<std::vector<double>> ptilde_r; //projector functions in real space


    //some helper functions for reading the xml file
    //scan for line containing certain pattern from file
    std::string scan_file(std::ifstream &ifs, std::string pattern);

    //reset input buffer to the beginning
    void reset_buffer(std::ifstream &ifs);

    //extracting values from line
    //this part should be written with template; will try later
    std::string extract_string(std::string line, std::string key);
    double extract_double(std::string line, std::string key);
    int extract_int(std::string line, std::string key);

    int count_nstates(std::ifstream &ifs); //count nstates
    void nstates_to_mstates(); //from nstates to mstates

    //find grid point corresponding to rcut
    void get_nrcut();

};

#endif