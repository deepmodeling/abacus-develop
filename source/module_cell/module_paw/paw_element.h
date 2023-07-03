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

//===================================================
// In paw_element.cpp : subroutines for reading paw xml file
//===================================================

    void read_paw_xml(std::string filename); //read info from paw file

    private:

    //some helper functions for reading the xml file
    //scan for line containing certain pattern from file
    std::string scan_file(std::ifstream &ifs, std::string pattern);

    //reset input buffer to the beginning
    void reset_buffer(std::ifstream &ifs);

    //extracting values from line
    //this part should be written with template; will try later
    std::string extract_string(std::string line, std::string key);
    double      extract_double(std::string line, std::string key);
    int         extract_int   (std::string line, std::string key);

    int count_nstates(std::ifstream &ifs); //count nstates
    void nstates_to_mstates(); //from nstates to mstates

    void get_nrcut(); //find grid point corresponding to certain radius r

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

    //for log grid, r_i = rstep * exp[(lstep * i)-1]
    //rstep <-> a, lstep <-> d from xml file
    double lstep, rstep;
    std::vector<double> rr, dr; //radial grid and increments
    std::vector<std::vector<double>> ptilde_r, ptilde_q; //projector functions in real and reciprocal space
    std::vector<std::vector<double>> d2ptilde_q; //2nd derivative; for spline
    std::vector<double> qgrid; //the q grid on which ptilde_q and d2ptilde_q are defined

//===================================================
// In paw_sphbes.cpp : subroutines for carrying out spherical bessel
// transformation of the projector functions
//===================================================

    public:

    // some of the functions should be put in a math lib
    // but later

    // spherical bessel function and derivative,
    // adapted from m_special_funcs/paw_jbessel_4spline and m_paw_numeric/paw_jbessel of ABINIT
    static void spherical_bessel_function(const int l, const double xx, 
        double & bes, double & besp, const bool calc_der);

    //Note as grid information and rcut is already in this class
    //I have chosen not to pass them around

    double spherical_bessel_transform(const int l, std::vector<double> & fr, const double q) const;

    //some helper functions for carrying out spherical bessel transformation
    //will switch to the one in math lib later

    // simpson integration
    double simpson_integration(std::vector<double> & f) const;

    // some preparation for simpson integration
    void prepare_simpson_integration(const double r_for_intg, int & meshsz, std::vector<double> & simp_fact) const;

    // helper function : spline and splint, for cubic interpolation
    // adapted from m_paw_numeric/paw_spline of ABINIT
    void spline(const std::vector<double> & r, const std::vector<double> & f,
        std::vector<double> & d2f, const double yp1, const double ypn) const;

    // f is the function stored on grid rr, x is the value at which to evaluate
    // adapted from m_splines/splfit from ABINIT
    double splint(const std::vector<double> & r, const std::vector<double> & f,
        const std::vector<double> & d2f, const double x) const;

    private:

    // converts projectors from real to reciprocal space
    // ptilde_l(q) = int_0^{rc} dr r^2 ptilde_l(r) j_l(qr)
    // this is doing the same work as calculating ffspl in libpaw
    void transform_ptilde();

//===================================================
// In paw_interface.cpp : communication with other classes
//===================================================

    public:

    // ecutwfc_in : unit in Rydberg
    void init_paw_element(const double ecutwfc_in, const double cell_factor_in);

    int get_mstates(){return mstates;}
    int get_nstates(){return nstates;}
    std::vector<int> get_lstate(){return lstate;}
    std::vector<int> get_mstate(){return mstate;}
    std::vector<int> get_im_to_istate(){return im_to_istate;}

    // return ptilde_q for a given channel at a given q_in, using spline
    double get_ptilde(const int istate_in, const double q_in);

    private:

    double ecutwfc;
    double cell_factor;
};

#endif