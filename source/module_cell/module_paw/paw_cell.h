// The Paw_Cell class stores PAW-related information in unitcell
// including number of atoms, number of element types, and so on
// it basically serves as an interface between the unitcell class
// and the libpaw library from ABINIT

#ifndef PAW_CELL
#define PAW_CELL

#include <vector>
#include <complex>
#include <string>

#include "paw_element.h"

class Paw_Cell
{
    public:

    Paw_Cell(){};
    ~Paw_Cell(){};

    void init_paw_cell(
        const double ecutwfc_in, const double cell_factor_in,
        const int nat_in, const int ntyp_in,
        const int * atom_type_in, const double ** atom_coord_in,
        const std::vector<std::string> & filename_list_in,
        const int nx_in, const int ny_in, const int nz_in,
        const std::complex<double> * eigts1_in, const std::complex<double> * eigts2_in, const std::complex<double> * eigts3_in);

    // Given a list of k points, calculate the structure factors
    // exp(-i(k+G)R_I) = exp(-ikR_I) exp(-iG_xR_Ix) exp(-iG_yR_Iy) exp(-iG_zR_Iz)
    // as well as the spherical harmonics Ylm(k+G)
    void set_paw_k(
        const int npw, const double * kpt,
        const int * ig_to_ix, const int * ig_to_iy, const int * ig_to_iz,
        const double ** kpg);

    int get_nproj_tot(){return nproj_tot;}
    // map projector to atom
    std::vector<int> get_iprj_to_ia(){return iprj_to_ia;}
    // map projector to mstate of that element
    std::vector<int> get_iprj_to_im(){return iprj_to_im;}
    // map projector to lstate of that element
    std::vector<int> get_iprj_to_il(){return iprj_to_il;}
    // map projector to l quantum number of that element
    std::vector<int> get_iprj_to_l() {return iprj_to_l;}
    // map projector to m quantum number of that element
    std::vector<int> get_iprj_to_m() {return iprj_to_m;}
    // max l quantum number of all elements
    int get_lmax(){return lmax;}

    // ylm(r), adapted from initylmg of ABINIT
    static std::vector<double> calc_ylm(const int lmax, const double * r);

    // helper function for calc_ylm: Legendre polynomial
    static double ass_leg_pol(const int l, const int m, const double arg);

    private:

    // based on list of atom type, calculate total number of projectors
    // of this system, then map each projector to atom index, mstate and lstate
    // record in iproj_to_ia/im/il
    void map_paw_proj();

    // array of paw_element
    std::vector<Paw_Element> paw_element_list;

    int nproj_tot; // total number of projectors
    std::vector<int> iprj_to_ia; // map projector to atom
    std::vector<int> iprj_to_im; // map projector to mstate of that element
    std::vector<int> iprj_to_il; // map projector to lstate of that element
    std::vector<int> iprj_to_l;  // map projector to l quantum number of that element
    std::vector<int> iprj_to_m;  // map projector to m quantum number of that element

    // max l quantum number of all elements
    int lmax;

    // atomic positions and types
    int nat;
    int ntyp;
    std::vector<int> atom_type; // the element type of each atom
    std::vector<std::vector<double>> atom_coord; // Cartesian coordinate of each atom (in Bohr)

    // FFT grid
    int nx, ny, nz;

    // structure factor ('eigts1-3' from structure_factor class)
    // stores exp(- i G R_I) where G = (Gx,0,0), (0,Gy,0) and (0,0,Gz)
    std::vector<std::vector<std::complex<double>>> eigts1;
    std::vector<std::vector<std::complex<double>>> eigts2;
    std::vector<std::vector<std::complex<double>>> eigts3;

    // structure factor of (k+G) for current k point
    std::vector<std::vector<std::complex<double>>> struc_fact;

    // spherical harmonics Y_lm (k+G) for current k point
    std::vector<std::vector<double>> ylm_k;

    void set_ylm(const int npw, const double ** kpg);

};

#endif