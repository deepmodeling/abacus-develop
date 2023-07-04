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
        const double * eigts1_in, const double * eigts2_in, const double * eigts3_in);

    void set_paw_k(
        const int npw, const double * kpt,
        const int * ig_to_ix, const int * ig_to_iy, const int * ig_to_iz);

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

    // structure factor of the current k points
    std::vector<std::vector<std::complex<double>>> struc_fact;
};

#endif