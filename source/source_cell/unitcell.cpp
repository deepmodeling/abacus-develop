/**
 * @file unitcell.cpp
 * @brief Implementation of UnitCell class.
 */
#include <cstdlib>
#include <cstring> // Peize Lin fix bug about strcmp 2016-08-02

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "unitcell.h"
#include "bcast_cell.h"
#include "source_base/tool_quit.h"
#include "source_base/output.h"

#include "source_cell/read_stru.h"
#include "source_base/atom_in.h"
#include "source_base/element_elec_config.h"
#include "source_base/global_file.h"
#include "source_base/parallel_common.h"
#include "source_cell/sep_cell.h"

#ifdef __MPI
#include "mpi.h"
#endif

#include "update_cell.h"

UnitCell::UnitCell()
{
    itia2iat.create(1, 1);
}

UnitCell::~UnitCell()
{
    if (set_atom_flag)
    {
        delete[] atoms;
    }
}


void UnitCell::set_iat2itia() {
    assert(nat > 0);
    delete[] iat2it;
    delete[] iat2ia;
    this->iat2it = new int[nat];
    this->iat2ia = new int[nat];
    int iat = 0;
    for (int it = 0; it < ntype; it++) {
        for (int ia = 0; ia < atoms[it].na; ia++) {
            this->iat2it[iat] = it;
            this->iat2ia[iat] = ia;
            ++iat;
        }
    }
    return;
}

std::map<int, int> UnitCell::get_atom_Counts() const {
    std::map<int, int> atomCounts;
    for (int it = 0; it < this->ntype; it++) {
        atomCounts.insert(std::pair<int, int>(it, this->atoms[it].na));
    }
    return atomCounts;
}

std::map<int, int> UnitCell::get_orbital_Counts() const {
    std::map<int, int> orbitalCounts;
    for (int it = 0; it < this->ntype; it++) {
        orbitalCounts.insert(std::pair<int, int>(it, this->atoms[it].nw));
    }
    return orbitalCounts;
}

std::map<int, std::map<int, int>> UnitCell::get_lnchi_Counts() const {
    std::map<int, std::map<int, int>> lnchiCounts;
    for (int it = 0; it < this->ntype; it++) {
        for (int L = 0; L < this->atoms[it].nwl + 1; L++) {
            // Check if the key 'it' exists in the outer map
            if (lnchiCounts.find(it) == lnchiCounts.end()) {
                // If it doesn't exist, initialize an empty inner map
                lnchiCounts[it] = std::map<int, int>();
            }
            int l_nchi = this->atoms[it].l_nchi[L];
            // Insert the key-value pair into the inner map
            lnchiCounts[it].insert(std::pair<int, int>(L, l_nchi));
        }
    }
    return lnchiCounts;
}

//==============================================================
// Calculate various lattice related quantities for given latvec
//==============================================================
void UnitCell::setup_cell(const std::string& fn, std::ofstream& log, const double symmetry_prec, const int dfthalf_type, const std::string& pseudo_dir, const int nspin,
    const std::string& basis_type, const std::string& orbital_dir, const std::string& init_wfc,
    const double onsite_radius, const bool deepks_setorb, const bool rpa,
    const bool fixed_atoms, const bool noncolin, const std::string& calculation, const std::string& esolver_type,
    const int symmetry)
{
    ModuleBase::TITLE("UnitCell", "setup_cell");

    assert(ntype > 0);

    // (1) init *Atom class array.
    this->atoms = new Atom[this->ntype]; // atom species.
    this->set_atom_flag = true;

    this->symm.epsilon = symmetry_prec;
    this->symm.epsilon_input = symmetry_prec;

    bool ok = true;
    bool ok2 = true;

    bool ok3 = true; // for sep potential in DFT-1/2

    // (3) read in atom information
    this->atom_mass.resize(ntype);
    this->atom_label.resize(ntype);
    this->pseudo_fn.resize(ntype);
    this->pseudo_type.resize(ntype);
    this->orbital_fn.resize(ntype);

    if (GlobalV::MY_RANK == 0)
    {
        // open "atom_unitcell" file.
        std::ifstream ifa(fn.c_str(), std::ios::in);
        if (!ifa)
        {
            GlobalV::ofs_warning << fn;
            ok = false;
        }

        if (ok)
        {
            log << "\n\n";
            log << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            log << " |                                                                    |" << std::endl;
            log << " |                        #Setup Unitcell#                            |" << std::endl;
            log << " | From the input file and the structure file we know the number of   |" << std::endl;
            log << " | different elments in this unitcell, then we list the detail        |" << std::endl;
            log << " | information for each element, especially the zeta and polar atomic |" << std::endl;
            log << " | orbital number for each element. The total atom number is counted. |" << std::endl;
            log << " | We calculate the nearest atom distance for each atom and show the  |" << std::endl;
            log << " | Cartesian and Direct coordinates for each atom. We list the file   |" << std::endl;
            log << " | address for atomic orbitals. The volume and the lattice vectors    |" << std::endl;
            log << " | in real and reciprocal space is also shown.                        |" << std::endl;
            log << " |                                                                    |" << std::endl;
            log << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
            log << "\n";

            log << " READING UNITCELL INFORMATION" << std::endl;
            //========================
            // call read_atom_species
            //========================
            const bool read_atom_species = unitcell::read_atom_species(ifa, log, *this,
                basis_type, orbital_dir, init_wfc, onsite_radius, deepks_setorb, rpa);
            //========================
            // call read_lattice_constant
            //========================
            const bool read_lattice_constant = unitcell::read_lattice_constant(ifa, log ,this->lat);
            //==========================
            // readl sep potential, currently using the pseudopotential folder (pseudo_dir in INPUT)
            //==========================
            if (dfthalf_type > 0) {
                sep_cell.init(this->ntype);
                ok3 = sep_cell.read_sep_potentials(ifa, pseudo_dir, GlobalV::ofs_warning, this->atom_label);
            }
            //==========================
            // call read_atom_positions
            //==========================
            ok2 = unitcell::read_atom_positions(*this, ifa, log, GlobalV::ofs_warning, nspin,
                basis_type, orbital_dir, init_wfc, onsite_radius, fixed_atoms, noncolin,
                calculation, esolver_type, symmetry);
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_bool(ok);
    Parallel_Common::bcast_bool(ok2);
    Parallel_Common::bcast_bool(ok3);
#endif
    if (!ok) {
        ModuleBase::WARNING_QUIT(
            "UnitCell::setup_cell",
            "Can not find the file containing atom positions.!");
    }
    if (!ok2) {
        ModuleBase::WARNING_QUIT("UnitCell::setup_cell",
                                 "Something wrong during read_atom_positions.");
    }
    if (!ok3) {
        ModuleBase::WARNING_QUIT("UnitCell::setup_cell", "Something wrong during read_sep_potentials");
    }

#ifdef __MPI
    unitcell::bcast_unitcell(*this, nspin);
    sep_cell.bcast_sep_cell();
#endif

    //========================================================
    // Calculate unit cell volume
    // the reason to calculate volume here is
    // Firstly, latvec must be read in.
    //========================================================
    assert(lat0 > 0.0);
    this->omega = latvec.Det() * this->lat0 * this->lat0 * this->lat0;


    if (this->omega < 0)
    {
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        std::cout << " Warning: The lattice vector is left-handed; a right-handed vector is prefered." << std::endl;
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        GlobalV::ofs_warning <<
        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        GlobalV::ofs_warning <<
        " Warning: The lattice vector is left-handed; a right-handed vector is prefered." << std::endl;
        GlobalV::ofs_warning <<
        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        this->omega = std::abs(this->omega);
    }
    else if (this->omega == 0)
    {
        ModuleBase::WARNING_QUIT("setup_cell", "The volume is zero.");
    }
    else
    {
        ModuleBase::GlobalFunc::OUT(log, "Cell volume (Bohr^3)", this->omega);
        ModuleBase::GlobalFunc::OUT(log, "Cell volume (A^3)", this->omega * pow(ModuleBase::BOHR_TO_A, 3));
    }

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec have the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    this->GT = latvec.Inverse();
    this->G = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    log << std::endl;
    output::printM3(log,
                    "Lattice vectors: (Cartesian coordinate: in unit of a_0)",
                    latvec);
    output::printM3(
        log,
        "Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",
        G);

    //===================================
    // set index for iat2it, iat2ia
    //===================================
    this->set_iat2itia();

    sep_cell.set_omega(this->omega, this->tpiba2);

    return;
}


void UnitCell::set_iat2iwt(const int& npol_in)
{
#ifdef __DEBUG
    assert(npol_in == 1 || npol_in == 2);
    assert(this->nat > 0);
    assert(this->ntype > 0);
#endif
    this->iat2iwt.resize(this->nat);
    this->npol = npol_in;
    int iat = 0;
    int iwt = 0;

    for (int it = 0; it < this->ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            this->iat2iwt[iat] = iwt;
            iwt += atoms[it].nw * this->npol;
            ++iat;
        }
    }
    return;
}



void UnitCell::setup(const std::string& latname_in,
                     const int& ntype_in,
                     const int& lmaxmax_in,
                     const bool& init_vel_in,
                     const std::string& fixed_axes_in) {
    this->latName = latname_in;
    this->ntype = ntype_in;
    this->magnet.start_mag.resize(ntype_in, 0.0);
    this->lmaxmax = lmaxmax_in;
    this->init_vel = init_vel_in;
    // pengfei Li add 2018-11-11
    if (fixed_axes_in == "None") {
        this->lat_axis_free[0] = 1;
        this->lat_axis_free[1] = 1;
        this->lat_axis_free[2] = 1;
    } else if (fixed_axes_in == "volume") {
        this->lat_axis_free[0] = 1;
        this->lat_axis_free[1] = 1;
        this->lat_axis_free[2] = 1;
    } else if (fixed_axes_in == "shape") {
        this->lat_axis_free[0] = 1;
        this->lat_axis_free[1] = 1;
        this->lat_axis_free[2] = 1;
    } else if (fixed_axes_in == "a") {
        this->lat_axis_free[0] = 0;
        this->lat_axis_free[1] = 1;
        this->lat_axis_free[2] = 1;
    } else if (fixed_axes_in == "b") {
        this->lat_axis_free[0] = 1;
        this->lat_axis_free[1] = 0;
        this->lat_axis_free[2] = 1;
    } else if (fixed_axes_in == "c") {
        this->lat_axis_free[0] = 1;
        this->lat_axis_free[1] = 1;
        this->lat_axis_free[2] = 0;
    } else if (fixed_axes_in == "ab") {
        this->lat_axis_free[0] = 0;
        this->lat_axis_free[1] = 0;
        this->lat_axis_free[2] = 1;
    } else if (fixed_axes_in == "ac") {
        this->lat_axis_free[0] = 0;
        this->lat_axis_free[1] = 1;
        this->lat_axis_free[2] = 0;
    } else if (fixed_axes_in == "bc") {
        this->lat_axis_free[0] = 1;
        this->lat_axis_free[1] = 0;
        this->lat_axis_free[2] = 0;
    } else if (fixed_axes_in == "abc") {
        this->lat_axis_free[0] = 0;
        this->lat_axis_free[1] = 0;
        this->lat_axis_free[2] = 0;
    } else {
        ModuleBase::WARNING_QUIT(
            "Input",
            "fixed_axes should be none, volume, shape, a, b, c, ab, ac, bc or abc!");
    }
    return;
}
