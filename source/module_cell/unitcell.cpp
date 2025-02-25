#include <cstdlib>
#include <cstring> // Peize Lin fix bug about strcmp 2016-08-02

#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "unitcell.h"
#include "bcast_cell.h"
#include "module_parameter/parameter.h"
#include "module_cell/read_stru.h"
#include "module_base/atom_in.h"
#include "module_base/element_elec_config.h"
#include "module_base/global_file.h"
#include "module_base/parallel_common.h"
#include "module_parameter/parameter.h"

#ifdef __MPI
#include "mpi.h"
#endif
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif
#ifdef __LCAO
#include "../module_basis/module_ao/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif

#include "update_cell.h"
UnitCell::UnitCell() {
    itia2iat.create(1, 1);
}

UnitCell::~UnitCell() {
    if (set_atom_flag) {
        delete[] atoms;
    }
}


void UnitCell::print_cell(std::ofstream& ofs) const {

    ModuleBase::GlobalFunc::OUT(ofs, "print_unitcell()");

    ModuleBase::GlobalFunc::OUT(ofs, "latName", latName);
    ModuleBase::GlobalFunc::OUT(ofs, "ntype", ntype);
    ModuleBase::GlobalFunc::OUT(ofs, "nat", nat);
    ModuleBase::GlobalFunc::OUT(ofs, "lat0", lat0);
    ModuleBase::GlobalFunc::OUT(ofs, "lat0_angstrom", lat0_angstrom);
    ModuleBase::GlobalFunc::OUT(ofs, "tpiba", tpiba);
    ModuleBase::GlobalFunc::OUT(ofs, "omega", omega);

    output::printM3(ofs, "Lattices Vector (R) : ", latvec);
    output::printM3(ofs, "Supercell lattice vector : ", latvec_supercell);
    output::printM3(ofs, "Reciprocal lattice Vector (G): ", G);
    output::printM3(ofs, "GGT : ", GGT);

    ofs << std::endl;
    return;
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

std::vector<std::string> UnitCell::get_atomLabels() const {
    std::vector<std::string> atomLabels(this->ntype);
    for (int it = 0; it < this->ntype; it++) {
        atomLabels[it] = this->atoms[it].label;
    }
    return atomLabels;
}

std::vector<int> UnitCell::get_atomCounts() const {
    std::vector<int> atomCounts(this->ntype);
    for (int it = 0; it < this->ntype; it++) {
        atomCounts[it] = this->atoms[it].na;
    }
    return atomCounts;
}

std::vector<std::vector<int>> UnitCell::get_lnchiCounts() const {
    std::vector<std::vector<int>> lnchiCounts(this->ntype);
    for (int it = 0; it < this->ntype; it++) {
        lnchiCounts[it].resize(this->atoms[it].nwl + 1);
        for (int L = 0; L < this->atoms[it].nwl + 1; L++) {
            lnchiCounts[it][L] = this->atoms[it].l_nchi[L];
        }
    }
    return lnchiCounts;
}

std::vector<ModuleBase::Vector3<double>> UnitCell::get_target_mag() const
{
	std::vector<ModuleBase::Vector3<double>> target_mag(this->nat);
	for (int it = 0; it < this->ntype; it++)
	{
		for (int ia = 0; ia < this->atoms[it].na; ia++)
		{
			int iat = itia2iat(it, ia);
			target_mag[iat] = this->atoms[it].m_loc_[ia];
		}
	}
	return target_mag;
}

std::vector<ModuleBase::Vector3<double>> UnitCell::get_lambda() const
{
	std::vector<ModuleBase::Vector3<double>> lambda(this->nat);
	for (int it = 0; it < this->ntype; it++)
	{
		for (int ia = 0; ia < this->atoms[it].na; ia++)
		{
			int iat = itia2iat(it, ia);
			lambda[iat] = this->atoms[it].lambda[ia];
		}
	}
	return lambda;
}

std::vector<ModuleBase::Vector3<int>> UnitCell::get_constrain() const
{
	std::vector<ModuleBase::Vector3<int>> constrain(this->nat);
	for (int it = 0; it < this->ntype; it++)
	{
		for (int ia = 0; ia < this->atoms[it].na; ia++)
		{
			int iat = itia2iat(it, ia);
			constrain[iat] = this->atoms[it].constrain[ia];
		}
	}
	return constrain;
}

//==============================================================
// Calculate various lattice related quantities for given latvec
//==============================================================
void UnitCell::setup_cell(const std::string& fn, std::ofstream& log) {
    ModuleBase::TITLE("UnitCell", "setup_cell");
    // (1) init mag
    assert(ntype > 0);
    delete[] magnet.start_magnetization;
    magnet.start_magnetization = new double[this->ntype];

    // (2) init *Atom class array.
    this->atoms = new Atom[this->ntype]; // atom species.
    this->set_atom_flag = true;

    this->symm.epsilon = PARAM.inp.symmetry_prec;
    this->symm.epsilon_input = PARAM.inp.symmetry_prec;

    bool ok = true;
    bool ok2 = true;

    // (3) read in atom information
    this->atom_mass.resize(ntype);
    this->atom_label.resize(ntype);
    this->pseudo_fn.resize(ntype);
    this->pseudo_type.resize(ntype);
    this->orbital_fn.resize(ntype);
    if (GlobalV::MY_RANK == 0) {
        // open "atom_unitcell" file.
        std::ifstream ifa(fn.c_str(), std::ios::in);
        if (!ifa) 
        {
            GlobalV::ofs_warning << fn;
            ok = false;
        }

        if (ok) 
        {
            log << "\n\n\n\n";
            log << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                   ">>>>>>>>>>>>"
                << std::endl;
            log << " |                                                         "
                   "           |"
                << std::endl;
            log << " | Reading atom information in unitcell:                   "
                   "           |"
                << std::endl;
            log << " | From the input file and the structure file we know the "
                   "number of   |"
                << std::endl;
            log << " | different elments in this unitcell, then we list the "
                   "detail        |"
                << std::endl;
            log << " | information for each element, especially the zeta and "
                   "polar atomic |"
                << std::endl;
            log << " | orbital number for each element. The total atom number "
                   "is counted. |"
                << std::endl;
            log << " | We calculate the nearest atom distance for each atom "
                   "and show the  |"
                << std::endl;
            log << " | Cartesian and Direct coordinates for each atom. We list "
                   "the file   |"
                << std::endl;
            log << " | address for atomic orbitals. The volume and the lattice "
                   "vectors    |"
                << std::endl;
            log << " | in real and reciprocal space is also shown.             "
                   "           |"
                << std::endl;
            log << " |                                                         "
                   "           |"
                << std::endl;
            log << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                   "<<<<<<<<<<<<"
                << std::endl;
            log << "\n\n\n\n";

            log << " READING UNITCELL INFORMATION" << std::endl;
            //========================
            // call read_atom_species
            //========================
            const bool read_atom_species = unitcell::read_atom_species(ifa, log ,*this);
            //========================
            // call read_lattice_constant
            //========================
            const bool read_lattice_constant = unitcell::read_lattice_constant(ifa, log ,this->lat);
            //==========================
            // call read_atom_positions
            //==========================
            ok2 = unitcell::read_atom_positions(*this,ifa, log, GlobalV::ofs_warning);
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_bool(ok);
    Parallel_Common::bcast_bool(ok2);
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

#ifdef __MPI
    unitcell::bcast_unitcell(*this);
#endif

    //========================================================
    // Calculate unit cell volume
    // the reason to calculate volume here is
    // Firstly, latvec must be read in.
    //========================================================
    assert(lat0 > 0.0);
    this->omega = latvec.Det() * this->lat0 * lat0 * lat0;
    if (this->omega < 0)
    {
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        std::cout << " Warning: The lattice vector is left-handed; a right-handed vector is prefered." << std::endl;
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        GlobalV::ofs_warning << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        GlobalV::ofs_warning << " Warning: The lattice vector is left-handed; a right-handed vector is prefered." << std::endl;
        GlobalV::ofs_warning << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        this->omega = std::abs(this->omega);
    }
    else if (this->omega == 0)
    {
        ModuleBase::WARNING_QUIT("setup_cell", "The volume is zero.");
    }
    else
    {
        log << std::endl;
        ModuleBase::GlobalFunc::OUT(log, "Volume (Bohr^3)", this->omega);
        ModuleBase::GlobalFunc::OUT(log, "Volume (A^3)", this->omega * pow(ModuleBase::BOHR_TO_A, 3));
    }

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec have the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    this->GT = latvec.Inverse();
    this->G = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    // LiuXh add 20180515
    this->GT0 = latvec.Inverse();
    this->G0 = GT.Transpose();
    this->GGT0 = G * GT;
    this->invGGT0 = GGT.Inverse();

    log << std::endl;
    output::printM3(log,
                    "Lattice vectors: (Cartesian coordinate: in unit of a_0)",
                    latvec);
    output::printM3(
        log,
        "Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",
        G);
    //	OUT(log,"lattice center x",latcenter.x);
    //	OUT(log,"lattice center y",latcenter.y);
    //	OUT(log,"lattice center z",latcenter.z);

    //===================================
    // set index for iat2it, iat2ia
    //===================================
    this->set_iat2itia();

#ifdef USE_PAW
    if (PARAM.inp.use_paw) {
        GlobalC::paw_cell.set_libpaw_cell(latvec, lat0);

        int* typat;
        double* xred;

        typat = new int[nat];
        xred = new double[nat * 3];

        int iat = 0;
        for (int it = 0; it < ntype; it++) {
            for (int ia = 0; ia < atoms[it].na; ia++) {
                typat[iat] = it + 1; // Fortran index starts from 1 !!!!
                xred[iat * 3 + 0] = atoms[it].taud[ia].x;
                xred[iat * 3 + 1] = atoms[it].taud[ia].y;
                xred[iat * 3 + 2] = atoms[it].taud[ia].z;
                iat++;
            }
        }

        GlobalC::paw_cell.set_libpaw_atom(nat, ntype, typat, xred);
        delete[] typat;
        delete[] xred;

        GlobalC::paw_cell.set_libpaw_files();

        GlobalC::paw_cell.set_nspin(PARAM.inp.nspin);
    }
#endif

    return;
}


void UnitCell::set_iat2iwt(const int& npol_in) {
#ifdef __DEBUG
    assert(npol_in == 1 || npol_in == 2);
    assert(this->nat > 0);
    assert(this->ntype > 0);
#endif
    this->iat2iwt.resize(this->nat);
    this->npol = npol_in;
    int iat = 0;
    int iwt = 0;
    for (int it = 0; it < this->ntype; it++) {
        for (int ia = 0; ia < atoms[it].na; ia++) {
            this->iat2iwt[iat] = iwt;
            iwt += atoms[it].nw * this->npol;
            ++iat;
        }
    }
    return;
}



// check if any atom can be moved
bool UnitCell::if_atoms_can_move() const {
    for (int it = 0; it < this->ntype; it++) {
        Atom* atom = &atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            if (atom->mbl[ia].x || atom->mbl[ia].y || atom->mbl[ia].z) {
                return true;
}
        }
    }
    return false;
}

// check if lattice vector can be changed
bool UnitCell::if_cell_can_change() const {
    // need to be fixed next
    if (this->lc[0] || this->lc[1] || this->lc[2]) {
        return true;
    }
    return false;
}

void UnitCell::setup(const std::string& latname_in,
                     const int& ntype_in,
                     const int& lmaxmax_in,
                     const bool& init_vel_in,
                     const std::string& fixed_axes_in) {
    this->latName = latname_in;
    this->ntype = ntype_in;
    this->lmaxmax = lmaxmax_in;
    this->init_vel = init_vel_in;
    // pengfei Li add 2018-11-11
    if (fixed_axes_in == "None") {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "volume") {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
        if (!PARAM.inp.relax_new) {
            ModuleBase::WARNING_QUIT(
                "Input",
                "there are bugs in the old implementation; set relax_new to be "
                "1 for fixed_volume relaxation");
        }
    } else if (fixed_axes_in == "shape") {
        if (!PARAM.inp.relax_new) {
            ModuleBase::WARNING_QUIT(
                "Input",
                "set relax_new to be 1 for fixed_shape relaxation");
        }
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "a") {
        this->lc[0] = 0;
        this->lc[1] = 1;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "b") {
        this->lc[0] = 1;
        this->lc[1] = 0;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "c") {
        this->lc[0] = 1;
        this->lc[1] = 1;
        this->lc[2] = 0;
    } else if (fixed_axes_in == "ab") {
        this->lc[0] = 0;
        this->lc[1] = 0;
        this->lc[2] = 1;
    } else if (fixed_axes_in == "ac") {
        this->lc[0] = 0;
        this->lc[1] = 1;
        this->lc[2] = 0;
    } else if (fixed_axes_in == "bc") {
        this->lc[0] = 1;
        this->lc[1] = 0;
        this->lc[2] = 0;
    } else if (fixed_axes_in == "abc") {
        this->lc[0] = 0;
        this->lc[1] = 0;
        this->lc[2] = 0;
    } else {
        ModuleBase::WARNING_QUIT(
            "Input",
            "fixed_axes should be None,volume,shape,a,b,c,ab,ac,bc or abc!");
    }
    return;
}


void UnitCell::compare_atom_labels(std::string label1, std::string label2) {
    if (label1!= label2) //'!( "Ag" == "Ag" || "47" == "47" || "Silver" == Silver" )'
    {
        atom_in ai;
        if (!(std::to_string(ai.atom_Z[label1]) == label2
              ||                                  // '!( "Ag" == "47" )'
              ai.atom_symbol[label1] == label2 || // '!( "Ag" == "Silver" )'
              label1 == std::to_string(ai.atom_Z[label2])
              || // '!( "47" == "Ag" )'
              label1 == std::to_string(ai.symbol_Z[label2])
              ||                                  // '!( "47" == "Silver" )'
              label1 == ai.atom_symbol[label2] || // '!( "Silver" == "Ag" )'
              std::to_string(ai.symbol_Z[label1])
                  == label2)) // '!( "Silver" == "47" )'
        {
            std::string stru_label = "";
            std::string psuedo_label = "";
            for (int ip = 0; ip < label1.length(); ip++) {
                if (!(isdigit(label1[ip]) || label1[ip] == '_')) {
                    stru_label += label1[ip];
                } else {
                    break;
                }
            }
            stru_label[0] = toupper(stru_label[0]);

            for (int ip = 0; ip < label2.length(); ip++) {
                if (!(isdigit(label2[ip]) || label2[ip] == '_')) {
                    psuedo_label += label2[ip];
                } else {
                    break;
                }
            }
            psuedo_label[0] = toupper(psuedo_label[0]);

            if (!(stru_label == psuedo_label
                  || //' !("Ag1" == "ag_locpsp" || "47" == "47" || "Silver" ==
                     //Silver" )'
                  std::to_string(ai.atom_Z[stru_label]) == psuedo_label
                  || // ' !("Ag1" == "47" )'
                  ai.atom_symbol[stru_label] == psuedo_label
                  || // ' !("Ag1" == "Silver")'
                  stru_label == std::to_string(ai.atom_Z[psuedo_label])
                  || // ' !("47" == "Ag1" )'
                  stru_label == std::to_string(ai.symbol_Z[psuedo_label])
                  || // ' !("47" == "Silver1" )'
                  stru_label == ai.atom_symbol[psuedo_label]
                  || // ' !("Silver1" == "Ag" )'
                  std::to_string(ai.symbol_Z[stru_label])
                      == psuedo_label)) // ' !("Silver1" == "47" )'

            {
                std::string atom_label_in_orbtial
                    = "atom label in orbital file ";
                std::string mismatch_with_pseudo
                    = " mismatch with pseudo file of ";
                ModuleBase::WARNING_QUIT("UnitCell::read_pseudo",
                                         atom_label_in_orbtial + label1
                                             + mismatch_with_pseudo + label2);
            }
        }
    }
}