#include <cstdlib>
#ifdef __MPI
#include "mpi.h"
#endif

#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "unitcell.h"
#include "module_parameter/parameter.h"
#include "cal_atoms_info.h"

#ifdef __LCAO
#include "../module_basis/module_ao/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif
#include "module_base/atom_in.h"
#include "module_base/element_elec_config.h"
#include "module_base/global_file.h"
#include "module_base/parallel_common.h"

#include <cstring> // Peize Lin fix bug about strcmp 2016-08-02
#include "module_parameter/parameter.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif
#ifdef __EXX
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_ri/serialization_cereal.h"
#endif

UnitCell::UnitCell() {
    if (test_unitcell) {
        ModuleBase::TITLE("unitcell", "Constructor");
}
    Coordinate = "Direct";
    latName = "none";
    lat0 = 0.0;
    lat0_angstrom = 0.0;

    ntype = 0;
    nat = 0;
    namax = 0;
    nwmax = 0;

    iat2it = nullptr;
    iat2ia = nullptr;
    iwt2iat = nullptr;
    iwt2iw = nullptr;

    itia2iat.create(1, 1);
    lc = new int[3];

    latvec = ModuleBase::Matrix3();
    latvec_supercell = ModuleBase::Matrix3();
    G = ModuleBase::Matrix3();
    GT = ModuleBase::Matrix3();
    GGT = ModuleBase::Matrix3();
    invGGT = ModuleBase::Matrix3();

    tpiba = 0.0;
    tpiba2 = 0.0;
    omega = 0.0;

    atom_label = new std::string[1];
    atom_mass = nullptr;
    pseudo_fn = new std::string[1];
    pseudo_type = new std::string[1];
    orbital_fn = new std::string[1];

    set_atom_flag = false;
}

UnitCell::~UnitCell() {
    delete[] atom_label;
    delete[] atom_mass;
    delete[] pseudo_fn;
    delete[] pseudo_type;
    delete[] orbital_fn;
    delete[] iat2it;
    delete[] iat2ia;
    delete[] iwt2iat;
    delete[] iwt2iw;
    delete[] lc;
    if (set_atom_flag) {
        delete[] atoms;
    }
}

#include "module_base/parallel_common.h"
#ifdef __MPI
void UnitCell::bcast_unitcell() {
    if (test_unitcell) {
        ModuleBase::TITLE("UnitCell", "bcast_unitcell");
}
    Parallel_Common::bcast_string(Coordinate);
    Parallel_Common::bcast_int(nat);

    Parallel_Common::bcast_double(lat0);
    Parallel_Common::bcast_double(lat0_angstrom);
    Parallel_Common::bcast_double(tpiba);
    Parallel_Common::bcast_double(tpiba2);

    // distribute lattice vectors.
    Parallel_Common::bcast_double(latvec.e11);
    Parallel_Common::bcast_double(latvec.e12);
    Parallel_Common::bcast_double(latvec.e13);
    Parallel_Common::bcast_double(latvec.e21);
    Parallel_Common::bcast_double(latvec.e22);
    Parallel_Common::bcast_double(latvec.e23);
    Parallel_Common::bcast_double(latvec.e31);
    Parallel_Common::bcast_double(latvec.e32);
    Parallel_Common::bcast_double(latvec.e33);

    Parallel_Common::bcast_int(lc[0]);
    Parallel_Common::bcast_int(lc[1]);
    Parallel_Common::bcast_int(lc[2]);

    // distribute lattice vectors.
    Parallel_Common::bcast_double(a1.x);
    Parallel_Common::bcast_double(a1.y);
    Parallel_Common::bcast_double(a1.z);
    Parallel_Common::bcast_double(a2.x);
    Parallel_Common::bcast_double(a2.y);
    Parallel_Common::bcast_double(a2.z);
    Parallel_Common::bcast_double(a3.x);
    Parallel_Common::bcast_double(a3.y);
    Parallel_Common::bcast_double(a3.z);

    // distribute latcenter
    Parallel_Common::bcast_double(latcenter.x);
    Parallel_Common::bcast_double(latcenter.y);
    Parallel_Common::bcast_double(latcenter.z);

    // distribute superlattice vectors.
    Parallel_Common::bcast_double(latvec_supercell.e11);
    Parallel_Common::bcast_double(latvec_supercell.e12);
    Parallel_Common::bcast_double(latvec_supercell.e13);
    Parallel_Common::bcast_double(latvec_supercell.e21);
    Parallel_Common::bcast_double(latvec_supercell.e22);
    Parallel_Common::bcast_double(latvec_supercell.e23);
    Parallel_Common::bcast_double(latvec_supercell.e31);
    Parallel_Common::bcast_double(latvec_supercell.e32);
    Parallel_Common::bcast_double(latvec_supercell.e33);
    Parallel_Common::bcast_double(magnet.start_magnetization, ntype);

    if (PARAM.inp.nspin == 4) {
        Parallel_Common::bcast_double(magnet.ux_[0]);
        Parallel_Common::bcast_double(magnet.ux_[1]);
        Parallel_Common::bcast_double(magnet.ux_[2]);
    }

    for (int i = 0; i < ntype; i++) {
        atoms[i].bcast_atom(); // init tau and mbl array
    }

#ifdef __EXX
    ModuleBase::bcast_data_cereal(GlobalC::exx_info.info_ri.files_abfs,
                                  MPI_COMM_WORLD,
                                  0);
#endif
    return;
}

void UnitCell::bcast_unitcell2() {
    for (int i = 0; i < ntype; i++) {
        atoms[i].bcast_atom2();
    }
    return;
}
#endif

void UnitCell::print_cell(std::ofstream& ofs) const {
    if (test_unitcell) {
        ModuleBase::TITLE("UnitCell", "print_cell");
}

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

/*
void UnitCell::print_cell_xyz(const std::string& fn) const
{
    if (test_unitcell)
        ModuleBase::TITLE("UnitCell", "print_cell_xyz");

    if (GlobalV::MY_RANK != 0)
        return; // xiaohui add 2015-03-15

    std::stringstream ss;
    ss << PARAM.globalv.global_out_dir << fn;

    std::ofstream ofs(ss.str().c_str());

    ofs << nat << std::endl;
    ofs << latName << std::endl;
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            ofs << atoms[it].label << " " << atoms[it].tau[ia].x * lat0 *
0.529177 << " "
                << atoms[it].tau[ia].y * lat0 * 0.529177 << " " <<
atoms[it].tau[ia].z * lat0 * 0.529177 << std::endl;
        }
    }

    ofs.close();
    return;
}
*/

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

void UnitCell::update_pos_tau(const double* pos) {
    int iat = 0;
    for (int it = 0; it < this->ntype; it++) {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            for (int ik = 0; ik < 3; ++ik) {
                if (atom->mbl[ia][ik]) {
                    atom->dis[ia][ik]
                        = pos[3 * iat + ik] / this->lat0 - atom->tau[ia][ik];
                    atom->tau[ia][ik] = pos[3 * iat + ik] / this->lat0;
                }
            }

            // the direct coordinates also need to be updated.
            atom->dis[ia] = atom->dis[ia] * this->GT;
            atom->taud[ia] = atom->tau[ia] * this->GT;
            iat++;
        }
    }
    assert(iat == this->nat);
    this->periodic_boundary_adjustment();
    this->bcast_atoms_tau();
}

void UnitCell::update_pos_taud(double* posd_in) {
    int iat = 0;
    for (int it = 0; it < this->ntype; it++) {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            for (int ik = 0; ik < 3; ++ik) {
                atom->taud[ia][ik] += posd_in[3 * iat + ik];
                atom->dis[ia][ik] = posd_in[3 * iat + ik];
            }
            iat++;
        }
    }
    assert(iat == this->nat);
    this->periodic_boundary_adjustment();
    this->bcast_atoms_tau();
}

// posd_in is atomic displacements here  liuyu 2023-03-22
void UnitCell::update_pos_taud(const ModuleBase::Vector3<double>* posd_in) {
    int iat = 0;
    for (int it = 0; it < this->ntype; it++) {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            for (int ik = 0; ik < 3; ++ik) {
                atom->taud[ia][ik] += posd_in[iat][ik];
                atom->dis[ia][ik] = posd_in[iat][ik];
            }
            iat++;
        }
    }
    assert(iat == this->nat);
    this->periodic_boundary_adjustment();
    this->bcast_atoms_tau();
}

void UnitCell::update_vel(const ModuleBase::Vector3<double>* vel_in) {
    int iat = 0;
    for (int it = 0; it < this->ntype; ++it) {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ++ia) {
            this->atoms[it].vel[ia] = vel_in[iat];
            ++iat;
        }
    }
    assert(iat == this->nat);
}

void UnitCell::periodic_boundary_adjustment() {
    //----------------------------------------------
    // because of the periodic boundary condition
    // we need to adjust the atom positions,
    // first adjust direct coordinates,
    // then update them into cartesian coordinates,
    //----------------------------------------------
    for (int it = 0; it < this->ntype; it++) {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            // mohan update 2011-03-21
            if (atom->taud[ia].x < 0) {
                atom->taud[ia].x += 1.0;
}
            if (atom->taud[ia].y < 0) {
                atom->taud[ia].y += 1.0;
}
            if (atom->taud[ia].z < 0) {
                atom->taud[ia].z += 1.0;
}
            if (atom->taud[ia].x >= 1.0) {
                atom->taud[ia].x -= 1.0;
}
            if (atom->taud[ia].y >= 1.0) {
                atom->taud[ia].y -= 1.0;
}
            if (atom->taud[ia].z >= 1.0) {
                atom->taud[ia].z -= 1.0;
}

            if (atom->taud[ia].x < 0 || atom->taud[ia].y < 0
                || atom->taud[ia].z < 0 || atom->taud[ia].x >= 1.0
                || atom->taud[ia].y >= 1.0 || atom->taud[ia].z >= 1.0) {
                GlobalV::ofs_warning << " it=" << it + 1 << " ia=" << ia + 1
                                     << std::endl;
                GlobalV::ofs_warning << "d=" << atom->taud[ia].x << " "
                                     << atom->taud[ia].y << " "
                                     << atom->taud[ia].z << std::endl;
                ModuleBase::WARNING_QUIT(
                    "Ions_Move_Basic::move_ions",
                    "the movement of atom is larger than the length of cell.");
            }

            atom->tau[ia] = atom->taud[ia] * this->latvec;
        }
    }
    return;
}

void UnitCell::bcast_atoms_tau() {
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < ntype; i++) {
        atoms[i].bcast_atom(); // bcast tau array
    }
#endif
}

void UnitCell::cal_ux() {
    double amag, uxmod;
    int starting_it = 0;
    int starting_ia = 0;
    bool is_paraller;
    // do not sign feature in teh general case
    magnet.lsign_ = false;
    ModuleBase::GlobalFunc::ZEROS(magnet.ux_, 3);

    for (int it = 0; it < ntype; it++) {
        for (int ia = 0; ia < atoms[it].na; ia++) {
            amag = pow(atoms[it].m_loc_[ia].x, 2)
                   + pow(atoms[it].m_loc_[ia].y, 2)
                   + pow(atoms[it].m_loc_[ia].z, 2);
            if (amag > 1e-6) {
                magnet.ux_[0] = atoms[it].m_loc_[ia].x;
                magnet.ux_[1] = atoms[it].m_loc_[ia].y;
                magnet.ux_[2] = atoms[it].m_loc_[ia].z;
                starting_it = it;
                starting_ia = ia;
                magnet.lsign_ = true;
                break;
            }
        }
        if (magnet.lsign_) {
            break;
}
    }
    // whether the initial magnetizations is parallel
    for (int it = starting_it; it < ntype; it++) {
        for (int ia = 0; ia < atoms[it].na; ia++) {
            if (it > starting_it || ia > starting_ia) {
                magnet.lsign_
                    = magnet.lsign_
                      && judge_parallel(magnet.ux_, atoms[it].m_loc_[ia]);
            }
        }
    }
    if (magnet.lsign_) {
        uxmod = pow(magnet.ux_[0], 2) + pow(magnet.ux_[1], 2)
                + pow(magnet.ux_[2], 2);
        if (uxmod < 1e-6) {
            ModuleBase::WARNING_QUIT("cal_ux", "wrong uxmod");
        }
        for (int i = 0; i < 3; i++) {
            magnet.ux_[i] *= 1 / sqrt(uxmod);
        }
        //       std::cout<<"    Fixed quantization axis for GGA: "
        //<<std::setw(10)<<ux[0]<<"  "<<std::setw(10)<<ux[1]<<"
        //"<<std::setw(10)<<ux[2]<<std::endl;
    }
    return;
}

bool UnitCell::judge_parallel(double a[3], ModuleBase::Vector3<double> b) {
    bool jp = false;
    double cross;
    cross = pow((a[1] * b.z - a[2] * b.y), 2)
            + pow((a[2] * b.x - a[0] * b.z), 2)
            + pow((a[0] * b.y - a[1] * b.x), 2);
    jp = (fabs(cross) < 1e-6);
    return jp;
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

    bool ok = true;
    bool ok2 = true;

    // (3) read in atom information
    if (GlobalV::MY_RANK == 0) {
        // open "atom_unitcell" file.
        std::ifstream ifa(fn.c_str(), std::ios::in);
        if (!ifa) {
            GlobalV::ofs_warning << fn;
            ok = false;
        }

        if (ok) {

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
            const int error = this->read_atom_species(ifa, log);

            //==========================
            // call read_atom_positions
            //==========================
            ok2 = this->read_atom_positions(ifa, log, GlobalV::ofs_warning);
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
    this->bcast_unitcell();
#endif

    //========================================================
    // Calculate unit cell volume
    // the reason to calculate volume here is
    // Firstly, latvec must be read in.
    //========================================================
    assert(lat0 > 0.0);
    this->omega = std::abs(latvec.Det()) * this->lat0 * lat0 * lat0;
    if (this->omega <= 0) {
        std::cout << "The volume is negative: " << this->omega << std::endl;
        ModuleBase::WARNING_QUIT("setup_cell", "omega <= 0 .");
    } else {
        log << std::endl;
        ModuleBase::GlobalFunc::OUT(log, "Volume (Bohr^3)", this->omega);
        ModuleBase::GlobalFunc::OUT(log,
                                    "Volume (A^3)",
                                    this->omega
                                        * pow(ModuleBase::BOHR_TO_A, 3));
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

void UnitCell::read_pseudo(std::ofstream& ofs) {
    // read in non-local pseudopotential and ouput the projectors.
    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
           ">>>>"
        << std::endl;
    ofs << " |                                                                 "
           "   |"
        << std::endl;
    ofs << " | Reading pseudopotentials files:                                 "
           "   |"
        << std::endl;
    ofs << " | The pseudopotential file is in UPF format. The 'NC' indicates "
           "that |"
        << std::endl;
    ofs << " | the type of pseudopotential is 'norm conserving'. Functional of "
           "   |"
        << std::endl;
    ofs << " | exchange and correlation is decided by 4 given parameters in "
           "UPF   |"
        << std::endl;
    ofs << " | file.  We also read in the 'core correction' if there exists.   "
           "   |"
        << std::endl;
    ofs << " | Also we can read the valence electrons number and the maximal   "
           "   |"
        << std::endl;
    ofs << " | angular momentum used in this pseudopotential. We also read in "
           "the |"
        << std::endl;
    ofs << " | trail wave function, trail atomic density and "
           "local-pseudopotential|"
        << std::endl;
    ofs << " | on logrithmic grid. The non-local pseudopotential projector is "
           "also|"
        << std::endl;
    ofs << " | read in if there is any.                                        "
           "   |"
        << std::endl;
    ofs << " |                                                                 "
           "   |"
        << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
           "<<<<"
        << std::endl;
    ofs << "\n\n\n\n";

    read_cell_pseudopots(PARAM.inp.pseudo_dir, ofs);

    if (GlobalV::MY_RANK == 0) {
        for (int it = 0; it < this->ntype; it++) {
            Atom* atom = &atoms[it];
            if (!(atom->label_orb.empty())) {
                compare_atom_labels(atom->label_orb, atom->ncpp.psd);
            }
        }

        if (PARAM.inp.out_element_info) {
            for (int i = 0; i < this->ntype; i++) {
                ModuleBase::Global_File::make_dir_atom(this->atoms[i].label);
            }
            for (int it = 0; it < ntype; it++) {
                Atom* atom = &atoms[it];
                std::stringstream ss;
                ss << PARAM.globalv.global_out_dir << atom->label << "/"
                   << atom->label << ".NONLOCAL";
                std::ofstream ofs(ss.str().c_str());

                ofs << "<HEADER>" << std::endl;
                ofs << std::setw(10) << atom->label << "\t"
                    << "label" << std::endl;
                ofs << std::setw(10) << atom->ncpp.pp_type << "\t"
                    << "pseudopotential type" << std::endl;
                ofs << std::setw(10) << atom->ncpp.lmax << "\t"
                    << "lmax" << std::endl;
                ofs << "</HEADER>" << std::endl;

                ofs << "\n<DIJ>" << std::endl;
                ofs << std::setw(10) << atom->ncpp.nbeta << "\t"
                    << "nummber of projectors." << std::endl;
                for (int ib = 0; ib < atom->ncpp.nbeta; ib++) {
                    for (int ib2 = 0; ib2 < atom->ncpp.nbeta; ib2++) {
                        ofs << std::setw(10) << atom->ncpp.lll[ib] << " "
                            << atom->ncpp.lll[ib2] << " "
                            << atom->ncpp.dion(ib, ib2) << std::endl;
                    }
                }
                ofs << "</DIJ>" << std::endl;

                for (int i = 0; i < atom->ncpp.nbeta; i++) {
                    ofs << "<PP_BETA>" << std::endl;
                    ofs << std::setw(10) << i << "\t"
                        << "the index of projectors." << std::endl;
                    ofs << std::setw(10) << atom->ncpp.lll[i] << "\t"
                        << "the angular momentum." << std::endl;

                    // mohan add
                    // only keep the nonzero part.
                    int cut_mesh = atom->ncpp.mesh;
                    for (int j = atom->ncpp.mesh - 1; j >= 0; --j) {
                        if (std::abs(atom->ncpp.betar(i, j)) > 1.0e-10) {
                            cut_mesh = j;
                            break;
                        }
                    }
                    if (cut_mesh % 2 == 0) {
                        ++cut_mesh;
}

                    ofs << std::setw(10) << cut_mesh << "\t"
                        << "the number of mesh points." << std::endl;

                    for (int j = 0; j < cut_mesh; ++j) {
                        ofs << std::setw(15) << atom->ncpp.r[j] << std::setw(15)
                            << atom->ncpp.betar(i, j) << std::setw(15)
                            << atom->ncpp.rab[j] << std::endl;
                    }
                    ofs << "</PP_BETA>" << std::endl;
                }

                ofs.close();
            }
        }
    }

#ifdef __MPI
    bcast_unitcell2();
#endif

    for (int it = 0; it < ntype; it++) {
        if (atoms[0].ncpp.xc_func != atoms[it].ncpp.xc_func) {
            GlobalV::ofs_warning << "\n type " << atoms[0].label
                                 << " functional is " << atoms[0].ncpp.xc_func;

            GlobalV::ofs_warning << "\n type " << atoms[it].label
                                 << " functional is " << atoms[it].ncpp.xc_func
                                 << std::endl;

            ModuleBase::WARNING_QUIT("setup_cell",
                                     "All DFT functional must consistent.");
        }
    }

    // setup the total number of PAOs
    cal_natomwfc(ofs);

    // Calculate the information of atoms from the pseudopotential to set PARAM
    CalAtomsInfo ca;
    ca.cal_atoms_info(this->atoms, this->ntype, PARAM);

    // setup PARAM.globalv.nlocal
    cal_nwfc(ofs);

    // Check whether the number of valence is minimum
    if (GlobalV::MY_RANK == 0) {
        int abtype = 0;
        for (int it = 0; it < ntype; it++) {
            if (ModuleBase::MinZval.find(atoms[it].ncpp.psd)
                != ModuleBase::MinZval.end()) {
                if (atoms[it].ncpp.zv
                    > ModuleBase::MinZval.at(atoms[it].ncpp.psd)) {
                    abtype += 1;
                    if (abtype == 1) {
                        std::cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                                     "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                                     "%%%%%%%%%%%%%%%%%%%%%%%%%%"
                                  << std::endl;
                        ofs << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                               "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                               "%%%%%%%%%%%%%%%%%%%%%"
                            << std::endl;
                    }
                    std::cout << " Warning: the number of valence electrons in "
                                 "pseudopotential > "
                              << ModuleBase::MinZval.at(atoms[it].ncpp.psd);
                    std::cout << " for " << atoms[it].ncpp.psd << ": "
                              << ModuleBase::EleConfig.at(atoms[it].ncpp.psd)
                              << std::endl;
                    ofs << " Warning: the number of valence electrons in "
                           "pseudopotential > "
                        << ModuleBase::MinZval.at(atoms[it].ncpp.psd);
                    ofs << " for " << atoms[it].ncpp.psd << ": "
                        << ModuleBase::EleConfig.at(atoms[it].ncpp.psd)
                        << std::endl;
                }
            }
        }
        if (abtype > 0) {
            std::cout << " Pseudopotentials with additional electrons can "
                         "yield (more) accurate outcomes, but may be "
                         "less efficient."
                      << std::endl;
            std::cout
                << " If you're confident that your chosen pseudopotential is "
                   "appropriate, you can safely ignore "
                   "this warning."
                << std::endl;
            std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%\n"
                      << std::endl;
            ofs << " Pseudopotentials with additional electrons can yield "
                   "(more) accurate outcomes, but may be less "
                   "efficient."
                << std::endl;
            ofs << " If you're confident that your chosen pseudopotential is "
                   "appropriate, you can safely ignore this "
                   "warning."
                << std::endl;
            ofs << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                   "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                   "%%%%%%%";
            ModuleBase::GlobalFunc::OUT(ofs, "");
        }
    }

    cal_meshx();

#ifdef __MPI
    Parallel_Common::bcast_int(meshx);
    Parallel_Common::bcast_int(natomwfc);
    Parallel_Common::bcast_int(lmax);
    Parallel_Common::bcast_int(lmax_ppwf);
#endif
}

//===========================================
// calculate the total number of local basis
// Target : nwfc, lmax,
// 			atoms[].stapos_wf
// 			PARAM.inp.nbands
//===========================================
void UnitCell::cal_nwfc(std::ofstream& log) {
    ModuleBase::TITLE("UnitCell", "cal_nwfc");
    assert(ntype > 0);
    assert(nat > 0);

    //===========================
    // (1) set iw2l, iw2n, iw2m
    //===========================
    for (int it = 0; it < ntype; it++) {
        this->atoms[it].set_index();
    }

    //===========================
    // (2) set namax and nwmax
    //===========================
    this->namax = 0;
    this->nwmax = 0;
    for (int it = 0; it < ntype; it++) {
        this->namax = std::max(atoms[it].na, namax);
        this->nwmax = std::max(atoms[it].nw, nwmax);
    }
    assert(namax > 0);
    // for tests
    //		OUT(GlobalV::ofs_running,"max input atom number",namax);
    //		OUT(GlobalV::ofs_running,"max wave function number",nwmax);

    //===========================
    // (3) set nwfc and stapos_wf
    //===========================
    int nlocal_tmp = 0;
    for (int it = 0; it < ntype; it++) {
        atoms[it].stapos_wf = nlocal_tmp;
        const int nlocal_it = atoms[it].nw * atoms[it].na;
        if (PARAM.inp.nspin != 4) {
            nlocal_tmp += nlocal_it;
        } else {
            nlocal_tmp += nlocal_it * 2; // zhengdy-soc
        }

        // for tests
        //		OUT(GlobalV::ofs_running,ss1.str(),nlocal_it);
        //		OUT(GlobalV::ofs_running,"start position of local
        //orbitals",atoms[it].stapos_wf);
    }

    // OUT(GlobalV::ofs_running,"NLOCAL",PARAM.globalv.nlocal);
    log << " " << std::setw(40) << "NLOCAL"
        << " = " << nlocal_tmp << std::endl;
    //========================================================
    // (4) set index for itia2iat, itiaiw2iwt
    //========================================================

    // mohan add 2010-09-26
    assert(nlocal_tmp > 0);
    assert(nlocal_tmp == PARAM.globalv.nlocal);
    delete[] iwt2iat;
    delete[] iwt2iw;
    this->iwt2iat = new int[nlocal_tmp];
    this->iwt2iw = new int[nlocal_tmp];

    this->itia2iat.create(ntype, namax);
    // this->itiaiw2iwt.create(ntype, namax, nwmax*PARAM.globalv.npol);
    this->set_iat2iwt(PARAM.globalv.npol);
    int iat = 0;
    int iwt = 0;
    for (int it = 0; it < ntype; it++) {
        for (int ia = 0; ia < atoms[it].na; ia++) {
            this->itia2iat(it, ia) = iat;
            // this->iat2ia[iat] = ia;
            for (int iw = 0; iw < atoms[it].nw * PARAM.globalv.npol; iw++) {
                // this->itiaiw2iwt(it, ia, iw) = iwt;
                this->iwt2iat[iwt] = iat;
                this->iwt2iw[iwt] = iw;
                ++iwt;
            }
            ++iat;
        }
    }

    //========================
    // (5) set lmax and nmax
    //========================
    this->lmax = 0;
    this->nmax = 0;
    for (int it = 0; it < ntype; it++) {
        lmax = std::max(lmax, atoms[it].nwl);
        for (int l = 0; l < atoms[it].nwl + 1; l++) {
            nmax = std::max(nmax, atoms[it].l_nchi[l]);
        }

        int nchi = 0;
        for (int l = 0; l < atoms[it].nwl + 1; l++) {
            nchi += atoms[it].l_nchi[l];
        }
        this->nmax_total = std::max(nmax_total, nchi);
    }

    //=======================
    // (6) set lmax_ppwf
    //=======================
    this->lmax_ppwf = 0;
    for (int it = 0; it < ntype; it++) {
        for (int ic = 0; ic < atoms[it].ncpp.nchi; ic++) {
            if (lmax_ppwf < atoms[it].ncpp.lchi[ic]) {
                this->lmax_ppwf = atoms[it].ncpp.lchi[ic];
            }
        }
    }

    /*
    for(int it=0; it< ntype; it++)
    {
        std::cout << " label=" << it << " nbeta=" << atoms[it].nbeta <<
    std::endl; for(int ic=0; ic<atoms[it].nbeta; ic++)
        {
            std::cout << " " << atoms[it].lll[ic] << std::endl;
        }
    }
    */

    //	OUT(GlobalV::ofs_running,"lmax between L(pseudopotential)",lmax_ppwf);

    //=====================
    // Use localized basis
    //=====================
    if ((PARAM.inp.basis_type == "lcao") || (PARAM.inp.basis_type == "lcao_in_pw")
        || ((PARAM.inp.basis_type == "pw") && (PARAM.inp.psi_initializer)
            && (PARAM.inp.init_wfc.substr(0, 3) == "nao")
            && (PARAM.inp.esolver_type == "ksdft"))) // xiaohui add 2013-09-02
    {
        ModuleBase::GlobalFunc::AUTO_SET("NBANDS", PARAM.inp.nbands);
    } else // plane wave basis
    {
        // if(winput::after_iter && winput::sph_proj)
        //{
        //	if(PARAM.inp.nbands < PARAM.globalv.nlocal)
        //	{
        //		ModuleBase::WARNING_QUIT("cal_nwfc","NBANDS must > PARAM.globalv.nlocal
        //!");
        //	}
        // }
    }

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

//======================
// Target : meshx
// Demand : atoms[].msh
//======================
void UnitCell::cal_meshx() {
    if (PARAM.inp.test_pseudo_cell) {
        ModuleBase::TITLE("UnitCell", "cal_meshx");
}
    this->meshx = 0;
    for (int it = 0; it < this->ntype; it++) {
        const int mesh = this->atoms[it].ncpp.msh;
        if (mesh > this->meshx) {
            this->meshx = mesh;
        }
    }
    return;
}

//=========================
// Target : natomwfc
// Demand : atoms[].nchi
// 			atoms[].lchi
// 			atoms[].oc
// 			atoms[].na
//=========================
void UnitCell::cal_natomwfc(std::ofstream& log) {
    if (PARAM.inp.test_pseudo_cell) {
        ModuleBase::TITLE("UnitCell", "cal_natomwfc");
}

    this->natomwfc = 0;
    for (int it = 0; it < ntype; it++) {
        //============================
        // Use pseudo-atomic orbitals
        //============================
        int tmp = 0;
        for (int l = 0; l < atoms[it].ncpp.nchi; l++) {
            if (atoms[it].ncpp.oc[l] >= 0) {
                if (PARAM.inp.nspin == 4) {
                    if (atoms[it].ncpp.has_so) {
                        tmp += 2 * atoms[it].ncpp.lchi[l];
                        if (fabs(atoms[it].ncpp.jchi[l] - atoms[it].ncpp.lchi[l]
                                 - 0.5)
                            < 1e-6) {
                            tmp += 2;
}
                    } else {
                        tmp += 2 * (2 * atoms[it].ncpp.lchi[l] + 1);
                    }
                } else {
                    tmp += 2 * atoms[it].ncpp.lchi[l] + 1;
}
            }
        }
        natomwfc += tmp * atoms[it].na;
    }
    ModuleBase::GlobalFunc::OUT(log,
                                "initial pseudo atomic orbital number",
                                natomwfc);
    return;
}

// LiuXh add a new function here,
// 20180515
void UnitCell::setup_cell_after_vc(std::ofstream& log) {
    ModuleBase::TITLE("UnitCell", "setup_cell_after_vc");
    assert(lat0 > 0.0);
    this->omega = std::abs(latvec.Det()) * this->lat0 * lat0 * lat0;
    if (this->omega <= 0) {
        ModuleBase::WARNING_QUIT("setup_cell_after_vc", "omega <= 0 .");
    } else {
        log << std::endl;
        ModuleBase::GlobalFunc::OUT(log, "Volume (Bohr^3)", this->omega);
        ModuleBase::GlobalFunc::OUT(log,
                                    "Volume (A^3)",
                                    this->omega
                                        * pow(ModuleBase::BOHR_TO_A, 3));
    }

    lat0_angstrom = lat0 * 0.529177;
    tpiba = ModuleBase::TWO_PI / lat0;
    tpiba2 = tpiba * tpiba;

    // lattice vectors in another form.
    a1.x = latvec.e11;
    a1.y = latvec.e12;
    a1.z = latvec.e13;

    a2.x = latvec.e21;
    a2.y = latvec.e22;
    a2.z = latvec.e23;

    a3.x = latvec.e31;
    a3.y = latvec.e32;
    a3.z = latvec.e33;

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec has the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    this->GT = latvec.Inverse();
    this->G = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    for (int it = 0; it < ntype; it++) {
        Atom* atom = &atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            atom->tau[ia] = atom->taud[ia] * latvec;
        }
    }

#ifdef __MPI
    this->bcast_unitcell();
#endif

    log << std::endl;
    output::printM3(log,
                    "Lattice vectors: (Cartesian coordinate: in unit of a_0)",
                    latvec);
    output::printM3(
        log,
        "Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",
        G);

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

void UnitCell::remake_cell() {
    ModuleBase::TITLE("UnitCell", "rmake_cell");

    // The idea is as follows: for each type of lattice, first calculate
    // from current latvec the lattice parameters, then use the parameters
    // to reconstruct latvec

    if (latName == "none") {
        ModuleBase::WARNING_QUIT(
            "UnitCell",
            "to use fixed_ibrav, latname must be provided");
    } else if (latName == "sc") // ibrav = 1
    {
        double celldm = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                  + pow(latvec.e13, 2));

        latvec.Zero();
        latvec.e11 = latvec.e22 = latvec.e33 = celldm;
    } else if (latName == "fcc") // ibrav = 2
    {
        double celldm = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                  + pow(latvec.e13, 2))
                        / std::sqrt(2.0);

        latvec.e11 = -celldm;
        latvec.e12 = 0.0;
        latvec.e13 = celldm;
        latvec.e21 = 0.0;
        latvec.e22 = celldm;
        latvec.e23 = celldm;
        latvec.e31 = -celldm;
        latvec.e32 = celldm;
        latvec.e33 = 0.0;
    } else if (latName == "bcc") // ibrav = 3
    {
        double celldm = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                  + pow(latvec.e13, 2))
                        / std::sqrt(3.0);

        latvec.e11 = celldm;
        latvec.e12 = celldm;
        latvec.e13 = celldm;
        latvec.e21 = -celldm;
        latvec.e22 = celldm;
        latvec.e23 = celldm;
        latvec.e31 = -celldm;
        latvec.e32 = -celldm;
        latvec.e33 = celldm;
    } else if (latName == "hexagonal") // ibrav = 4
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));
        double e22 = sqrt(3.0) / 2.0;

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = -0.5 * celldm1;
        latvec.e22 = celldm1 * e22;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "trigonal") // ibrav = 5
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22
                           + latvec.e13 * latvec.e23);
        double cos12 = celldm12 / celldm1 / celldm2;

        if (cos12 <= -0.5 || cos12 >= 1.0) {
            ModuleBase::WARNING_QUIT("unitcell", "wrong cos12!");
        }
        double t1 = sqrt(1.0 + 2.0 * cos12);
        double t2 = sqrt(1.0 - cos12);

        double e11 = celldm1 * t2 / sqrt(2.0);
        double e12 = -celldm1 * t2 / sqrt(6.0);
        double e13 = celldm1 * t1 / sqrt(3.0);
        double e22 = celldm1 * sqrt(2.0) * t2 / sqrt(3.0);

        latvec.e11 = e11;
        latvec.e12 = e12;
        latvec.e13 = e13;
        latvec.e21 = 0.0;
        latvec.e22 = e22;
        latvec.e23 = e13;
        latvec.e31 = -e11;
        latvec.e32 = e12;
        latvec.e33 = e13;
    } else if (latName == "st") // ibrav = 6
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));
        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = 0.0;
        latvec.e22 = celldm1;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "bct") // ibrav = 7
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::abs(latvec.e13);

        latvec.e11 = celldm1;
        latvec.e12 = -celldm1;
        latvec.e13 = celldm2;
        latvec.e21 = celldm1;
        latvec.e22 = celldm1;
        latvec.e23 = celldm2;
        latvec.e31 = -celldm1;
        latvec.e32 = -celldm1;
        latvec.e33 = celldm2;
    } else if (latName == "so") // ibrav = 8
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = 0.0;
        latvec.e22 = celldm2;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "baco") // ibrav = 9
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::abs(latvec.e22);
        double celldm3 = std::abs(latvec.e33);

        latvec.e11 = celldm1;
        latvec.e12 = celldm2;
        latvec.e13 = 0.0;
        latvec.e21 = -celldm1;
        latvec.e22 = celldm2;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "fco") // ibrav = 10
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::abs(latvec.e22);
        double celldm3 = std::abs(latvec.e33);

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = celldm3;
        latvec.e21 = celldm1;
        latvec.e22 = celldm2;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = celldm2;
        latvec.e33 = celldm3;
    } else if (latName == "bco") // ibrav = 11
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::abs(latvec.e12);
        double celldm3 = std::abs(latvec.e13);

        latvec.e11 = celldm1;
        latvec.e12 = celldm2;
        latvec.e13 = celldm3;
        latvec.e21 = -celldm1;
        latvec.e22 = celldm2;
        latvec.e23 = celldm3;
        latvec.e31 = -celldm1;
        latvec.e32 = -celldm2;
        latvec.e33 = celldm3;
    } else if (latName == "sm") // ibrav = 12
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));
        double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22
                           + latvec.e13 * latvec.e23);
        double cos12 = celldm12 / celldm1 / celldm2;

        double e21 = celldm2 * cos12;
        double e22 = celldm2 * std::sqrt(1.0 - cos12 * cos12);

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = e21;
        latvec.e22 = e22;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "bacm") // ibrav = 13
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm3 = std::abs(latvec.e13);

        double cos12 = latvec.e21 / celldm2;
        if (cos12 >= 1.0) {
            ModuleBase::WARNING_QUIT("unitcell", "wrong cos12!");
        }

        double e21 = celldm2 * cos12;
        double e22 = celldm2 * std::sqrt(1.0 - cos12 * cos12);

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = -celldm3;
        latvec.e21 = e21;
        latvec.e22 = e22;
        latvec.e23 = 0.0;
        latvec.e31 = celldm1;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "triclinic") // ibrav = 14
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));
        double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22
                           + latvec.e13 * latvec.e23);
        double cos12 = celldm12 / celldm1 / celldm2;
        double celldm13 = (latvec.e11 * latvec.e31 + latvec.e12 * latvec.e32
                           + latvec.e13 * latvec.e33);
        double cos13 = celldm13 / celldm1 / celldm3;
        double celldm23 = (latvec.e21 * latvec.e31 + latvec.e22 * latvec.e32
                           + latvec.e23 * latvec.e33);
        double cos23 = celldm23 / celldm2 / celldm3;

        double sin12 = std::sqrt(1.0 - cos12 * cos12);
        if (cos12 >= 1.0) {
            ModuleBase::WARNING_QUIT("unitcell", "wrong cos12!");
        }

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = celldm2 * cos12;
        latvec.e22 = celldm2 * sin12;
        latvec.e23 = 0.0;
        latvec.e31 = celldm3 * cos13;
        latvec.e32 = celldm3 * (cos23 - cos13 * cos12) / sin12;
        double term = 1.0 + 2.0 * cos12 * cos13 * cos23 - cos12 * cos12
                      - cos13 * cos13 - cos23 * cos23;
        term = sqrt(term) / sin12;
        latvec.e33 = celldm3 * term;
    } else {
        std::cout << "latname is : " << latName << std::endl;
        ModuleBase::WARNING_QUIT("UnitCell::read_atom_species",
                                 "latname not supported!");
    }
}

void cal_nelec(const Atom* atoms, const int& ntype, double& nelec)
{
    ModuleBase::TITLE("UnitCell", "cal_nelec");
    GlobalV::ofs_running << "\n SETUP THE ELECTRONS NUMBER" << std::endl;

    if (nelec == 0)
    {
        if (PARAM.inp.use_paw)
        {
#ifdef USE_PAW
            for (int it = 0; it < ntype; it++)
            {
                std::stringstream ss1, ss2;
                ss1 << " electron number of element " << GlobalC::paw_cell.get_zat(it) << std::endl;
                const int nelec_it = GlobalC::paw_cell.get_val(it) * atoms[it].na;
                nelec += nelec_it;
                ss2 << "total electron number of element " << GlobalC::paw_cell.get_zat(it);

                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss1.str(), GlobalC::paw_cell.get_val(it));
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss2.str(), nelec_it);
            }
#endif
        }
        else
        {
            for (int it = 0; it < ntype; it++)
            {
                std::stringstream ss1, ss2;
                ss1 << "electron number of element " << atoms[it].label;
                const double nelec_it = atoms[it].ncpp.zv * atoms[it].na;
                nelec += nelec_it;
                ss2 << "total electron number of element " << atoms[it].label;

                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss1.str(), atoms[it].ncpp.zv);
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss2.str(), nelec_it);
            }
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "AUTOSET number of electrons: ", nelec);
        }
    }
    if (PARAM.inp.nelec_delta != 0)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                    "nelec_delta is NOT zero, please make sure you know what you are "
                                    "doing! nelec_delta: ",
                                    PARAM.inp.nelec_delta);
        nelec += PARAM.inp.nelec_delta;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec now: ", nelec);
    }
    return;
}

void cal_nbands(const int& nelec, const int& nlocal, const std::vector<double>& nelec_spin, int& nbands)
{
    if (PARAM.inp.esolver_type == "sdft") // qianrui 2021-2-20
    {
        return;
    }
    //=======================================
    // calculate number of bands (setup.f90)
    //=======================================
    double occupied_bands = static_cast<double>(nelec / ModuleBase::DEGSPIN);
    if (PARAM.inp.lspinorb == 1) {
        occupied_bands = static_cast<double>(nelec);
    }

    if ((occupied_bands - std::floor(occupied_bands)) > 0.0)
    {
        occupied_bands = std::floor(occupied_bands) + 1.0; // mohan fix 2012-04-16
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "occupied bands", occupied_bands);

    if (nbands == 0)
    {
        if (PARAM.inp.nspin == 1)
        {
            const int nbands1 = static_cast<int>(occupied_bands) + 10;
            const int nbands2 = static_cast<int>(1.2 * occupied_bands) + 1;
            nbands = std::max(nbands1, nbands2);
            if (PARAM.inp.basis_type != "pw") {
                nbands = std::min(nbands, nlocal);
            }
        }
        else if (PARAM.inp.nspin == 4)
        {
            const int nbands3 = nelec + 20;
            const int nbands4 = static_cast<int>(1.2 * nelec) + 1;
            nbands = std::max(nbands3, nbands4);
            if (PARAM.inp.basis_type != "pw") {
                nbands = std::min(nbands, nlocal);
            }
        }
        else if (PARAM.inp.nspin == 2)
        {
            const double max_occ = std::max(nelec_spin[0], nelec_spin[1]);
            const int nbands3 = static_cast<int>(max_occ) + 11;
            const int nbands4 = static_cast<int>(1.2 * max_occ) + 1;
            nbands = std::max(nbands3, nbands4);
            if (PARAM.inp.basis_type != "pw") {
                nbands = std::min(nbands, nlocal);
            }
        }
        ModuleBase::GlobalFunc::AUTO_SET("NBANDS", nbands);
    }
    // else if ( PARAM.inp.calculation=="scf" || PARAM.inp.calculation=="md" || PARAM.inp.calculation=="relax") //pengfei
    // 2014-10-13
    else
    {
        if (nbands < occupied_bands) {
            ModuleBase::WARNING_QUIT("unitcell", "Too few bands!");
        }
        if (PARAM.inp.nspin == 2)
        {
            if (nbands < nelec_spin[0])
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec_up", nelec_spin[0]);
                ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "Too few spin up bands!");
            }
            if (nbands < nelec_spin[1])
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec_down", nelec_spin[1]);
                ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "Too few spin down bands!");
            }
        }
    }

    // mohan add 2010-09-04
    // std::cout << "nbands(this-> = " <<nbands <<std::endl;
    if (nbands == occupied_bands)
    {
        if (PARAM.inp.smearing_method != "fixed")
        {
            ModuleBase::WARNING_QUIT("ElecState::cal_nbands", "for smearing, num. of bands > num. of occupied bands");
        }
    }

    // mohan update 2021-02-19
    // mohan add 2011-01-5
    if (PARAM.inp.basis_type == "lcao" || PARAM.inp.basis_type == "lcao_in_pw")
    {
        if (nbands > nlocal)
        {
            ModuleBase::WARNING_QUIT("ElecState::cal_nbandsc", "NLOCAL < NBANDS");
        }
        else
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NLOCAL", nlocal);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NBANDS", nbands);
        }
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "NBANDS", nbands);
}

void UnitCell::compare_atom_labels(std::string label1, std::string label2) {
    if (label1
        != label2) //'!( "Ag" == "Ag" || "47" == "47" || "Silver" == Silver" )'
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