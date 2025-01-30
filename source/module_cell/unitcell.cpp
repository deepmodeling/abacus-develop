#include <cstdlib>
#ifdef __MPI
#include "mpi.h"
#endif

#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "unitcell.h"
#include "bcast_cell.h"
#include "module_parameter/parameter.h"
#include "read_stru.h"
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

/*
void UnitCell::print_cell_xyz(const std::string& fn) const
{

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
    this->nmax_total = 0;
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
        || ((PARAM.inp.basis_type == "pw") && (PARAM.inp.init_wfc.substr(0, 3) == "nao")
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