#include <cstdlib>
#ifdef __MPI
#include "mpi.h"
#endif

#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "unitcell.h"

#ifdef __LCAO
#include "../module_basis/module_ao/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02
#include "module_base/global_file.h"
#include "module_base/parallel_common.h"
#include "module_base/element_elec_config.h"
#include "module_base/element_covalent_radius.h"
#include "module_base/atom_in.h"

#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif
#ifdef __EXX
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_ri/serialization_cereal.h"
#endif

UnitCell::UnitCell()
{
    if (GlobalV::test_unitcell)
        ModuleBase::TITLE("unitcell", "Constructor");
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

UnitCell::~UnitCell()
{
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
    if(set_atom_flag)
	{
		delete[] atoms;
	}
}

#include "module_base/parallel_common.h"
#ifdef __MPI
void UnitCell::bcast_unitcell(void)
{
    if (GlobalV::test_unitcell)
        ModuleBase::TITLE("UnitCell", "bcast_unitcell");
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

    if (GlobalV::NSPIN == 4)
    {
        Parallel_Common::bcast_double(magnet.ux_[0]);
        Parallel_Common::bcast_double(magnet.ux_[1]);
        Parallel_Common::bcast_double(magnet.ux_[2]);
    }

    for (int i = 0; i < ntype; i++)
    {
        atoms[i].bcast_atom(); // init tau and mbl array
    }

#ifdef __EXX
    ModuleBase::bcast_data_cereal(GlobalC::exx_info.info_ri.files_abfs, MPI_COMM_WORLD, 0);
#endif
    return;
}

void UnitCell::bcast_unitcell2(void)
{
    for (int i = 0; i < ntype; i++)
    {
        atoms[i].bcast_atom2();
    }
    return;
}
#endif

void UnitCell::print_cell(std::ofstream& ofs) const
{
    if (GlobalV::test_unitcell)
        ModuleBase::TITLE("UnitCell", "print_cell");

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

void UnitCell::print_cell_cif(const std::string& fn) const
{
    if (GlobalV::test_unitcell)
        ModuleBase::TITLE("UnitCell", "print_cell_cif");

    if (GlobalV::MY_RANK != 0)
        return; // xiaohui add 2015-03-15

    std::stringstream ss;
    ss << GlobalV::global_out_dir << fn;

    std::ofstream ofs(ss.str().c_str());
    ofs << "data_" << latName << std::endl;
    ofs << std::endl;
    // xiaohui modify 2015-03-25
    // ofs << "_audit_creation_method" << " generated by MESIA" << std::endl;
    ofs << "_audit_creation_method"
        << " generated by ABACUS" << std::endl;
    ofs << std::endl;
    ofs << "_cell_length_a " << a1.norm() * lat0 * 0.529177 << std::endl;
    ofs << "_cell_length_b " << a2.norm() * lat0 * 0.529177 << std::endl;
    ofs << "_cell_length_c " << a3.norm() * lat0 * 0.529177 << std::endl;

    // xiaohui modify and add 2014-12-21
    // ofs << "_cell_angle_alpha " << "90.00" << std::endl;
    // ofs << "_cell_angle_beta " << "90.00" << std::endl;
    // ofs << "_cell_angle_gamma " << "90.00" << std::endl;
    // xiaohui modify alpha, beta and gamma 2015-09-29
    double angle_alpha = acos((a2 * a3) / (a2.norm() * a3.norm())) / ModuleBase::PI * 180.0;
    double angle_beta = acos((a1 * a3) / (a1.norm() * a3.norm())) / ModuleBase::PI * 180.0;
    double angle_gamma = acos((a1 * a2) / (a1.norm() * a2.norm())) / ModuleBase::PI * 180.0;
    ofs << "_cell_angle_alpha " << angle_alpha << std::endl;
    ofs << "_cell_angle_beta " << angle_beta << std::endl;
    ofs << "_cell_angle_gamma " << angle_gamma << std::endl;
    ofs << std::endl;
    // ofs << "_symmetry_space_group_name_H-M"
    //     << " " << std::endl;
    // ofs << "_symmetry_Int_Tables_number"
    //     << " " << std::endl;
    // ofs << std::endl;
    ofs << "loop_" << std::endl;
    ofs << "_atom_site_label" << std::endl;
    ofs << "_atom_site_fract_x" << std::endl;
    ofs << "_atom_site_fract_y" << std::endl;
    ofs << "_atom_site_fract_z" << std::endl;
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            ofs << atoms[it].label << " " << atoms[it].taud[ia].x << " " << atoms[it].taud[ia].y << " "
                << atoms[it].taud[ia].z << std::endl;
        }
    }
    ofs.close();
}

/*
void UnitCell::print_cell_xyz(const std::string& fn) const
{
    if (GlobalV::test_unitcell)
        ModuleBase::TITLE("UnitCell", "print_cell_xyz");

    if (GlobalV::MY_RANK != 0)
        return; // xiaohui add 2015-03-15

    std::stringstream ss;
    ss << GlobalV::global_out_dir << fn;

    std::ofstream ofs(ss.str().c_str());

    ofs << nat << std::endl;
    ofs << latName << std::endl;
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            ofs << atoms[it].label << " " << atoms[it].tau[ia].x * lat0 * 0.529177 << " "
                << atoms[it].tau[ia].y * lat0 * 0.529177 << " " << atoms[it].tau[ia].z * lat0 * 0.529177 << std::endl;
        }
    }

    ofs.close();
    return;
}
*/

void UnitCell::set_iat2itia(void)
{
    assert(nat > 0);
    delete[] iat2it;
    delete[] iat2ia;
    this->iat2it = new int[nat];
    this->iat2ia = new int[nat];
    int iat = 0;
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            this->iat2it[iat] = it;
            this->iat2ia[iat] = ia;
            ++iat;
        }
    }
    return;
}

std::map<int, int> UnitCell::get_atomCounts() const
{
	std::map<int, int> atomCounts;
	for (int it = 0; it < this->ntype; it++)
	{
		atomCounts.insert(std::pair<int, int>(it, this->atoms[it].na));
	}
	return atomCounts;
}

std::map<int, int> UnitCell::get_orbitalCounts() const
{
	std::map<int, int> orbitalCounts;
	for (int it = 0; it < this->ntype; it++)
	{
		orbitalCounts.insert(std::pair<int, int>(it, this->atoms[it].nw));
	}
	return orbitalCounts;
}

std::map<int, std::map<int, int>> UnitCell::get_lnchiCounts() const
{
    std::map<int, std::map<int, int>> lnchiCounts;
    for (int it = 0; it < this->ntype; it++)
    {
        for (int L = 0; L < this->atoms[it].nwl + 1; L++)
        {
            // Check if the key 'it' exists in the outer map
            if (lnchiCounts.find(it) == lnchiCounts.end())
            {
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

void UnitCell::update_pos_tau(const double* pos)
{
    int iat = 0;
    for (int it = 0; it < this->ntype; it++)
    {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++)
        {
            for ( int ik = 0; ik < 3; ++ik)
            {
                if (atom->mbl[ia][ik])
                {
                    atom->dis[ia][ik] = pos[3 * iat + ik] / this->lat0 - atom->tau[ia][ik];
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

void UnitCell::update_pos_taud(double* posd_in)
{
    int iat = 0;
    for (int it = 0; it < this->ntype; it++)
    {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++)
        {
            for ( int ik = 0; ik < 3; ++ik)
            {
                atom->taud[ia][ik] += posd_in[3*iat + ik];
                atom->dis[ia][ik] = posd_in[3*iat + ik];
            }
            iat++;
        }
    }
    assert(iat == this->nat);
    this->periodic_boundary_adjustment();
    this->bcast_atoms_tau();
}

// posd_in is atomic displacements here  liuyu 2023-03-22
void UnitCell::update_pos_taud(const ModuleBase::Vector3<double>* posd_in)
{
    int iat = 0;
    for (int it = 0; it < this->ntype; it++)
    {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++)
        {
            for ( int ik = 0; ik < 3; ++ik)
            {
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

void UnitCell::update_vel(const ModuleBase::Vector3<double>* vel_in)
{
    int iat = 0;
    for (int it = 0; it < this->ntype; ++it)
    {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ++ia)
        {
            this->atoms[it].vel[ia] = vel_in[iat];
            ++iat;
        }
    }
    assert(iat == this->nat);
}

void UnitCell::periodic_boundary_adjustment()
{
    //----------------------------------------------
    // because of the periodic boundary condition
    // we need to adjust the atom positions,
    // first adjust direct coordinates,
    // then update them into cartesian coordinates,
    //----------------------------------------------
    for (int it = 0; it < this->ntype; it++)
    {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++)
        {
            // mohan update 2011-03-21
            if (atom->taud[ia].x < 0)
                atom->taud[ia].x += 1.0;
            if (atom->taud[ia].y < 0)
                atom->taud[ia].y += 1.0;
            if (atom->taud[ia].z < 0)
                atom->taud[ia].z += 1.0;
            if (atom->taud[ia].x >= 1.0)
                atom->taud[ia].x -= 1.0;
            if (atom->taud[ia].y >= 1.0)
                atom->taud[ia].y -= 1.0;
            if (atom->taud[ia].z >= 1.0)
                atom->taud[ia].z -= 1.0;

            if (atom->taud[ia].x < 0 || atom->taud[ia].y < 0 || atom->taud[ia].z < 0 || atom->taud[ia].x >= 1.0
                || atom->taud[ia].y >= 1.0 || atom->taud[ia].z >= 1.0)
            {
                GlobalV::ofs_warning << " it=" << it + 1 << " ia=" << ia + 1 << std::endl;
                GlobalV::ofs_warning << "d=" << atom->taud[ia].x << " " << atom->taud[ia].y << " " << atom->taud[ia].z
                                     << std::endl;
                ModuleBase::WARNING_QUIT("Ions_Move_Basic::move_ions",
                                         "the movement of atom is larger than the length of cell.");
            }

            atom->tau[ia] = atom->taud[ia] * this->latvec;
        }
    }
    return;
}

void UnitCell::bcast_atoms_tau()
{
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < ntype; i++)
    {
        atoms[i].bcast_atom(); // bcast tau array
    }
#endif
}

void UnitCell::cal_ux()
{
    double amag, uxmod;
    int starting_it = 0;
    int starting_ia = 0;
    bool is_paraller;
    // do not sign feature in teh general case
    magnet.lsign_ = false;
    ModuleBase::GlobalFunc::ZEROS(magnet.ux_, 3);

    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            amag = pow(atoms[it].m_loc_[ia].x, 2) + pow(atoms[it].m_loc_[ia].y, 2) + pow(atoms[it].m_loc_[ia].z, 2);
            if (amag > 1e-6)
            {
                magnet.ux_[0] = atoms[it].m_loc_[ia].x;
                magnet.ux_[1] = atoms[it].m_loc_[ia].y;
                magnet.ux_[2] = atoms[it].m_loc_[ia].z;
                starting_it = it;
                starting_ia = ia;
                magnet.lsign_ = true;
                break;
            }
        }
        if (magnet.lsign_)
            break;
    }
    // whether the initial magnetizations is parallel
    for (int it = starting_it; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            if (it > starting_it || ia > starting_ia)
            {
                magnet.lsign_ = magnet.lsign_ && judge_parallel(magnet.ux_, atoms[it].m_loc_[ia]);
            }
        }
    }
    if (magnet.lsign_)
    {
        uxmod = pow(magnet.ux_[0], 2) + pow(magnet.ux_[1], 2) + pow(magnet.ux_[2], 2);
        if (uxmod < 1e-6)
        {
            ModuleBase::WARNING_QUIT("cal_ux", "wrong uxmod");
        }
        for (int i = 0; i < 3; i++)
        {
            magnet.ux_[i] *= 1 / sqrt(uxmod);
        }
        //       std::cout<<"    Fixed quantization axis for GGA: "
        //<<std::setw(10)<<ux[0]<<"  "<<std::setw(10)<<ux[1]<<"  "<<std::setw(10)<<ux[2]<<std::endl;
    }
    return;
}

bool UnitCell::judge_parallel(double a[3], ModuleBase::Vector3<double> b)
{
    bool jp = false;
    double cross;
    cross = pow((a[1] * b.z - a[2] * b.y), 2) + pow((a[2] * b.x - a[0] * b.z), 2) + pow((a[0] * b.y - a[1] * b.x), 2);
    jp = (fabs(cross) < 1e-6);
    return jp;
}


//==============================================================
//Calculate various lattice related quantities for given latvec
//==============================================================
void UnitCell::setup_cell(const std::string &fn, std::ofstream &log)
{
	ModuleBase::TITLE("UnitCell","setup_cell");	
	// (1) init mag
	assert(ntype>0);
	delete[] magnet.start_magnetization;
	magnet.start_magnetization = new double[this->ntype];

	// (2) init *Atom class array.
	this->atoms = new Atom[this->ntype]; // atom species.
	this->set_atom_flag = true;


	bool ok = true;
	bool ok2 = true;

	// (3) read in atom information
	if(GlobalV::MY_RANK == 0)
	{
		// open "atom_unitcell" file.
		std::ifstream ifa(fn.c_str(), std::ios::in);
		if (!ifa)
		{
			GlobalV::ofs_warning << fn;
			ok = false;
		}

		if(ok)
		{

			log << "\n\n\n\n";
			log << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			log << " |                                                                    |" << std::endl;
			log << " | Reading atom information in unitcell:                              |" << std::endl;
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
	if(GlobalV::NSPIN==4) 
	{
		Parallel_Common::bcast_bool(GlobalV::DOMAG);
		Parallel_Common::bcast_bool(GlobalV::DOMAG_Z);
	}
#endif
	if(!ok)
	{
		ModuleBase::WARNING_QUIT("UnitCell::setup_cell","Can not find the file containing atom positions.!");
	}
	if(!ok2)
	{
		ModuleBase::WARNING_QUIT("UnitCell::setup_cell","Something wrong during read_atom_positions.");
	}

#ifdef __MPI
	this->bcast_unitcell();
#endif

	//after read STRU, calculate initial total magnetization when NSPIN=2
	if(GlobalV::NSPIN == 2 && !GlobalV::TWO_EFERMI)
	{
		for(int it = 0;it<this->ntype; it++)
		{
			for(int ia = 0; ia<this->atoms[it].na; ia++)
			{
				GlobalV::nupdown += this->atoms[it].mag[ia];
			}
		}
		GlobalV::ofs_running<<" The readin total magnetization is "<<GlobalV::nupdown<<std::endl;
	}
	
	//========================================================
	// Calculate unit cell volume
	// the reason to calculate volume here is 
	// Firstly, latvec must be read in.
	//========================================================
	assert(lat0 > 0.0);
	this->omega = std::abs( latvec.Det() ) * this->lat0 * lat0 * lat0 ;
	if(this->omega<=0)
	{	
		std::cout << "The volume is negative: " << this->omega<<std::endl;
		ModuleBase::WARNING_QUIT("setup_cell","omega <= 0 .");
	}
	else
	{
		log << std::endl;
		ModuleBase::GlobalFunc::OUT(log,"Volume (Bohr^3)", this->omega);
		ModuleBase::GlobalFunc::OUT(log,"Volume (A^3)", this->omega * pow(ModuleBase::BOHR_TO_A, 3));
	}
		
	//==========================================================
	// Calculate recip. lattice vectors and dot products
	// latvec have the unit of lat0, but G has the unit 2Pi/lat0
	//==========================================================
	this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
	this->invGGT = GGT.Inverse();

    //LiuXh add 20180515
    this->GT0 = latvec.Inverse();
    this->G0  = GT.Transpose();
    this->GGT0 = G * GT;
    this->invGGT0 = GGT.Inverse();

	log << std::endl;
	output::printM3(log,"Lattice vectors: (Cartesian coordinate: in unit of a_0)",latvec); 
	output::printM3(log,"Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",G);
//	OUT(log,"lattice center x",latcenter.x);
//	OUT(log,"lattice center y",latcenter.y);
//	OUT(log,"lattice center z",latcenter.z);

    //===================================
    // set index for iat2it, iat2ia
    //===================================
    this->set_iat2itia();

#ifdef USE_PAW
	if(GlobalV::use_paw)
	{
		GlobalC::paw_cell.set_libpaw_cell(latvec, lat0);

		int * typat;
		double * xred;

		typat = new int[nat];
		xred = new double[nat*3];

		int iat = 0;
		for(int it = 0; it < ntype; it ++)
		{
			for(int ia = 0; ia < atoms[it].na; ia ++)
			{
				typat[iat] = it + 1; //Fortran index starts from 1 !!!!
				xred[iat*3+0] = atoms[it].taud[ia].x;
				xred[iat*3+1] = atoms[it].taud[ia].y;
				xred[iat*3+2] = atoms[it].taud[ia].z;
				iat ++;
			}
		}

		GlobalC::paw_cell.set_libpaw_atom(nat,ntype,typat,xred);
		delete[] typat;
		delete[] xred;

		GlobalC::paw_cell.set_libpaw_files();

		GlobalC::paw_cell.set_nspin(GlobalV::NSPIN);
	}
#endif

    return;
}

void UnitCell::read_pseudo(std::ofstream &ofs)
{
    // read in non-local pseudopotential and ouput the projectors.
    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " | Reading pseudopotentials files:                                    |" << std::endl;
    ofs << " | The pseudopotential file is in UPF format. The 'NC' indicates that |" << std::endl;
    ofs << " | the type of pseudopotential is 'norm conserving'. Functional of    |" << std::endl;
    ofs << " | exchange and correlation is decided by 4 given parameters in UPF   |" << std::endl;
    ofs << " | file.  We also read in the 'core correction' if there exists.      |" << std::endl;
    ofs << " | Also we can read the valence electrons number and the maximal      |" << std::endl;
    ofs << " | angular momentum used in this pseudopotential. We also read in the |" << std::endl;
    ofs << " | trail wave function, trail atomic density and local-pseudopotential|" << std::endl;
    ofs << " | on logrithmic grid. The non-local pseudopotential projector is also|" << std::endl;
    ofs << " | read in if there is any.                                           |" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    ofs << "\n\n\n\n";

    read_cell_pseudopots(GlobalV::global_pseudo_dir, ofs);

    if(GlobalV::MY_RANK == 0)
    {
        for (int it = 0; it < this->ntype; it++)
        {
            Atom* atom = &atoms[it];
            if (!(atom->label_orb.empty()))
            {
                compare_atom_labels(atom->label_orb, atom->ncpp.psd);
            }
        }

        if(GlobalV::out_element_info)
        { 
            for(int i=0;i<this->ntype;i++)
            {
            	ModuleBase::Global_File::make_dir_atom( this->atoms[i].label );
            }
            for(int it=0; it<ntype; it++)
            {
                Atom* atom = &atoms[it];
                std::stringstream ss;
                ss << GlobalV::global_out_dir << atom->label 
                    << "/" << atom->label
                    << ".NONLOCAL";
                std::ofstream ofs(ss.str().c_str());
    
                ofs << "<HEADER>" << std::endl;
                ofs << std::setw(10) << atom->label << "\t" << "label" << std::endl;
                ofs << std::setw(10) << atom->ncpp.pp_type << "\t" << "pseudopotential type" << std::endl;
                ofs << std::setw(10) << atom->ncpp.lmax << "\t" << "lmax" << std::endl;
                ofs << "</HEADER>" << std::endl;
    
                ofs << "\n<DIJ>" << std::endl;
                ofs << std::setw(10) << atom->ncpp.nbeta << "\t" << "nummber of projectors." << std::endl;
                for(int ib=0; ib<atom->ncpp.nbeta; ib++)
                {
                    for(int ib2=0; ib2<atom->ncpp.nbeta; ib2++)
                    {
                        ofs << std::setw(10) << atom->ncpp.lll[ib] 
                            << " " << atom->ncpp.lll[ib2]
                            << " " << atom->ncpp.dion(ib,ib2)<<std::endl;
                    }
                }
                ofs << "</DIJ>" << std::endl;
    
                for(int i=0; i<atom->ncpp.nbeta; i++)
                {
                    ofs << "<PP_BETA>" << std::endl;
                    ofs << std::setw(10) << i << "\t" << "the index of projectors." <<std::endl;
                    ofs << std::setw(10) << atom->ncpp.lll[i] << "\t" << "the angular momentum." <<std::endl;
    
                    // mohan add
                    // only keep the nonzero part.
                    int cut_mesh = atom->ncpp.mesh; 
                    for(int j=atom->ncpp.mesh-1; j>=0; --j)
                    {
                        if( std::abs( atom->ncpp.betar(i,j) ) > 1.0e-10 )
                        {
                            cut_mesh = j; 
                            break;
                        }
                    }
                    if(cut_mesh %2 == 0) ++cut_mesh;
    
                    ofs << std::setw(10) << cut_mesh << "\t" << "the number of mesh points." << std::endl;
    
                    for(int j=0; j<cut_mesh; ++j)
                    {
                        ofs << std::setw(15) << atom->ncpp.r[j]
                            << std::setw(15) << atom->ncpp.betar(i, j)
                            << std::setw(15) << atom->ncpp.rab[j] << std::endl;
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

    for(int it=0; it<ntype; it++)
    {
        if(atoms[0].ncpp.xc_func !=atoms[it].ncpp.xc_func)
        {
            GlobalV::ofs_warning << "\n type " << atoms[0].label << " functional is " 
                << atoms[0].ncpp.xc_func;

            GlobalV::ofs_warning << "\n type " << atoms[it].label << " functional is " 
                << atoms[it].ncpp.xc_func << std::endl;

            ModuleBase::WARNING_QUIT("setup_cell","All DFT functional must consistent.");
        }
        if (atoms[it].ncpp.tvanp)
        {
            GlobalV::use_uspp = true;
        }
    }

    check_structure(GlobalV::MIN_DIST_COEF);

    // setup the total number of PAOs
    cal_natomwfc(ofs);

    // setup GlobalV::NLOCAL
    cal_nwfc(ofs);

    // Check whether the number of valence is minimum 
    if(GlobalV::MY_RANK==0)
    {
        int abtype = 0;
        for(int it=0; it<ntype; it++)
        {
            if(ModuleBase::MinZval.find(atoms[it].ncpp.psd) != ModuleBase::MinZval.end())
            {
                if(atoms[it].ncpp.zv > ModuleBase::MinZval.at(atoms[it].ncpp.psd))
                {
                    abtype += 1;
                    if(abtype == 1)
                    {
                        std::cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
                        ofs << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
                    }
                    std::cout<<" Warning: the number of valence electrons in pseudopotential > " << ModuleBase::MinZval.at(atoms[it].ncpp.psd);
                    std::cout<<" for " << atoms[it].ncpp.psd << ": " << ModuleBase::EleConfig.at(atoms[it].ncpp.psd) << std::endl;
                    ofs << " Warning: the number of valence electrons in pseudopotential > " << ModuleBase::MinZval.at(atoms[it].ncpp.psd);
                    ofs << " for " << atoms[it].ncpp.psd << ": " << ModuleBase::EleConfig.at(atoms[it].ncpp.psd) << std::endl;
                }
            }
        }
        if(abtype>0)
        {
            std::cout<< " Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient." << std::endl;
			std::cout<< " If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning."<<std::endl;
            std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"<<std::endl;
            ofs << " Pseudopotentials with additional electrons can yield (more) accurate outcomes, but may be less efficient."<<std::endl;
            ofs << " If you're confident that your chosen pseudopotential is appropriate, you can safely ignore this warning."<<std::endl;
            ofs << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
            ModuleBase::GlobalFunc::OUT(ofs,"");
        }
    }

    cal_meshx();

#ifdef __MPI
    Parallel_Common::bcast_int( meshx );
    Parallel_Common::bcast_int( natomwfc );
    Parallel_Common::bcast_int( lmax );
    Parallel_Common::bcast_int( lmax_ppwf );
#endif
}

//===========================================
// calculate the total number of local basis
// Target : nwfc, lmax,
// 			atoms[].stapos_wf
// 			GlobalV::NBANDS
//===========================================
void UnitCell::cal_nwfc(std::ofstream &log)
{
	ModuleBase::TITLE("UnitCell","cal_nwfc");
	assert(ntype>0);
	assert(nat>0);

	//===========================
	// (1) set iw2l, iw2n, iw2m
	//===========================
	for(int it=0; it<ntype; it++)
	{
		this->atoms[it].set_index();
	}

	//===========================
	// (2) set namax and nwmax
	//===========================
	this->namax = 0;
	this->nwmax = 0;
	for(int it=0; it<ntype; it++)
	{
		this->namax = std::max( atoms[it].na, namax );
		this->nwmax = std::max( atoms[it].nw, nwmax );
	}
	assert(namax>0);
// for tests
//		OUT(GlobalV::ofs_running,"max input atom number",namax);
//		OUT(GlobalV::ofs_running,"max wave function number",nwmax);	

	//===========================
	// (3) set nwfc and stapos_wf
	//===========================
	GlobalV::NLOCAL = 0;
	for(int it=0; it<ntype; it++)
	{
		atoms[it].stapos_wf = GlobalV::NLOCAL;
		const int nlocal_it = atoms[it].nw * atoms[it].na;
		if(GlobalV::NSPIN!=4) 
		{
			GlobalV::NLOCAL += nlocal_it;
		}
		else 
		{
			GlobalV::NLOCAL += nlocal_it * 2;//zhengdy-soc
		}

// for tests
//		OUT(GlobalV::ofs_running,ss1.str(),nlocal_it);
//		OUT(GlobalV::ofs_running,"start position of local orbitals",atoms[it].stapos_wf);
	}
	
	//OUT(GlobalV::ofs_running,"NLOCAL",GlobalV::NLOCAL);
	log << " " << std::setw(40) << "NLOCAL" << " = " << GlobalV::NLOCAL <<std::endl;
	//========================================================
	// (4) set index for itia2iat, itiaiw2iwt
	//========================================================

	// mohan add 2010-09-26
	assert(GlobalV::NLOCAL>0);
	delete[] iwt2iat;
	delete[] iwt2iw;
	this->iwt2iat = new int[GlobalV::NLOCAL];
	this->iwt2iw = new int[GlobalV::NLOCAL];

	this->itia2iat.create(ntype, namax);
	//this->itiaiw2iwt.create(ntype, namax, nwmax*GlobalV::NPOL);
	this->set_iat2iwt(GlobalV::NPOL);
	int iat=0;
	int iwt=0;
	for(int it = 0;it < ntype;it++)
	{
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			this->itia2iat(it, ia) = iat;
			//this->iat2ia[iat] = ia;
			for(int iw=0; iw<atoms[it].nw * GlobalV::NPOL; iw++)
			{
				//this->itiaiw2iwt(it, ia, iw) = iwt;
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
	for(int it=0; it<ntype; it++)
	{
		lmax = std::max( lmax, atoms[it].nwl);
		for(int l=0; l<atoms[it].nwl+1; l++)
		{
			nmax = std::max( nmax, atoms[it].l_nchi[l] );
		}

		int nchi = 0;
		for(int l=0; l<atoms[it].nwl+1; l++)
		{
			nchi += atoms[it].l_nchi[l];
		}
		this->nmax_total = std::max(nmax_total, nchi);
	}

	//=======================
	// (6) set lmax_ppwf
	//=======================
	this->lmax_ppwf = 0;
	for(int it=0; it< ntype; it++)
	{
		for(int ic=0; ic<atoms[it].ncpp.nchi; ic++)
		{
			if( lmax_ppwf < atoms[it].ncpp.lchi[ic] )
			{
				this->lmax_ppwf = atoms[it].ncpp.lchi[ic]; 
			}
		}
	}

	/*
	for(int it=0; it< ntype; it++)
	{
		std::cout << " label=" << it << " nbeta=" << atoms[it].nbeta << std::endl;
		for(int ic=0; ic<atoms[it].nbeta; ic++)
		{
			std::cout << " " << atoms[it].lll[ic] << std::endl;
		}
	}
	*/

//	OUT(GlobalV::ofs_running,"lmax between L(pseudopotential)",lmax_ppwf);

	//=====================
	// Use localized basis
	//=====================
	if(
		(GlobalV::BASIS_TYPE=="lcao")
	  ||(GlobalV::BASIS_TYPE=="lcao_in_pw")
	  ||(
            (GlobalV::BASIS_TYPE=="pw")
		  &&(GlobalV::psi_initializer)
		  &&(GlobalV::init_wfc.substr(0, 3)=="nao")
		  &&(GlobalV::ESOLVER_TYPE == "ksdft")
	    )
	 ) //xiaohui add 2013-09-02
	{
		ModuleBase::GlobalFunc::AUTO_SET("NBANDS",GlobalV::NBANDS);
	}
	else // plane wave basis
	{
		//if(winput::after_iter && winput::sph_proj)
		//{
		//	if(GlobalV::NBANDS < GlobalV::NLOCAL)
		//	{
		//		ModuleBase::WARNING_QUIT("cal_nwfc","NBANDS must > GlobalV::NLOCAL !");
		//	}
		//}
	}

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
	int iat=0;
	int iwt=0;
	for(int it = 0;it < this->ntype; it++)
	{
		for(int ia=0; ia<atoms[it].na; ia++)
		{
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
void UnitCell::cal_meshx()
{
	if(GlobalV::test_pseudo_cell) ModuleBase::TITLE("UnitCell","cal_meshx");
	this->meshx = 0;
	for (int it = 0;it < this->ntype;it++)
	{
		const int mesh = this->atoms[it].ncpp.msh;
		if (mesh > this->meshx)
		{
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
void UnitCell::cal_natomwfc(std::ofstream &log)
{
	if(GlobalV::test_pseudo_cell) ModuleBase::TITLE("UnitCell","cal_natomwfc");

	this->natomwfc = 0;
	for (int it = 0;it < ntype;it++)
	{
		//============================
		// Use pseudo-atomic orbitals
		//============================
		int tmp=0;
		for (int l = 0;l < atoms[it].ncpp.nchi;l++)
		{
			if (atoms[it].ncpp.oc[l] >= 0)
			{
				if(GlobalV::NSPIN==4)
				{
					if(atoms[it].ncpp.has_so)
					{
						tmp += 2 * atoms[it].ncpp.lchi[l];
						if(fabs(atoms[it].ncpp.jchi[l] - atoms[it].ncpp.lchi[l] - 0.5)<1e-6)
						tmp += 2 ;
					}
					else
					{
						tmp += 2 * (2 * atoms[it].ncpp.lchi[l] + 1);
					}
				}
				else
					tmp += 2 * atoms[it].ncpp.lchi[l] + 1;
			}
		}
		natomwfc += tmp * atoms[it].na;
	}
	ModuleBase::GlobalFunc::OUT(log,"initial pseudo atomic orbital number",natomwfc);
	return;
}



//LiuXh add a new function here,
//20180515
void UnitCell::setup_cell_after_vc(std::ofstream &log)
{
	ModuleBase::TITLE("UnitCell","setup_cell_after_vc");
    assert(lat0 > 0.0);
    this->omega = std::abs(latvec.Det()) * this->lat0 * lat0 * lat0;
    if(this->omega <= 0)
    {
        ModuleBase::WARNING_QUIT("setup_cell_after_vc", "omega <= 0 .");
    }
    else
    {
        log << std::endl;
        ModuleBase::GlobalFunc::OUT(log, "Volume (Bohr^3)", this->omega);
        ModuleBase::GlobalFunc::OUT(log, "Volume (A^3))", this->omega * pow(ModuleBase::BOHR_TO_A, 3));
    }

    lat0_angstrom = lat0 * 0.529177;
    tpiba  = ModuleBase::TWO_PI / lat0;
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
    this->G  = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    for(int it=0; it<ntype; it++)
    {
        Atom* atom = &atoms[it];
        for(int ia =0;ia< atom->na;ia++)
        {
            atom->tau[ia] = atom->taud[ia] * latvec;
        }
    }
    
#ifdef __MPI
    this->bcast_unitcell();
#endif

    log << std::endl;
    output::printM3(log,"Lattice vectors: (Cartesian coordinate: in unit of a_0)",latvec);
    output::printM3(log,"Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",G);

    return;
}

//check if any atom can be moved
bool UnitCell::if_atoms_can_move()const
{
	for(int it=0; it<this->ntype; it++)
    {
        Atom* atom = &atoms[it];
        for(int ia =0;ia< atom->na;ia++)
        {
            if(atom->mbl[ia].x||atom->mbl[ia].y||atom->mbl[ia].z) return 1;
		}
	}
	return 0;
}

//check if lattice vector can be changed
bool UnitCell::if_cell_can_change()const
{
	//need to be fixed next
	if(this->lc[0]||this->lc[1]||this->lc[2])
	{
		return 1;
	}
	return 0;
}

void UnitCell::setup(const std::string &latname_in,
	const int &ntype_in,
	const int &lmaxmax_in,
	const bool &init_vel_in,
	const std::string &fixed_axes_in)
{
	this->latName = latname_in;
	this->ntype = ntype_in;
	this->lmaxmax = lmaxmax_in;
	this->init_vel = init_vel_in;
	// pengfei Li add 2018-11-11
	if (fixed_axes_in == "None")
	{
		this->lc[0] = 1;
		this->lc[1] = 1;
		this->lc[2] = 1;
	}
	else if (fixed_axes_in == "volume")
	{
		this->lc[0] = 1;
		this->lc[1] = 1;
		this->lc[2] = 1;
		if(!GlobalV::relax_new)
		{
			ModuleBase::WARNING_QUIT("Input","there are bugs in the old implementation; set relax_new to be 1 for fixed_volume relaxation");
		}
	}
	else if (fixed_axes_in == "shape")
	{
		if(!GlobalV::relax_new)
		{
			ModuleBase::WARNING_QUIT("Input","set relax_new to be 1 for fixed_shape relaxation");
		}
		this->lc[0] = 1;
		this->lc[1] = 1;
		this->lc[2] = 1;
	}
	else if (fixed_axes_in == "a")
	{
		this->lc[0] = 0;
		this->lc[1] = 1;
		this->lc[2] = 1;
	}
	else if (fixed_axes_in == "b")
	{
		this->lc[0] = 1;
		this->lc[1] = 0;
		this->lc[2] = 1;
	}
	else if (fixed_axes_in == "c")
	{
		this->lc[0] = 1;
		this->lc[1] = 1;
		this->lc[2] = 0;
	}
	else if (fixed_axes_in == "ab")
	{
		this->lc[0] = 0;
		this->lc[1] = 0;
		this->lc[2] = 1;
	}
	else if (fixed_axes_in == "ac")
	{
		this->lc[0] = 0;
		this->lc[1] = 1;
		this->lc[2] = 0;
	}
	else if (fixed_axes_in == "bc")
	{
		this->lc[0] = 1;
		this->lc[1] = 0;
		this->lc[2] = 0;
	}
	else if (fixed_axes_in == "abc")
	{
		this->lc[0] = 0;
		this->lc[1] = 0;
		this->lc[2] = 0;
	}
	else
	{
		ModuleBase::WARNING_QUIT("Input", "fixed_axes should be None,volume,shape,a,b,c,ab,ac,bc or abc!");
	}
	return;
}

void UnitCell::check_structure(double factor)
{
	//First we calculate all bond length in the structure,
	//and compare with the covalent_bond_length,
	//if there has bond length is shorter than covalent_bond_length * factor,
	//we think this structure is unreasonable.
	const double warning_coef = 0.6;
	assert(ntype>0);
	std::stringstream errorlog;
	bool all_pass = true;
	bool no_warning = true;
	for (int it1 = 0;it1 < ntype; it1++)
	{ 
		std::string symbol1 = this->atoms[it1].ncpp.psd;
		double symbol1_covalent_radius;
		if (ModuleBase::CovalentRadius.find(symbol1) != ModuleBase::CovalentRadius.end())
		{
			symbol1_covalent_radius = ModuleBase::CovalentRadius.at(symbol1);
		}
		else
		{
			std::stringstream mess;
			mess << "Notice: symbol '" << symbol1 << "' is not an element symbol!!!! ";
			mess << "set the covalent radius to be 0." << std::endl;
			GlobalV::ofs_running << mess.str() ;
			std::cout << mess.str() ;
			symbol1_covalent_radius = 0.0;
		}
 
		for (int ia1 =0;ia1 <this->atoms[it1].na;ia1++)
		{
			double x1 = this->atoms[it1].taud[ia1].x;
			double y1 = this->atoms[it1].taud[ia1].y;
			double z1 = this->atoms[it1].taud[ia1].z;

			for(int it2=0;it2 <ntype;it2++)
			{
				std::string symbol2 = this->atoms[it2].ncpp.psd;
				double symbol2_covalent_radius;
				if (ModuleBase::CovalentRadius.find(symbol2) != ModuleBase::CovalentRadius.end())
				{
					symbol2_covalent_radius = ModuleBase::CovalentRadius.at(symbol2);
				}
				else
				{
					symbol2_covalent_radius = 0.0;
				}

				double covalent_length = (symbol1_covalent_radius + symbol2_covalent_radius) / ModuleBase::BOHR_TO_A;

				for (int ia2 =0;ia2 <this->atoms[it2].na;ia2++)
				{
					for(int a=-1; a<2; a++)
					{
						for(int b=-1;b<2;b++)
						{
							for (int c=-1;c<2;c++)
							{
								if (it1 > it2)
									continue;
								else if (it1==it2 && ia1 > ia2)
									continue;
								else if(it1==it2 && ia1==ia2 && a==0 && b==0 && c==0)
									continue;	

								double x2 = this->atoms[it2].taud[ia2].x + a;
								double y2 = this->atoms[it2].taud[ia2].y + b;
								double z2 = this->atoms[it2].taud[ia2].z + c;

								double bond_length = sqrt(pow((x2-x1)*this->a1.x + (y2-y1)*this->a2.x + (z2-z1)*this->a3.x,2) + 
														  pow((x2-x1)*this->a1.y + (y2-y1)*this->a2.y + (z2-z1)*this->a3.y,2) +
														  pow((x2-x1)*this->a1.z + (y2-y1)*this->a2.z + (z2-z1)*this->a3.z,2) ) * this->lat0;

								if (bond_length < covalent_length*factor || bond_length < covalent_length*warning_coef)
								{
									errorlog.setf(std::ios_base::fixed, std::ios_base::floatfield);
									errorlog << std::setw(3) << ia1+1 << "-th " << std::setw(3) << this->atoms[it1].label << ", "; 
									errorlog << std::setw(3) << ia2+1 << "-th " << std::setw(3) << this->atoms[it2].label;
									errorlog << " (cell:" << std::setw(2) << a << " " << std::setw(2) << b << " " << std::setw(2) << c << ")"; 
									errorlog << ", distance= " << std::setprecision(3) << bond_length << " Bohr (";
									errorlog << bond_length*ModuleBase::BOHR_TO_A << " Angstrom)" << std::endl;

									if (bond_length < covalent_length*factor)
									{
										all_pass = false;
									}
									else
									{
										no_warning = false;
									}
								}
							}//c
						}//b
					}//a
				}//ia2
			}//it2
		}//ia1	
	}//it1

	if (!all_pass || !no_warning)
	{
		std::stringstream mess;
		mess << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
		mess <<   "%%%%%% WARNING  WARNING  WARNING  WARNING  WARNING  %%%%%%" << std::endl;
		mess <<   "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
		mess << "!!! WARNING: Some atoms are too close!!!" << std::endl;
		mess << "!!! Please check the nearest-neighbor list in log file." << std::endl;
		mess <<   "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
		mess <<   "%%%%%% WARNING  WARNING  WARNING  WARNING  WARNING  %%%%%%" << std::endl;
		mess <<   "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

		GlobalV::ofs_running << mess.str() << mess.str() << mess.str() << errorlog.str();
		std::cout << mess.str() << mess.str() << mess.str() << std::endl;


		if (!all_pass)
		{
			mess.clear();
			mess.str("");
			mess << "If this structure is what you want, you can set 'min_dist_coef'" << std::endl;
			mess << "as a smaller value (the current value is " << factor << ") in INPUT file." << std::endl;
			GlobalV::ofs_running << mess.str();
			std::cout << mess.str();
			ModuleBase::WARNING_QUIT("Input", "The structure is unreasonable!");
		}
		
		
	}
}

void UnitCell::remake_cell()
{
	ModuleBase::TITLE("UnitCell","rmake_cell");

	//The idea is as follows: for each type of lattice, first calculate
	//from current latvec the lattice parameters, then use the parameters
	//to reconstruct latvec

	if(latName == "none")
	{
		ModuleBase::WARNING_QUIT("UnitCell","to use fixed_ibrav, latname must be provided");
	}
	else if(latName == "sc") //ibrav = 1
	{
		double celldm = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2));
		
		latvec.Zero();
		latvec.e11 = latvec.e22 = latvec.e33 = celldm;
	}
	else if(latName == "fcc") //ibrav = 2
	{
		double celldm = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2)) / std::sqrt(2.0);
		
		latvec.e11 =-celldm; latvec.e12 = 0.0;    latvec.e13 = celldm;
		latvec.e21 = 0.0;    latvec.e22 = celldm; latvec.e23 = celldm;
		latvec.e31 =-celldm; latvec.e32 = celldm; latvec.e33 = 0.0;
	}
	else if(latName == "bcc") //ibrav = 3
	{
		double celldm = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2)) / std::sqrt(3.0);
		
		latvec.e11 = celldm; latvec.e12 = celldm; latvec.e13 = celldm;
		latvec.e21 =-celldm; latvec.e22 = celldm; latvec.e23 = celldm;
		latvec.e31 =-celldm; latvec.e32 =-celldm; latvec.e33 = celldm;
	}
	else if(latName == "hexagonal") //ibrav = 4
	{
		double celldm1 = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2));
		double celldm3 = std::sqrt(pow(latvec.e31,2)+pow(latvec.e32,2)+pow(latvec.e33,2));
		double e22 = sqrt(3.0) / 2.0;

		latvec.e11 = celldm1;     latvec.e12 = 0.0;           latvec.e13 = 0.0;
		latvec.e21 =-0.5*celldm1; latvec.e22 = celldm1 * e22; latvec.e23 = 0.0;
		latvec.e31 = 0.0;         latvec.e32 = 0.0;	          latvec.e33 = celldm3;
	}
	else if(latName == "trigonal") //ibrav = 5
	{
		double celldm1  = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2));
		double celldm2  = std::sqrt(pow(latvec.e21,2)+pow(latvec.e22,2)+pow(latvec.e23,2));
		double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22 + latvec.e13 * latvec.e23);
		double cos12 = celldm12 / celldm1 / celldm2;

		if(cos12 <= -0.5 || cos12 >= 1.0)
		{
			ModuleBase::WARNING_QUIT("unitcell","wrong cos12!");
		}
		double t1 = sqrt(1.0 + 2.0*cos12);
		double t2 = sqrt(1.0 - cos12);

		double e11 =  celldm1 * t2 / sqrt(2.0);
		double e12 = -celldm1 * t2 / sqrt(6.0);
		double e13 =  celldm1 * t1 / sqrt(3.0);
		double e22 =  celldm1 * sqrt(2.0) * t2 / sqrt(3.0);
	
		latvec.e11 = e11; latvec.e12 = e12; latvec.e13 = e13;
		latvec.e21 = 0.0; latvec.e22 = e22;	latvec.e23 = e13;
		latvec.e31 =-e11; latvec.e32 = e12;	latvec.e33 = e13;		
	}
	else if(latName == "st") //ibrav = 6
	{
		double celldm1 = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2));
		double celldm3 = std::sqrt(pow(latvec.e31,2)+pow(latvec.e32,2)+pow(latvec.e33,2));
		latvec.e11 = celldm1; latvec.e12 = 0.0;     latvec.e13 = 0.0;
		latvec.e21 = 0.0;     latvec.e22 = celldm1; latvec.e23 = 0.0;
		latvec.e31 = 0.0;      latvec.e32 = 0.0;	latvec.e33 = celldm3;
	}
	else if(latName == "bct") //ibrav = 7
	{
		double celldm1 = std::abs(latvec.e11);
		double celldm2 = std::abs(latvec.e13);
			
		latvec.e11 = celldm1; latvec.e12 =-celldm1; latvec.e13 = celldm2;
		latvec.e21 = celldm1; latvec.e22 = celldm1; latvec.e23 = celldm2;
		latvec.e31 =-celldm1; latvec.e32 =-celldm1;	latvec.e33 = celldm2;	
	}
	else if(latName == "so") //ibrav = 8
	{
		double celldm1 = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2));
		double celldm2 = std::sqrt(pow(latvec.e21,2)+pow(latvec.e22,2)+pow(latvec.e23,2));
		double celldm3 = std::sqrt(pow(latvec.e31,2)+pow(latvec.e32,2)+pow(latvec.e33,2));
		
		latvec.e11 = celldm1; latvec.e12 = 0.0;     latvec.e13 = 0.0;
		latvec.e21 = 0.0;     latvec.e22 = celldm2;	latvec.e23 = 0.0;
		latvec.e31 = 0.0;     latvec.e32 = 0.0;     latvec.e33 = celldm3;
	}
	else if(latName == "baco") //ibrav = 9
	{
		double celldm1 = std::abs(latvec.e11);
		double celldm2 = std::abs(latvec.e22);
		double celldm3 = std::abs(latvec.e33);

		latvec.e11 = celldm1; latvec.e12 = celldm2; latvec.e13 = 0.0;
		latvec.e21 =-celldm1; latvec.e22 = celldm2;	latvec.e23 = 0.0;
		latvec.e31 = 0.0;     latvec.e32 = 0.0;   	latvec.e33 = celldm3;
	}
	else if(latName == "fco") //ibrav = 10
	{
		double celldm1 = std::abs(latvec.e11);
		double celldm2 = std::abs(latvec.e22);
		double celldm3 = std::abs(latvec.e33);

		latvec.e11 = celldm1; latvec.e12 = 0.0;     latvec.e13 = celldm3;
		latvec.e21 = celldm1; latvec.e22 = celldm2;	latvec.e23 = 0.0;
		latvec.e31 = 0.0;     latvec.e32 = celldm2;	latvec.e33 = celldm3;
	}
	else if(latName == "bco") //ibrav = 11
	{
		double celldm1 = std::abs(latvec.e11);
		double celldm2 = std::abs(latvec.e12);
		double celldm3 = std::abs(latvec.e13);

		latvec.e11 = celldm1; latvec.e12 = celldm2; latvec.e13 = celldm3;
		latvec.e21 =-celldm1; latvec.e22 = celldm2;	latvec.e23 = celldm3;
		latvec.e31 =-celldm1; latvec.e32 =-celldm2;	latvec.e33 = celldm3;		
	}
	else if(latName == "sm") //ibrav = 12
	{
		double celldm1 = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2));
		double celldm2 = std::sqrt(pow(latvec.e21,2)+pow(latvec.e22,2)+pow(latvec.e23,2));
		double celldm3 = std::sqrt(pow(latvec.e31,2)+pow(latvec.e32,2)+pow(latvec.e33,2));
		double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22 + latvec.e13 * latvec.e23);
		double cos12 = celldm12 / celldm1 / celldm2;

		double e21 = celldm2 * cos12;
		double e22 = celldm2 * std::sqrt(1.0 - cos12 * cos12);

		latvec.e11 = celldm1; latvec.e12 = 0.0; latvec.e13 = 0.0;
		latvec.e21 = e21;     latvec.e22 = e22;	latvec.e23 = 0.0;
		latvec.e31 = 0.0;     latvec.e32 = 0.0;	latvec.e33 = celldm3;
	}
	else if (latName == "bacm") //ibrav = 13
	{
		double celldm1 = std::abs(latvec.e11);
		double celldm2 = std::sqrt(pow(latvec.e21,2)+pow(latvec.e22,2)+pow(latvec.e23,2));
		double celldm3 = std::abs(latvec.e13);

		double cos12 = latvec.e21 / celldm2;
		if(cos12 >= 1.0)
		{
			ModuleBase::WARNING_QUIT("unitcell","wrong cos12!");
		}

		double e21 = celldm2 * cos12;
		double e22 = celldm2 * std::sqrt(1.0 - cos12 * cos12);

		latvec.e11 = celldm1; latvec.e12 = 0.0; latvec.e13 =-celldm3;
		latvec.e21 = e21;     latvec.e22 = e22;	latvec.e23 = 0.0;
		latvec.e31 = celldm1; latvec.e32 = 0.0;	latvec.e33 = celldm3;		
	}
	else if(latName == "triclinic") //ibrav = 14
	{
		double celldm1 = std::sqrt(pow(latvec.e11,2)+pow(latvec.e12,2)+pow(latvec.e13,2));
		double celldm2 = std::sqrt(pow(latvec.e21,2)+pow(latvec.e22,2)+pow(latvec.e23,2));
		double celldm3 = std::sqrt(pow(latvec.e31,2)+pow(latvec.e32,2)+pow(latvec.e33,2));
		double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22 + latvec.e13 * latvec.e23);
		double cos12 = celldm12 / celldm1 / celldm2;
		double celldm13 = (latvec.e11 * latvec.e31 + latvec.e12 * latvec.e32 + latvec.e13 * latvec.e33);
		double cos13 = celldm13 / celldm1 / celldm3;
		double celldm23 = (latvec.e21 * latvec.e31 + latvec.e22 * latvec.e32 + latvec.e23 * latvec.e33);
		double cos23 = celldm23 / celldm2 / celldm3;

		double sin12 = std::sqrt(1.0 - cos12 * cos12);
		if(cos12 >= 1.0)
		{
			ModuleBase::WARNING_QUIT("unitcell","wrong cos12!");
		}

		latvec.e11 = celldm1; latvec.e12 = 0.0; latvec.e13 = 0.0;
		latvec.e21 = celldm2 * cos12;
		latvec.e22 = celldm2 * sin12;
		latvec.e23 = 0.0;
		latvec.e31 = celldm3 * cos13;
		latvec.e32 = celldm3 * (cos23 - cos13*cos12) / sin12;
		double term = 1.0 + 2.0 * cos12*cos13*cos23 - cos12*cos12 - cos13*cos13 - cos23*cos23;
		term = sqrt(term)/sin12;
		latvec.e33 = celldm3 * term;
	}
	else{ 
		std::cout << "latname is : " << latName << std::endl;
		ModuleBase::WARNING_QUIT("UnitCell::read_atom_species","latname not supported!");
	}
}

void UnitCell::cal_nelec(double& nelec)
{
    ModuleBase::TITLE("UnitCell", "cal_nelec");
    GlobalV::ofs_running << "\n SETUP THE ELECTRONS NUMBER" << std::endl;

    if (nelec == 0)
    {
		if(GlobalV::use_paw)
		{
#ifdef USE_PAW
			for(int it = 0; it < this->ntype; it ++)
			{
				std::stringstream ss1, ss2;
				ss1 << " electron number of element " << GlobalC::paw_cell.get_zat(it) << std::endl;
				const int nelec_it = GlobalC::paw_cell.get_val(it) * this->atoms[it].na;
				nelec += nelec_it;
				ss2 << "total electron number of element " << GlobalC::paw_cell.get_zat(it);

				ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss1.str(), GlobalC::paw_cell.get_val(it));
				ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss2.str(), nelec_it);				
			}
#endif
		}
		else
		{
			for (int it = 0; it < this->ntype; it++)
			{
				std::stringstream ss1, ss2;
				ss1 << "electron number of element " << this->atoms[it].label;
				const int nelec_it = this->atoms[it].ncpp.zv * this->atoms[it].na;
				nelec += nelec_it;
				ss2 << "total electron number of element " << this->atoms[it].label;

				ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss1.str(), this->atoms[it].ncpp.zv);
				ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, ss2.str(), nelec_it);
			}
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "AUTOSET number of electrons: ", nelec);
		}
    }
	if (GlobalV::nelec_delta != 0)
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec_delta is NOT zero, please make sure you know what you are doing! nelec_delta: ", GlobalV::nelec_delta);
		nelec += GlobalV::nelec_delta;
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nelec now: ", nelec);
	}
    return;
}

void UnitCell::compare_atom_labels(std::string label1, std::string label2)
{
    if (label1 != label2) //'!( "Ag" == "Ag" || "47" == "47" || "Silver" == Silver" )'
    {	
        atom_in ai;
        if (!(std::to_string(ai.atom_Z[label1]) == label2 ||   // '!( "Ag" == "47" )'
			  ai.atom_symbol[label1] == label2 ||              // '!( "Ag" == "Silver" )'
			  label1 == std::to_string(ai.atom_Z[label2]) ||   // '!( "47" == "Ag" )'
		      label1 == std::to_string(ai.symbol_Z[label2]) || // '!( "47" == "Silver" )'
			  label1 == ai.atom_symbol[label2] ||              // '!( "Silver" == "Ag" )'
			  std::to_string(ai.symbol_Z[label1]) == label2 )) // '!( "Silver" == "47" )'
	    {		
	    	std::string stru_label = "";
            std::string psuedo_label = "";
            for (int ip = 0; ip < label1.length(); ip++)
            {
                if (!(isdigit(label1[ip]) || label1[ip]=='_'))
                {
                    stru_label += label1[ip];
                }
	    		else
	    		{
	    			break;
	    		}
            }
	    	stru_label[0] = toupper(stru_label[0]);
    
	    	for (int ip = 0; ip < label2.length(); ip++)
            {
                if (!(isdigit(label2[ip]) || label2[ip]=='_'))
                {
                    psuedo_label += label2[ip];
	    		}
	    		else
	    		{
	    			break;
	    		}
            }
	    	psuedo_label[0] = toupper(psuedo_label[0]);
    
            if (!(stru_label == psuedo_label || //' !("Ag1" == "ag_locpsp" || "47" == "47" || "Silver" == Silver" )'
			      std::to_string(ai.atom_Z[stru_label]) == psuedo_label ||   // ' !("Ag1" == "47" )'
			      ai.atom_symbol[stru_label] == psuedo_label ||              // ' !("Ag1" == "Silver")'
			      stru_label == std::to_string(ai.atom_Z[psuedo_label]) ||  // ' !("47" == "Ag1" )'
		          stru_label == std::to_string(ai.symbol_Z[psuedo_label]) || // ' !("47" == "Silver1" )'
			      stru_label == ai.atom_symbol[psuedo_label] ||              // ' !("Silver1" == "Ag" )'
			      std::to_string(ai.symbol_Z[stru_label]) == psuedo_label )) // ' !("Silver1" == "47" )'
            
				
			{	
				std::string atom_label_in_orbtial = "atom label in orbital file ";
				std::string mismatch_with_pseudo = " mismatch with pseudo file of ";
                ModuleBase::WARNING_QUIT("UnitCell::read_pseudo", atom_label_in_orbtial + label1 + mismatch_with_pseudo +label2);
            }
	    }
	}
}