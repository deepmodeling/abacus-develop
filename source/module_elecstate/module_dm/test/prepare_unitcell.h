#ifndef PREPARE_UNITCELL_H
#define PREPARE_UNITCELL_H
#include<map>
#include<string>
#include "module_base/mathzone.h"

class UcellTestPrepare
{
public:
	UcellTestPrepare()=default;
	UcellTestPrepare(std::string latname_in,
		int lmaxmax_in,
		bool init_vel_in,
		bool selective_dynamics_in,
		bool relax_new_in,
		std::string fixed_axes_in,
		double lat0_in,
		std::valarray<double> latvec_in,
		std::vector<std::string> elements_in,
		std::vector<std::string> pp_files_in,
		std::vector<std::string> pp_types_in,
		std::vector<std::string> orb_files_in,
		std::valarray<int> natom_in,
		std::vector<double> atomic_mass_in,
		std::string coor_type_in,
		std::valarray<double> coordinates_in);
	UcellTestPrepare(std::string latname_in,
		int lmaxmax_in,
		bool init_vel_in,
		bool selective_dynamics_in,
		bool relax_new_in,
		std::string fixed_axes_in,
		double lat0_in,
		std::valarray<double> latvec_in,
		std::vector<std::string> elements_in,
		std::vector<std::string> pp_files_in,
		std::vector<std::string> pp_types_in,
		std::vector<std::string> orb_files_in,
		std::valarray<int> natom_in,
		std::vector<double> atomic_mass_in,
		std::string coor_type_in,
		std::valarray<double> coordinates_in,
		std::valarray<double> mbl_in,
		std::valarray<double> velocity_in);
	UcellTestPrepare(const UcellTestPrepare &utp);

	std::string latname;
	int lmaxmax;
	bool init_vel;
	bool selective_dynamics;
	bool relax_new;
	std::string fixed_axes;
	double lat0;
	std::valarray<double> latvec;
	std::vector<std::string> elements;
	std::vector<std::string> pp_files;
	std::vector<std::string> pp_types;
	std::vector<std::string> orb_files;
	std::valarray<int> natom;
	std::vector<double> atomic_mass;
	std::string coor_type;
	std::valarray<double> coordinates;
	std::valarray<double> mbl;
	std::valarray<double> velocity;
	// ntype
	int ntype;
	int atomic_index;

    UnitCell* SetUcellInfo(const std::vector<int>& nw, int& nlocal)
    {
		//basic info
		this->ntype = this->elements.size();
        static UnitCell ucell;
        ucell.setup(this->latname, this->ntype, this->lmaxmax, this->init_vel, this->fixed_axes);
        delete[] ucell.atom_label;
        delete[] ucell.atom_mass;
        delete[] ucell.pseudo_fn;
        delete[] ucell.pseudo_type;
        delete[] ucell.orbital_fn;
        delete[] ucell.magnet.start_magnetization; // mag set here
        ucell.atom_label = new std::string[ucell.ntype];
        ucell.atom_mass = new double[ucell.ntype];
        ucell.pseudo_fn = new std::string[ucell.ntype];
        ucell.pseudo_type = new std::string[ucell.ntype];
        ucell.orbital_fn = new std::string[ucell.ntype];
        ucell.magnet.start_magnetization = new double[ucell.ntype]; // mag set here
        ucell.magnet.ux_[0] = 0.0;                                  // ux_ set here
        ucell.magnet.ux_[1] = 0.0;
        ucell.magnet.ux_[2] = 0.0;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            ucell.atom_label[it] = this->elements[it];
            ucell.atom_mass[it] = this->atomic_mass[it];
            ucell.pseudo_fn[it] = this->pp_files[it];
            ucell.pseudo_type[it] = this->pp_types[it];
            ucell.orbital_fn[it] = this->orb_files[it];
            ucell.magnet.start_magnetization[it] = 0.0; // mag set here
        }
        // lattice info
        ucell.lat0 = this->lat0;
        ucell.lat0_angstrom = ucell.lat0 * 0.529177;
        ucell.tpiba = ModuleBase::TWO_PI / ucell.lat0;
        ucell.tpiba2 = ucell.tpiba * ucell.tpiba;
        ucell.latvec.e11 = this->latvec[0];
        ucell.latvec.e12 = this->latvec[1];
        ucell.latvec.e13 = this->latvec[2];
        ucell.latvec.e21 = this->latvec[3];
        ucell.latvec.e22 = this->latvec[4];
        ucell.latvec.e23 = this->latvec[5];
        ucell.latvec.e31 = this->latvec[6];
        ucell.latvec.e32 = this->latvec[7];
        ucell.latvec.e33 = this->latvec[8];
        ucell.a1.x = ucell.latvec.e11;
        ucell.a1.y = ucell.latvec.e12;
        ucell.a1.z = ucell.latvec.e13;
        ucell.a2.x = ucell.latvec.e21;
        ucell.a2.y = ucell.latvec.e22;
        ucell.a2.z = ucell.latvec.e23;
        ucell.a3.x = ucell.latvec.e31;
        ucell.a3.y = ucell.latvec.e32;
        ucell.a3.z = ucell.latvec.e33;
        ucell.GT = ucell.latvec.Inverse();
        ucell.G = ucell.GT.Transpose();
        ucell.GGT = ucell.G * ucell.GT;
        ucell.invGGT = ucell.GGT.Inverse();
        ucell.omega = std::abs(ucell.latvec.Det()) * (ucell.lat0) * (ucell.lat0) * (ucell.lat0);
        // atomic info
        ucell.Coordinate = this->coor_type;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.set_atom_flag = true;
        this->atomic_index = 0;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            ucell.atoms[it].label = this->elements[it];
            /*
            ucell.atoms[it].nw = 0;
            ucell.atoms[it].nwl = 2;
            delete[] ucell.atoms[it].l_nchi;
            ucell.atoms[it].l_nchi = new int[ ucell.atoms[it].nwl+1];
            for(int L=0; L<ucell.atoms[it].nwl+1; L++)
            {
                ucell.atoms[it].l_nchi[L] = 1;
                ucell.atoms[it].nw += (2*L + 1) * ucell.atoms[it].l_nchi[L];
            }
            */
            ucell.atoms[it].na = this->natom[it];
            // coordinates and related physical quantities
            delete[] ucell.atoms[it].tau;
            delete[] ucell.atoms[it].dis;
            delete[] ucell.atoms[it].taud;
            delete[] ucell.atoms[it].vel;
            delete[] ucell.atoms[it].mag;
            delete[] ucell.atoms[it].angle1;
            delete[] ucell.atoms[it].angle2;
            delete[] ucell.atoms[it].m_loc_;
            delete[] ucell.atoms[it].mbl;
            ucell.atoms[it].tau = new ModuleBase::Vector3<double>[ucell.atoms[it].na];
            ucell.atoms[it].dis = new ModuleBase::Vector3<double>[ucell.atoms[it].na];
            ucell.atoms[it].taud = new ModuleBase::Vector3<double>[ucell.atoms[it].na];
            ucell.atoms[it].vel = new ModuleBase::Vector3<double>[ucell.atoms[it].na];
            ucell.atoms[it].mag = new double[ucell.atoms[it].na];
            ucell.atoms[it].angle1 = new double[ucell.atoms[it].na];
            ucell.atoms[it].angle2 = new double[ucell.atoms[it].na];
            ucell.atoms[it].m_loc_ = new ModuleBase::Vector3<double>[ucell.atoms[it].na];
            ucell.atoms[it].mbl = new ModuleBase::Vector3<int>[ucell.atoms[it].na];
            ucell.atoms[it].mass = ucell.atom_mass[it]; // mass set here
            for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
            {
                if (ucell.Coordinate == "Direct")
                {
                    ucell.atoms[it].taud[ia].x = this->coordinates[this->atomic_index * 3 + 0];
                    ucell.atoms[it].taud[ia].y = this->coordinates[this->atomic_index * 3 + 1];
                    ucell.atoms[it].taud[ia].z = this->coordinates[this->atomic_index * 3 + 2];
                    ucell.atoms[it].tau[ia] = ucell.atoms[it].taud[ia] * ucell.latvec;
                }
                else if (ucell.Coordinate == "Cartesian")
                {
                    ucell.atoms[it].tau[ia].x = this->coordinates[this->atomic_index * 3 + 0];
                    ucell.atoms[it].tau[ia].y = this->coordinates[this->atomic_index * 3 + 1];
                    ucell.atoms[it].tau[ia].z = this->coordinates[this->atomic_index * 3 + 2];
                    ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[it].tau[ia].x,
                                                              ucell.atoms[it].tau[ia].y,
                                                              ucell.atoms[it].tau[ia].z,
                                                              ucell.latvec.e11,
                                                              ucell.latvec.e12,
                                                              ucell.latvec.e13,
                                                              ucell.latvec.e21,
                                                              ucell.latvec.e22,
                                                              ucell.latvec.e23,
                                                              ucell.latvec.e31,
                                                              ucell.latvec.e32,
                                                              ucell.latvec.e33,
                                                              ucell.atoms[it].taud[ia].x,
                                                              ucell.atoms[it].taud[ia].y,
                                                              ucell.atoms[it].taud[ia].z);
                }
                ucell.atoms[it].dis[ia].set(0, 0, 0);
                if (this->init_vel)
                {
                    ucell.atoms[it].vel[ia].x = this->velocity[this->atomic_index * 3 + 0];
                    ucell.atoms[it].vel[ia].y = this->velocity[this->atomic_index * 3 + 1];
                    ucell.atoms[it].vel[ia].z = this->velocity[this->atomic_index * 3 + 2];
                }
                else
                {
                    ucell.atoms[it].vel[ia].set(0, 0, 0);
                }
                ucell.atoms[it].m_loc_[ia].set(0, 0, 0);
                ucell.atoms[it].angle1[ia] = 0;
                ucell.atoms[it].angle2[ia] = 0;
                if (this->selective_dynamics)
                {
                    ucell.atoms[it].mbl[ia].x = this->mbl[this->atomic_index * 3 + 0];
                    ucell.atoms[it].mbl[ia].y = this->mbl[this->atomic_index * 3 + 1];
                    ucell.atoms[it].mbl[ia].z = this->mbl[this->atomic_index * 3 + 2];
                }
                else
                {
                    ucell.atoms[it].mbl[ia] = {1, 1, 1};
                }
                ++(this->atomic_index);
            }
        }
        ucell.nat = this->natom.sum();
        // set_nw
        assert(nw.size() == ucell.ntype);
        for (int it = 0; it < ucell.ntype; ++it)
        {
            ucell.atoms[it].nw = nw[it];
        }
        // cal_nloc
        for (int it = 0; it < ucell.ntype; ++it)
        {
            nlocal += ucell.atoms[it].na * ucell.atoms[it].nw;
        }
        // cal_namax
        int namax = 0;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            namax = std::max(namax, ucell.atoms[it].na);
        }
        ucell.namax = namax;
        // cal_index
        assert(nlocal > 0);
        delete[] ucell.iwt2iat;
        delete[] ucell.iwt2iw;
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
        ucell.iwt2iat = new int[nlocal];
        ucell.iwt2iw = new int[nlocal];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat]; // set_iat2itia
        ucell.itia2iat.create(ucell.ntype, ucell.namax);
        ucell.iat2iwt.resize(ucell.nat);
        int iat = 0;
        int iwt = 0;
        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                ucell.iat2it[iat] = it;
                ucell.iat2ia[iat] = ia;
                ucell.itia2iat(it, ia) = iat;
                ucell.iat2iwt[iat] = iwt;
                for (int iw = 0; iw < ucell.atoms[it].nw; iw++)
                {
                    ucell.iwt2iat[iwt] = iat;
                    ucell.iwt2iw[iwt] = iw;
                    ++iwt;
                }
                ++iat;
            }
        }
        return &ucell;
    }
};

UcellTestPrepare::UcellTestPrepare(std::string latname_in,
		int lmaxmax_in,
		bool init_vel_in,
		bool selective_dynamics_in,
		bool relax_new_in,
		std::string fixed_axes_in,
		double lat0_in,
		std::valarray<double> latvec_in,
		std::vector<std::string> elements_in,
		std::vector<std::string> pp_files_in,
		std::vector<std::string> pp_types_in,
		std::vector<std::string> orb_files_in,
		std::valarray<int> natom_in,
		std::vector<double> atomic_mass_in,
		std::string coor_type_in,
		std::valarray<double> coordinates_in):
	latname(latname_in),
	lmaxmax(lmaxmax_in),
	init_vel(init_vel_in),
	selective_dynamics(selective_dynamics_in),
	relax_new(relax_new_in),
	fixed_axes(fixed_axes_in),
	lat0(lat0_in),
	latvec(latvec_in),
	elements(elements_in),
	pp_files(pp_files_in),
	pp_types(pp_types_in),
	orb_files(orb_files_in),
	natom(natom_in),
	atomic_mass(atomic_mass_in),
	coor_type(coor_type_in),
	coordinates(coordinates_in)
{
	mbl = {0};
	velocity = {0};
}

UcellTestPrepare::UcellTestPrepare(std::string latname_in,
		int lmaxmax_in,
		bool init_vel_in,
		bool selective_dynamics_in,
		bool relax_new_in,
		std::string fixed_axes_in,
		double lat0_in,
		std::valarray<double> latvec_in,
		std::vector<std::string> elements_in,
		std::vector<std::string> pp_files_in,
		std::vector<std::string> pp_types_in,
		std::vector<std::string> orb_files_in,
		std::valarray<int> natom_in,
		std::vector<double> atomic_mass_in,
		std::string coor_type_in,
		std::valarray<double> coordinates_in,
		std::valarray<double> mbl_in,
		std::valarray<double> velocity_in):
	latname(latname_in),
	lmaxmax(lmaxmax_in),
	init_vel(init_vel_in),
	selective_dynamics(selective_dynamics_in),
	relax_new(relax_new_in),
	fixed_axes(fixed_axes_in),
	lat0(lat0_in),
	latvec(latvec_in),
	elements(elements_in),
	pp_files(pp_files_in),
	pp_types(pp_types_in),
	orb_files(orb_files_in),
	natom(natom_in),
	atomic_mass(atomic_mass_in),
	coor_type(coor_type_in),
	coordinates(coordinates_in),
	mbl(mbl_in),
	velocity(velocity_in) // velocity assume the existence of mbl in print_stru_file()
{}

UcellTestPrepare::UcellTestPrepare(const UcellTestPrepare &utp):
	latname(utp.latname),
	lmaxmax(utp.lmaxmax),
	init_vel(utp.init_vel),
	selective_dynamics(utp.selective_dynamics),
	relax_new(utp.relax_new),
	fixed_axes(utp.fixed_axes),
	lat0(utp.lat0),
	latvec(utp.latvec),
	elements(utp.elements),
	pp_files(utp.pp_files),
	pp_types(utp.pp_types),
	orb_files(utp.orb_files),
	natom(utp.natom),
	atomic_mass(utp.atomic_mass),
	coor_type(utp.coor_type),
	coordinates(utp.coordinates),
	mbl(utp.mbl),
	velocity(utp.velocity) // velocity assume the existence of mbl in print_stru_file()
{}

std::map<std::string,UcellTestPrepare> UcellTestLib
{
	{"Si", UcellTestPrepare(
				"fcc",		//latname
				2,		//lmaxmax
				true,		//init_vel
				true,		//selective_dyanmics
				true,		//relax_new
				"volume",	//fixed_axes
				10.2, //lat0
				{0.5,0.5,0.0,	//latvec
				 0.5,0.0,0.5,
				 0.0,0.5,0.5},
				{"Si"},	//elements
				{"Si.upf"},	//upf file
				{"upf201"},	//upf types
				{"Si.orb"},	//orb file
				{2},		//number of each elements
				{28.0},	//atomic mass
				"Cartesian",	//coordination type
				{0.0,0.0,0.0,	//atomic coordinates
				 0.25,0.25,0.25})}
};
#endif
