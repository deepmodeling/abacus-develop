#include <cstdlib>
#ifdef __MPI
#include "mpi.h"
#endif

//#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/constants.h"
#include "unitcell.h"
using namespace std;

UnitCell::UnitCell()
{
    if (GlobalV::test_unitcell) ModuleBase::TITLE("unitcell","Constructor");
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

    itia2iat.create(1,1);
	lc = new int[3];
    itiaiw2iwt.create(1,1,1);

    latvec = ModuleBase::Matrix3();
    latvec_supercell = ModuleBase::Matrix3();
    G = ModuleBase::Matrix3();
    GT = ModuleBase::Matrix3();
    GGT = ModuleBase::Matrix3();
    invGGT = ModuleBase::Matrix3();

    tpiba = 0.0;
    tpiba2 = 0.0;
    omega = 0.0;

    atom_label = new string[1];
    atom_mass = new double[1];
    pseudo_fn = new string[1];
}

UnitCell::~UnitCell()
{
    delete[] atom_label;
    delete[] atom_mass;
    delete[] pseudo_fn;
    delete[] iat2it;
    delete[] iat2ia;
	delete[] iwt2iat;
	delete[] iwt2iw;
	delete[] lc;
}

#include "../src_parallel/parallel_common.h"
#ifdef __MPI
void UnitCell::bcast_unitcell(void)
{
    if (GlobalV::test_unitcell)ModuleBase::TITLE("UnitCell","bcast_unitcell");
    Parallel_Common::bcast_string( Coordinate );
    Parallel_Common::bcast_int( nat );

    Parallel_Common::bcast_double( lat0 );
    Parallel_Common::bcast_double( lat0_angstrom );
	Parallel_Common::bcast_double( tpiba );
	Parallel_Common::bcast_double( tpiba2 );

	// distribute lattice vectors.
    Parallel_Common::bcast_double( latvec.e11 );
    Parallel_Common::bcast_double( latvec.e12 );
    Parallel_Common::bcast_double( latvec.e13 );
    Parallel_Common::bcast_double( latvec.e21 );
    Parallel_Common::bcast_double( latvec.e22 );
    Parallel_Common::bcast_double( latvec.e23 );
    Parallel_Common::bcast_double( latvec.e31 );
    Parallel_Common::bcast_double( latvec.e32 );
    Parallel_Common::bcast_double( latvec.e33 );

    Parallel_Common::bcast_int( lc[0] );
    Parallel_Common::bcast_int( lc[1] );
    Parallel_Common::bcast_int( lc[2] );
	
	// distribute lattice vectors.
    Parallel_Common::bcast_double( a1.x );
    Parallel_Common::bcast_double( a1.y );
    Parallel_Common::bcast_double( a1.z );
    Parallel_Common::bcast_double( a2.x );
    Parallel_Common::bcast_double( a2.y );
    Parallel_Common::bcast_double( a2.z );
    Parallel_Common::bcast_double( a3.x );
    Parallel_Common::bcast_double( a3.y );
    Parallel_Common::bcast_double( a3.z );

	// distribute latcenter
    Parallel_Common::bcast_double( latcenter.x );
    Parallel_Common::bcast_double( latcenter.y );
    Parallel_Common::bcast_double( latcenter.z );

	// distribute superlattice vectors.
    Parallel_Common::bcast_double( latvec_supercell.e11 );
    Parallel_Common::bcast_double( latvec_supercell.e12 );
    Parallel_Common::bcast_double( latvec_supercell.e13 );
    Parallel_Common::bcast_double( latvec_supercell.e21 );
    Parallel_Common::bcast_double( latvec_supercell.e22 );
    Parallel_Common::bcast_double( latvec_supercell.e23 );
    Parallel_Common::bcast_double( latvec_supercell.e31 );
    Parallel_Common::bcast_double( latvec_supercell.e32 );
    Parallel_Common::bcast_double( latvec_supercell.e33 );

#ifndef __CMD
    Parallel_Common::bcast_double( magnet.start_magnetization, ntype );

    if(GlobalV::NSPIN==4)
    {
         Parallel_Common::bcast_double( magnet.ux_[0] );
         Parallel_Common::bcast_double( magnet.ux_[1] );
         Parallel_Common::bcast_double( magnet.ux_[2] );
    }
#endif

    for (int i=0;i<ntype;i++)
    {
        atoms[i].bcast_atom(); // init tau and mbl array
    }
    return;
}

void UnitCell::bcast_unitcell2(void)
{
    for (int i=0;i<ntype;i++)
    {
        atoms[i].bcast_atom2();
    }
    return;
}
#endif

void UnitCell::print_cell(std::ofstream &ofs)const
{
    if (GlobalV::test_unitcell) ModuleBase::TITLE("UnitCell","print_cell");

    ModuleBase::GlobalFunc::OUT(ofs,"print_unitcell()");

    ModuleBase::GlobalFunc::OUT(ofs,"latName",latName);
    ModuleBase::GlobalFunc::OUT(ofs,"ntype",ntype);
    ModuleBase::GlobalFunc::OUT(ofs,"nat",nat);
    ModuleBase::GlobalFunc::OUT(ofs,"lat0",lat0);
    ModuleBase::GlobalFunc::OUT(ofs,"lat0_angstrom",lat0_angstrom);
    ModuleBase::GlobalFunc::OUT(ofs,"tpiba",tpiba);
    ModuleBase::GlobalFunc::OUT(ofs,"omega",omega);

    output::printM3(ofs,"Lattices Vector (R) : ", latvec);
    output::printM3(ofs ,"Supercell lattice vector : ", latvec_supercell);
    output::printM3(ofs, "Reciprocal lattice Vector (G): ", G);
    output::printM3(ofs, "GGT : ", GGT);

    ofs<<std::endl;
    return;
}

void UnitCell::print_cell_cif(const std::string &fn)const
{
	if (GlobalV::test_unitcell) ModuleBase::TITLE("UnitCell","print_cell_cif");
	
	if(GlobalV::MY_RANK!=0) return;//xiaohui add 2015-03-15

	std::stringstream ss;
	ss << GlobalV::global_out_dir << fn;
	
	std::ofstream ofs( ss.str().c_str() );
	ofs << "data_" << latName << std::endl;
	ofs << std::endl;
	//xiaohui modify 2015-03-25
	//ofs << "_audit_creation_method" << " generated by MESIA" << std::endl;
	ofs << "_audit_creation_method" << " generated by ABACUS" << std::endl;
	ofs << std::endl;
	ofs << "_cell_length_a " << a1.norm()*lat0*0.529177 << std::endl;
	ofs << "_cell_length_b " << a2.norm()*lat0*0.529177 << std::endl;
	ofs << "_cell_length_c " << a3.norm()*lat0*0.529177 << std::endl;

	//xiaohui modify and add 2014-12-21
	//ofs << "_cell_angle_alpha " << "90.00" << std::endl; 
	//ofs << "_cell_angle_beta " << "90.00" << std::endl; 
	//ofs << "_cell_angle_gamma " << "90.00" << std::endl; 
	//xiaohui modify alpha, beta and gamma 2015-09-29
	double angle_alpha = acos((a2 * a3) / (a2.norm() * a3.norm())) /ModuleBase::PI*180.0;
	double angle_beta = acos((a1 * a3) / (a1.norm() * a3.norm())) /ModuleBase::PI*180.0;
	double angle_gamma = acos((a1 * a2) / (a1.norm() * a2.norm())) /ModuleBase::PI*180.0;
	ofs << "_cell_angle_alpha " << angle_alpha << std::endl;
	ofs << "_cell_angle_beta " << angle_beta << std::endl;
	ofs << "_cell_angle_gamma " << angle_gamma << std::endl;
	ofs << std::endl;
	ofs << "_symmetry_space_group_name_H-M" << " " << std::endl;
	ofs << "_symmetry_Int_Tables_number" << " " << std::endl;
	ofs << std::endl;
	ofs << "loop_" << std::endl;
	ofs << "_atom_site_label" << std::endl;
	ofs << "_atom_site_fract_x" << std::endl;
	ofs << "_atom_site_fract_y" << std::endl;
	ofs << "_atom_site_fract_z" << std::endl;
    for (int it=0; it<ntype; it++)
    {
        for (int ia=0; ia<atoms[it].na; ia++)
        {
            ofs << atoms[it].label
            << " " << atoms[it].taud[ia].x
            << " " << atoms[it].taud[ia].y
            << " " << atoms[it].taud[ia].z << std::endl;
        }
    }
    ofs.close();
}

void UnitCell::print_cell_xyz(const std::string &fn)const
{
    if (GlobalV::test_unitcell) ModuleBase::TITLE("UnitCell","print_cell_xyz");

	if(GlobalV::MY_RANK!=0) return;//xiaohui add 2015-03-15

    std::stringstream ss;
    ss << GlobalV::global_out_dir << fn;

    std::ofstream ofs( ss.str().c_str() );

    ofs << nat << std::endl;
    ofs << latName << std::endl;
    for (int it=0; it<ntype; it++)
    {
        for (int ia=0; ia<atoms[it].na; ia++)
        {
            ofs << atoms[it].label
            << " " << atoms[it].tau[ia].x * lat0 * 0.529177
            << " " << atoms[it].tau[ia].y * lat0 * 0.529177
            << " " << atoms[it].tau[ia].z * lat0 * 0.529177 << std::endl;
        }
    }

    ofs.close();
    return;
}

void UnitCell::set_iat2itia(void)
{
	assert(nat>0);
	delete[] iat2it;
    delete[] iat2ia;
	this->iat2it = new int[nat];
    this->iat2ia = new int[nat];
	int iat=0;
	for(int it = 0;it < ntype;it++)
	{
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			this->iat2it[iat] = it;
            this->iat2ia[iat] = ia;
			++iat;
		}	
	}
	return;
}

void UnitCell::update_pos_tau(const double* pos)
{
    int iat = 0;
	for(int it = 0;it < this->ntype;it++)
	{
		Atom* atom = &this->atoms[it];
		for(int ia =0;ia< atom->na;ia++)
		{		
			if(atom->mbl[ia].x!=0)
			{
				atom->tau[ia].x = pos[3*iat] / this->lat0;
			}
			if(atom->mbl[ia].y!=0)
			{
				atom->tau[ia].y = pos[3*iat+1] / this->lat0;
			}
			if(atom->mbl[ia].z!=0)
			{
				atom->tau[ia].z = pos[3*iat+2] / this->lat0;
			}

			// the direct coordinates also need to be updated.
			atom->taud[ia] = atom->tau[ia] * this->GT;
			iat++;
		}
	}
	assert(iat == this->nat);
    return;
}

void UnitCell::update_pos_tau(const ModuleBase::Vector3<double>* posd_in)
{
    int iat = 0;
	for(int it = 0; it < this->ntype; ++it)
	{
		Atom* atom = &this->atoms[it];
		for(int ia = 0; ia < atom->na; ++ia)
		{		
			if(atom->mbl[ia].x!=0)
			{
				atom->tau[ia].x = posd_in[iat].x / this->lat0;
			}
			if(atom->mbl[ia].y!=0)
			{
				atom->tau[ia].y = posd_in[iat].y / this->lat0;
			}
			if(atom->mbl[ia].z!=0)
			{
				atom->tau[ia].z = posd_in[iat].z / this->lat0;
			}

			// the direct coordinates also need to be updated.
			atom->taud[ia] = atom->tau[ia] * this->GT;
			iat++;
		}
	}
	assert(iat == this->nat);
    return;
}

void UnitCell::update_pos_taud(const ModuleBase::Vector3<double>* posd_in)
{
    int iat = 0;
    for(int it = 0;it < this->ntype;it++)
    {
        Atom* atom = &this->atoms[it];
        for(int ia =0;ia< atom->na;ia++)
        {
            this->atoms[it].taud[ia] += posd_in[iat];
            iat++;
        }
    }
    assert(iat == this->nat);
    this->periodic_boundary_adjustment();
}

void UnitCell::update_vel(const ModuleBase::Vector3<double>* vel_in)
{
    int iat = 0;
    for(int it=0; it<this->ntype; ++it)
    {
        Atom* atom = &this->atoms[it];
        for(int ia=0; ia<atom->na; ++ia)
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
	for(int it=0; it<this->ntype; it++)
	{
		Atom* atom = &this->atoms[it];
		for(int ia=0; ia<atom->na; ia++)
		{
			// mohan update 2011-03-21
			if(atom->taud[ia].x<0) atom->taud[ia].x += 1.0;
			if(atom->taud[ia].y<0) atom->taud[ia].y += 1.0;
			if(atom->taud[ia].z<0) atom->taud[ia].z += 1.0;
			if(atom->taud[ia].x>=1.0) atom->taud[ia].x -= 1.0;
			if(atom->taud[ia].y>=1.0) atom->taud[ia].y -= 1.0;
			if(atom->taud[ia].z>=1.0) atom->taud[ia].z -= 1.0;

			if(atom->taud[ia].x<0 || atom->taud[ia].y<0
				|| atom->taud[ia].z<0 ||
				atom->taud[ia].x>=1.0 ||
				atom->taud[ia].y>=1.0 ||
				atom->taud[ia].z>=1.0)
			{
				GlobalV::ofs_warning << " it=" << it+1 << " ia=" << ia+1 << std::endl;
				GlobalV::ofs_warning << "d=" << atom->taud[ia].x << " " << 
				atom->taud[ia].y << " " << atom->taud[ia].z << std::endl;
				ModuleBase::WARNING_QUIT("Ions_Move_Basic::move_ions","the movement of atom is larger than the length of cell.");
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
    for (int i=0;i<ntype;i++)
    {
        atoms[i].bcast_atom(); // bcast tau array
    }
#endif
}

void UnitCell::save_cartesian_position(double* pos)const
{
    int iat=0;
	for(int it = 0;it < this->ntype;it++)
	{
		Atom* atom = &this->atoms[it];
		for(int ia =0; ia<atom->na; ia++)
		{	
			pos[3*iat  ] = atom->tau[ia].x*this->lat0;
			pos[3*iat+1] = atom->tau[ia].y*this->lat0;
			pos[3*iat+2] = atom->tau[ia].z*this->lat0;
            iat++;
        }
    }
    assert(iat == this->nat);
    return;
}

bool UnitCell::judge_big_cell(void)const
{
	double diameter = 2*GlobalV::SEARCH_RADIUS;

	double dis_x = this->omega/cross(a1*lat0, a2*lat0).norm();
	double dis_y = this->omega/cross(a2*lat0, a3*lat0).norm();
	double dis_z = this->omega/cross(a3*lat0, a1*lat0).norm();

	if(dis_x>diameter && dis_y>diameter && dis_z>diameter)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


#ifndef __CMD
void UnitCell::cal_ux()
{
	double amag, uxmod;
	int starting_it = 0;
    int starting_ia = 0;
	bool is_paraller;
	//do not sign feature in teh general case
	magnet.lsign_ = false;
	ModuleBase::GlobalFunc::ZEROS(magnet.ux_, 3);

	for(int it = 0;it<ntype;it++)
	{
        for(int ia=0;ia<atoms[it].na;ia++)
        {
            amag = pow(atoms[it].m_loc_[ia].x,2) + pow(atoms[it].m_loc_[ia].y,2) + pow(atoms[it].m_loc_[ia].z,2);
            if(amag > 1e-6)
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
        if (magnet.lsign_) break;
	}
	//whether the initial magnetizations is parallel
	for(int it = starting_it; it<ntype;it++)
	{
        for(int ia=0;ia<atoms[it].na;ia++)
        {
            if(it>starting_it || ia>starting_ia)
            {
                magnet.lsign_ = magnet.lsign_ && judge_parallel(magnet.ux_, atoms[it].m_loc_[ia]);
            }
        }
		
	}
	if(magnet.lsign_)
	{
		uxmod =  pow(magnet.ux_[0],2) + pow(magnet.ux_[1],2) +pow(magnet.ux_[2],2);
		if(uxmod<1e-6) 
		{
			ModuleBase::WARNING_QUIT("cal_ux","wrong uxmod");
		}
		for(int i = 0;i<3;i++)
		{
			magnet.ux_[i] *= 1/sqrt(uxmod);
		}
		//       std::cout<<"    Fixed quantization axis for GGA: "
		//<<std::setw(10)<<ux[0]<<"  "<<std::setw(10)<<ux[1]<<"  "<<std::setw(10)<<ux[2]<<std::endl;
	}
	return;
}
#endif

bool UnitCell::judge_parallel(double a[3], ModuleBase::Vector3<double> b)
{
   bool jp=false;
   double cross;
   cross = pow((a[1]*b.z-a[2]*b.y),2) +  pow((a[2]*b.x-a[0]*b.z),2) + pow((a[0]*b.y-a[1]*b.x),2);
   jp = (fabs(cross)<1e-6);
   return jp;
}
