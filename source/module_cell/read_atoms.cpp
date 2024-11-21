#include "unitcell.h"
#include "module_parameter/parameter.h"
#ifdef __LCAO
#include "../module_basis/module_ao/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif
#include "../module_base/timer.h"
#include "../module_base/constants.h"
#include "../module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/formatter.h"
#include <cstring>        // Peize Lin fix bug about strcmp 2016-08-02
#include <cassert>
#include <regex>
int UnitCell::read_atom_species(std::ifstream &ifa, std::ofstream &ofs_running)
{
    ModuleBase::TITLE("UnitCell","read_atom_species");

    int error = 0;//0 for correct, >0 for warning and quit

    delete[] atom_label;
    delete[] atom_mass;
    delete[] pseudo_fn;
    delete[] pseudo_type;
    delete[] orbital_fn;
    this->atom_mass  = new double[ntype]; //atom masses
    this->atom_label = new std::string[ntype]; //atom labels
    this->pseudo_fn  = new std::string[ntype]; //file name of pseudopotential
    this->pseudo_type = new std::string[ntype]; // type of pseudopotential
    this->orbital_fn = new std::string[ntype]; // filename of orbitals

    std::string word;
    //==========================================
    // read in information of each type of atom
    //==========================================
    if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "ATOMIC_SPECIES") )
    {    
        ifa.ignore(300, '\n');
        ModuleBase::GlobalFunc::OUT(ofs_running,"ntype",ntype);
        for (int i = 0;i < ntype;i++)
        {
            std::string one_line, one_string;
            std::getline(ifa, one_line);
            std::stringstream ss;
            ss << one_line;
            ss >> atom_label[i] >> atom_mass[i];
            pseudo_fn[i] = "auto";
            pseudo_type[i] = "auto";

            if(!PARAM.inp.use_paw)
            {
                bool end = false;
                if (ss >> one_string)
                {
                    if (one_string[0] != '#')
                    {
                        pseudo_fn[i] = one_string;
                    }
                    else
                    {
                        end = true;
                    }
                }

                if (!end && ss >> one_string && one_string[0] != '#')
                {
                    if (one_string == "auto" || one_string == "upf" || one_string == "vwr" || one_string == "upf201" || one_string == "blps")
                    {
                        pseudo_type[i] = one_string;
                    }
                    else if (one_string == "1/r")
                    {
                        atoms[i].coulomb_potential = true;
                    }
                    else
                    {
                        GlobalV::ofs_warning << "unrecongnized pseudopotential type: " << one_string << ", check your STRU file." << std::endl;
                        ModuleBase::WARNING_QUIT("read_atom_species", "unrecongnized pseudo type.");
                    }
                }

                if(PARAM.inp.test_pseudo_cell==2) 
                {
                    ofs_running << "\n" << std::setw(6) << atom_label[i] 
                            << std::setw(12) << atom_mass[i] 
                            << std::setw(18) << pseudo_fn[i]
                            << std::setw(18) << pseudo_type[i];
                }

                // Peize Lin test for bsse 2021.04.07
                const std::string bsse_label = "empty";
                this->atoms[i].flag_empty_element = 
                    (search( atom_label[i].begin(), atom_label[i].end(), bsse_label.begin(), bsse_label.end() ) != atom_label[i].end())
                    ? true : false;
            }
        }
    }

    if(
        (PARAM.inp.basis_type == "lcao")
      ||(PARAM.inp.basis_type == "lcao_in_pw")
      ||(
          (PARAM.inp.basis_type == "pw")
        &&(PARAM.inp.psi_initializer)
        &&(PARAM.inp.init_wfc.substr(0, 3) == "nao")
        )
        || PARAM.inp.onsite_radius > 0.0
    )
    {
        if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "NUMERICAL_ORBITAL") )
        {
            for(int i=0; i<ntype; i++)
            {
                ifa >> orbital_fn[i];
            }
        }    
        // caoyu add 2021-03-16
        if(PARAM.globalv.deepks_setorb)
        {
            if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "NUMERICAL_DESCRIPTOR")) {
                ifa >> descriptor_file;
            }
        }
        else{
            descriptor_file = PARAM.inp.orbital_dir + orbital_fn[0];
        }
    }
#ifdef __LCAO
    // Peize Lin add 2016-09-23
#ifdef __MPI 
#ifdef __EXX
    if( GlobalC::exx_info.info_global.cal_exx || PARAM.inp.rpa )
    {
        if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "ABFS_ORBITAL") )
        {
            for(int i=0; i<ntype; i++)
            {
                std::string ofile;
                ifa >> ofile;
                GlobalC::exx_info.info_ri.files_abfs.push_back(ofile);
            }
        }
    }

#endif // __EXX
#endif // __MPI
#endif // __LCAO
    //==========================
    // read in lattice constant
    //==========================
    if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_CONSTANT") )
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifa, lat0);
        if(lat0<=0.0)
        {
            ModuleBase::WARNING_QUIT("read_atom_species","lat0<=0.0");
        }
        lat0_angstrom = lat0 * 0.529177 ;
        ModuleBase::GlobalFunc::OUT(ofs_running,"lattice constant (Bohr)",lat0);
        ModuleBase::GlobalFunc::OUT(ofs_running,"lattice constant (Angstrom)",lat0_angstrom);
        this->tpiba  = ModuleBase::TWO_PI / lat0;
        this->tpiba2 = tpiba * tpiba;
    }

    //===========================
    // Read in latticies vector
    //===========================
    if(latName=="none"){
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa,
                                               "LATTICE_PARAMETERS",
                                               true,
                                               false)) {
            ModuleBase::WARNING_QUIT("UnitCell::read_atom_species","do not use LATTICE_PARAMETERS without explicit specification of lattice type");
        }
        if( !ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_VECTORS") )
        {
            ModuleBase::WARNING_QUIT("UnitCell::read_atom_species","Please set LATTICE_VECTORS in STRU file");
        }
        else if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_VECTORS") )
        {
            // Reading lattice vectors. notice
            // here that only one cpu read these
            // parameters.
            ifa >> latvec.e11 >> latvec.e12;
            ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e13);
            ifa >> latvec.e21 >> latvec.e22;
            ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e23);
            ifa >> latvec.e31 >> latvec.e32;
            ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
        }
    }//supply lattice vectors
    else{
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa,
                                               "LATTICE_VECTORS",
                                               true,
                                               false)) {
            ModuleBase::WARNING_QUIT("UnitCell::read_atom_species","do not use LATTICE_VECTORS along with explicit specification of lattice type");
        }
        if(latName=="sc"){//simple-cubic, ibrav = 1
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = 0.0; latvec.e22 = 1.0;    latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0;    latvec.e33 = 1.0;
        }
        else if(latName=="fcc"){//face-centered cubic, ibrav = 2
            latvec.e11 =-0.5; latvec.e12 = 0.0; latvec.e13 = 0.5;
            latvec.e21 = 0.0; latvec.e22 = 0.5;    latvec.e23 = 0.5;
            latvec.e31 =-0.5; latvec.e32 = 0.5;    latvec.e33 = 0.0;
        }
        else if(latName=="bcc"){//body-centered cubic, ibrav = 3
            latvec.e11 = 0.5; latvec.e12 = 0.5; latvec.e13 = 0.5;
            latvec.e21 =-0.5; latvec.e22 = 0.5;    latvec.e23 = 0.5;
            latvec.e31 =-0.5; latvec.e32 =-0.5;    latvec.e33 = 0.5;
        }
        else if(latName=="hexagonal"){//hexagonal, ibrav = 4
            double e22 = sqrt(3.0) / 2.0;
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 =-0.5; latvec.e22 = e22; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0;    latvec.e33 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
            }
        }
        else if(latName=="trigonal"){//trigonal, ibrav = 5
            double t1 = 0.0;
            double t2 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                double cosab=0.0;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cosab);
                t1 = sqrt(1.0 + 2.0*cosab);
                t2 = sqrt(1.0 - cosab);
            }
            double e11 = t2 / sqrt(2.0);
            double e12 = -t2 / sqrt(6.0);
            double e13 = t1 / sqrt(3.0);
            double e22 = sqrt(2.0) * t2 / sqrt(3.0);
        
            latvec.e11 = e11; latvec.e12 = e12; latvec.e13 = e13;
            latvec.e21 = 0.0; latvec.e22 = e22;    latvec.e23 = e13;
            latvec.e31 =-e11; latvec.e32 = e12;    latvec.e33 = e13;
        }
        else if(latName=="st"){//simple tetragonal, ibrav= 6
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0;    latvec.e33 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
            }
        }
        else if(latName=="bct"){//body-centered tetragonal, ibrav = 7
            double cba = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cba);
                cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 =-0.5; latvec.e13 = cba;
            latvec.e21 = 0.5; latvec.e22 = 0.5; latvec.e23 = cba;
            latvec.e31 =-0.5; latvec.e32 =-0.5;    latvec.e33 = cba;
        }
        else if(latName=="so"){//simple orthorhombic, ibrav = 8
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = 0.0; latvec.e22 = 0.0;    latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0;    latvec.e33 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> latvec.e22;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
            }
        }
        else if(latName=="baco"){//base-centered orthorhombic, ibrav = 9
            latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 =-0.5; latvec.e22 = 0.0;    latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0;    latvec.e33 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> latvec.e12;
                latvec.e12 = latvec.e12 / 2.0;
                latvec.e22 = latvec.e12;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
            }
        }
        else if(latName=="fco"){//face-centered orthorhombic, ibrav = 10
            double bba = 0.0; double cba = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cba);
                bba = bba / 2.0; cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 = cba;
            latvec.e21 = 0.5; latvec.e22 = bba;    latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = bba;    latvec.e33 = cba;
        }
        else if(latName=="bco"){//body-centered orthorhombic, ibrav = 11
            double bba = 0.0; double cba = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cba);
                bba = bba / 2.0; cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 = bba; latvec.e13 = cba;
            latvec.e21 =-0.5; latvec.e22 = bba;    latvec.e23 = cba;
            latvec.e31 =-0.5; latvec.e32 =-bba;    latvec.e33 = cba;
        }
        else if(latName=="sm"){//simple monoclinic, ibrav = 12
            double bba = 0.0; double cba = 0.0;
            double cosab = 0.0;
            double e21 = 0.0; double e22 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cosab);
                e21 = bba * cosab;
                e22 = bba * sqrt(1.0-cosab*cosab);
            }
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = e21; latvec.e22 = e22;    latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0;    latvec.e33 = cba;
        }
        else if(latName=="bacm"){//base-centered monoclinic, ibrav = 13
            double bba = 0.0; double cba = 0.0;
            double cosab = 0.0;
            double e21 = 0.0; double e22 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cosab);
                e21 = bba * cosab;
                e22 = bba * sqrt(1.0-cosab*cosab);
                cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 =-cba;
            latvec.e21 = e21; latvec.e22 = e22;    latvec.e23 = 0.0;
            latvec.e31 = 0.5; latvec.e32 = 0.0;    latvec.e33 = cba;
        }
        else if(latName=="triclinic"){//triclinic, ibrav = 14
            double bba = 0.0; double cba = 0.0;
            double cosab = 0.0; double cosac = 0.0;
            double cosbc = 0.0; double sinab = 0.0;
            double term = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba >> cosab >> cosac;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cosbc);
                sinab = sqrt(1.0-cosab*cosab);
            }
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = bba * cosab;
            latvec.e22 = bba * sinab;
            latvec.e23 = 0.0;
            latvec.e31 = cba * cosac;
            latvec.e32 = cba * (cosbc - cosac*cosab) / sinab;
            term = 1.0 + 2.0 * cosab*cosac*cosbc - cosab*cosab - cosac*cosac - cosbc*cosbc;
            term = sqrt(term)/sinab;
            latvec.e33 = cba * term;
        }
        else{ 
            std::cout << "latname is : " << latName << std::endl;
            ModuleBase::WARNING_QUIT("UnitCell::read_atom_species","latname not supported!");
        }
    }

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
    return 0;
}

#include "../module_base/mathzone.h"
// Read atomic positions
// return 1: no problem.
// return 0: some problems.
bool UnitCell::read_atom_positions(std::ifstream &ifpos, std::ofstream &ofs_running, std::ofstream &ofs_warning)
{
    ModuleBase::TITLE("UnitCell","read_atom_positions");

    if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifpos, "ATOMIC_POSITIONS"))
    {
        ModuleBase::GlobalFunc::READ_VALUE( ifpos, Coordinate);
        if(Coordinate != "Cartesian" 
            && Coordinate != "Direct" 
            && Coordinate != "Cartesian_angstrom"
            && Coordinate != "Cartesian_au"
            && Coordinate != "Cartesian_angstrom_center_xy"
            && Coordinate != "Cartesian_angstrom_center_xz"
            && Coordinate != "Cartesian_angstrom_center_yz"
            && Coordinate != "Cartesian_angstrom_center_xyz"
            )
        {
            ModuleBase::WARNING("read_atom_position","Cartesian or Direct?");
            ofs_warning << " There are several options for you:" << std::endl;
            ofs_warning << " Direct" << std::endl;
            ofs_warning << " Cartesian_angstrom" << std::endl;
            ofs_warning << " Cartesian_au" << std::endl;
            ofs_warning << " Cartesian_angstrom_center_xy" << std::endl;
            ofs_warning << " Cartesian_angstrom_center_xz" << std::endl;
            ofs_warning << " Cartesian_angstrom_center_yz" << std::endl;
            ofs_warning << " Cartesian_angstrom_center_xyz" << std::endl;
            return false; // means something wrong
        }

        ModuleBase::Vector3<double> v;
        ModuleBase::Vector3<int> mv;
        int na = 0;
        this->nat = 0;

        //======================================
        // calculate total number of atoms
        // and adjust the order of atom species
        //======================================
        assert(ntype>0);
        for (int it = 0;it < ntype; it++)
        {
            ofs_running << "\n READING ATOM TYPE " << it+1 << std::endl;
            
            //=======================================
            // (1) read in atom label
            // start magnetization
            //=======================================
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, atoms[it].label);
            if(this->atoms[it].label != this->atom_label[it])
            {
                ofs_warning << " Label orders in ATOMIC_POSITIONS and ATOMIC_SPECIES sections do not match!" << std::endl;
                ofs_warning << " Label read from ATOMIC_POSITIONS is " << this->atoms[it].label << std::endl; 
                ofs_warning << " Label from ATOMIC_SPECIES is " << this->atom_label[it] << std::endl;
                return false;
            }
            ModuleBase::GlobalFunc::OUT(ofs_running, "atom label",atoms[it].label);

            bool set_element_mag_zero = false;
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, magnet.start_magnetization[it] );

#ifndef __SYMMETRY
            //===========================================
            // (2) read in numerical orbital information
            // int atoms[it].nwl
            // int* atoms[it].l_nchi;
            //===========================================

            if ((PARAM.inp.basis_type == "lcao")||(PARAM.inp.basis_type == "lcao_in_pw"))
            {
                std::string orbital_file = PARAM.inp.orbital_dir + orbital_fn[it];
                this->read_orb_file(it, orbital_file, ofs_running, &(atoms[it]));
            }
            else if(PARAM.inp.basis_type == "pw")
            {
                if ((PARAM.inp.psi_initializer)&&(PARAM.inp.init_wfc.substr(0, 3) == "nao"))
                {
                    std::string orbital_file = PARAM.inp.orbital_dir + orbital_fn[it];
                    this->read_orb_file(it, orbital_file, ofs_running, &(atoms[it]));
                }
                else
                {
                    this->atoms[it].nw = 0;
                    this->atoms[it].nwl = 2;
                    //std::cout << lmaxmax << std::endl;
                    if ( lmaxmax != 2 )
                    {
                        this->atoms[it].nwl = lmaxmax;
                    }
                    this->atoms[it].l_nchi.resize(this->atoms[it].nwl+1, 0);
                    for(int L=0; L<atoms[it].nwl+1; L++)
                    {
                        this->atoms[it].l_nchi[L] = 1;
                        // calculate the number of local basis(3D)
                        this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];
                        std::stringstream ss;
                        ss << "L=" << L << ", number of zeta";
                        ModuleBase::GlobalFunc::OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
                    }
                }
            } // end basis type
#endif

            //=========================
            // (3) read in atom number
            //=========================
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, na);
            this->atoms[it].na = na;

            ModuleBase::GlobalFunc::OUT(ofs_running,"number of atom for this type",na);

            this->nat += na;

            /**
             * liuyu update 2023-05-11
             * In order to employ the DP model as esolver,
             * all atom types must be specified in the `STRU` in the order consistent with that of the DP model,
             * even if the number of atoms is zero!
             */
            if (na < 0)
            {
                ModuleBase::WARNING("read_atom_positions", " atom number < 0.");
                return false;
            }
            if (na > 0)
            {
                atoms[it].tau.resize(na, ModuleBase::Vector3<double>(0,0,0));
                atoms[it].dis.resize(na, ModuleBase::Vector3<double>(0,0,0));
                atoms[it].taud.resize(na, ModuleBase::Vector3<double>(0,0,0));
                atoms[it].vel.resize(na, ModuleBase::Vector3<double>(0,0,0));
                atoms[it].mbl.resize(na, ModuleBase::Vector3<int>(0,0,0));
                atoms[it].mag.resize(na, 0);
                atoms[it].angle1.resize(na, 0);
                atoms[it].angle2.resize(na, 0);
                atoms[it].m_loc_.resize(na, ModuleBase::Vector3<double>(0,0,0));
                atoms[it].lambda.resize(na, ModuleBase::Vector3<double>(0,0,0));
                atoms[it].constrain.resize(na, ModuleBase::Vector3<int>(0,0,0));
                atoms[it].mass = this->atom_mass[it]; //mohan add 2011-11-07 
                for (int ia = 0;ia < na; ia++)
                {
                 // modify the reading of frozen ions and velocities  -- Yuanbo Li 2021/8/20
                    ifpos >> v.x >> v.y >> v.z;
                    mv.x = true ;
                    mv.y = true ;
                    mv.z = true ;
                    atoms[it].vel[ia].set(0,0,0);
                    atoms[it].mag[ia]=magnet.start_magnetization[it];//if this line is used, default startmag_type would be 2
                    atoms[it].angle1[ia]=0;
                    atoms[it].angle2[ia]=0;
                    atoms[it].m_loc_[ia].set(0,0,0);
                    atoms[it].lambda[ia].set(0,0,0);
                    atoms[it].constrain[ia].set(0,0,0);

                    std::string tmpid;
                    tmpid = ifpos.get();

                    if( (int)tmpid[0] < 0 )
                    {
                        std::cout << "read_atom_positions, mismatch in atom number for atom type: " << atoms[it].label << std::endl;
                        exit(1); 
                    }

                    bool input_vec_mag=false;
                    bool input_angle_mag=false;
                    // read if catch goodbit before "\n" and "#"
                    while ( (tmpid != "\n") && (ifpos.good()) && (tmpid !="#") )
                    {
                        tmpid = ifpos.get() ;
                        // old method of reading frozen ions
                        char tmp = (char)tmpid[0];
                        if ( tmp >= 48 && tmp <= 57 )
                        {
                                mv.x = std::stoi(tmpid);
                                ifpos >> mv.y >> mv.z ;
                        }
                        // new method of reading frozen ions and velocities
                        if ( tmp >= 'a' && tmp <='z')
                        {
                            ifpos.putback(tmp);
                            ifpos >> tmpid;
                        }
                        if ( tmpid == "m" )
                        {
                                ifpos >> mv.x >> mv.y >> mv.z ;
                        }
                        else if ( tmpid == "v" ||tmpid == "vel" || tmpid == "velocity" )
                        {
                                ifpos >> atoms[it].vel[ia].x >> atoms[it].vel[ia].y >> atoms[it].vel[ia].z;
                        }
                        else if ( tmpid == "mag" || tmpid == "magmom")
                        {
                            set_element_mag_zero = true;
                            double tmpamg=0;
                            ifpos >> tmpamg;
                            tmp=ifpos.get();
                            while (tmp==' ')
                            {
                                tmp=ifpos.get();
                            }
                            
                            if((tmp >= 48 && tmp <= 57) or tmp=='-')
                            {
                                ifpos.putback(tmp);
                                ifpos >> atoms[it].m_loc_[ia].y>>atoms[it].m_loc_[ia].z;
                                atoms[it].m_loc_[ia].x=tmpamg;
                                atoms[it].mag[ia]=sqrt(pow(atoms[it].m_loc_[ia].x,2)+pow(atoms[it].m_loc_[ia].y,2)+pow(atoms[it].m_loc_[ia].z,2));
                                input_vec_mag=true;
                                
                            }
                            else
                            {
                                ifpos.putback(tmp);
                                atoms[it].mag[ia]=tmpamg;
                            }
                            
                            // atoms[it].mag[ia];
                        }
                        else if ( tmpid == "angle1")
                        {
                                ifpos >> atoms[it].angle1[ia];
                                atoms[it].angle1[ia]=atoms[it].angle1[ia]/180 *ModuleBase::PI;
                                input_angle_mag=true;
                                set_element_mag_zero = true;
                        }
                        else if ( tmpid == "angle2")
                        {
                                ifpos >> atoms[it].angle2[ia];
                                atoms[it].angle2[ia]=atoms[it].angle2[ia]/180 *ModuleBase::PI;
                                input_angle_mag=true;
                                set_element_mag_zero = true;
                        }   
                        else if ( tmpid == "lambda")
                        {
                            double tmplam=0;
                            ifpos >> tmplam;
                            tmp=ifpos.get();
                            while (tmp==' ')
                            {
                                tmp=ifpos.get();
                            }
                            if((tmp >= 48 && tmp <= 57) or tmp=='-')
                            {
                                ifpos.putback(tmp);
                                ifpos >> atoms[it].lambda[ia].y>>atoms[it].lambda[ia].z;
                                atoms[it].lambda[ia].x=tmplam;
                            }
                            else
                            {
                                ifpos.putback(tmp);
                                atoms[it].lambda[ia].z=tmplam;
                            }
                            atoms[it].lambda[ia].x /= ModuleBase::Ry_to_eV;
                            atoms[it].lambda[ia].y /= ModuleBase::Ry_to_eV;
                            atoms[it].lambda[ia].z /= ModuleBase::Ry_to_eV;
                        }
                        else if ( tmpid == "sc")
                        {
                            double tmplam=0;
                            ifpos >> tmplam;
                            tmp=ifpos.get();
                            while (tmp==' ')
                            {
                                tmp=ifpos.get();
                            }
                            if((tmp >= 48 && tmp <= 57) or tmp=='-')
                            {
                                ifpos.putback(tmp);
                                ifpos >> atoms[it].constrain[ia].y>>atoms[it].constrain[ia].z;
                                atoms[it].constrain[ia].x=tmplam;
                            }
                            else
                            {
                                ifpos.putback(tmp);
                                atoms[it].constrain[ia].z=tmplam;
                            }
                        } 
                    }
                    // move to next line
                    while ( (tmpid != "\n") && (ifpos.good()) )
                    {
                            tmpid = ifpos.get();
                    }
                    std::string mags;
                    //cout<<"mag"<<atoms[it].mag[ia]<<"angle1"<<atoms[it].angle1[ia]<<"angle2"<<atoms[it].angle2[ia]<<'\n';

                    if(PARAM.inp.nspin==4)
                    {
                        if(PARAM.inp.noncolin)
                        {
                            if(input_angle_mag)
                            {
                                atoms[it].m_loc_[ia].z = atoms[it].mag[ia] *
                                    cos(atoms[it].angle1[ia]);
                                if(std::abs(sin(atoms[it].angle1[ia])) > 1e-10 )
                                {
                                    atoms[it].m_loc_[ia].x = atoms[it].mag[ia] *
                                        sin(atoms[it].angle1[ia]) * cos(atoms[it].angle2[ia]);
                                    atoms[it].m_loc_[ia].y = atoms[it].mag[ia] *
                                        sin(atoms[it].angle1[ia]) * sin(atoms[it].angle2[ia]);
                                }
                            }
                            else if (input_vec_mag)
                            {
                                double mxy=sqrt(pow(atoms[it].m_loc_[ia].x,2)+pow(atoms[it].m_loc_[ia].y,2));
                                atoms[it].angle1[ia]=atan2(mxy,atoms[it].m_loc_[ia].z);
                                if(mxy>1e-8)
                                {
                                    atoms[it].angle2[ia]=atan2(atoms[it].m_loc_[ia].y,atoms[it].m_loc_[ia].x);
                                }
                            }
                            else
                            {
                                atoms[it].m_loc_[ia].x = 0;
                                atoms[it].m_loc_[ia].y = 0;
                                atoms[it].m_loc_[ia].z = atoms[it].mag[ia];
                            }
                        }
                        else
                        {
                            atoms[it].m_loc_[ia].x = 0;
                            atoms[it].m_loc_[ia].y = 0;
                            atoms[it].m_loc_[ia].z = atoms[it].mag[ia];
                        }

                        //print only ia==0 && mag>0 to avoid too much output
                        //print when ia!=0 && mag[ia] != mag[0] to avoid too much output
                        if(ia==0 || (ia!=0 
                            && (atoms[it].m_loc_[ia].x != atoms[it].m_loc_[0].x 
                            || atoms[it].m_loc_[ia].y != atoms[it].m_loc_[0].y 
                            || atoms[it].m_loc_[ia].z != atoms[it].m_loc_[0].z)))
                        {
                            //use a stringstream to generate string: "concollinear magnetization of element it is:"
                            std::stringstream ss;
                            ss << "magnetization of element " << it+1;
                            if(ia!=0) 
                            {
                                ss<<" (atom"<<ia+1<<")";
                            }
                            ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(),atoms[it].m_loc_[ia].x, atoms[it].m_loc_[ia].y, atoms[it].m_loc_[ia].z);
                        }
                        ModuleBase::GlobalFunc::ZEROS(magnet.ux_ ,3);
                    }
                    else if(PARAM.inp.nspin==2)
                    {
                        atoms[it].m_loc_[ia].x = atoms[it].mag[ia];
                        //print only ia==0 && mag>0 to avoid too much output
                        //print when ia!=0 && mag[ia] != mag[0] to avoid too much output
                        if(ia==0 || (ia!=0 && atoms[it].mag[ia] != atoms[it].mag[0]))
                        {
                            //use a stringstream to generate string: "cocollinear magnetization of element it is:"
                            std::stringstream ss;
                            ss << "magnetization of element " << it+1;
                            if(ia!=0) 
                            {
                                ss<<" (atom"<<ia+1<<")";
                            }
                            ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(),atoms[it].mag[ia]);
                        }
                    }
            
                    if(Coordinate=="Direct")
                    {
                        // change v from direct to cartesian,
                        // the unit is GlobalC::sf.lat0
                        atoms[it].taud[ia] = v;
                        atoms[it].tau[ia] = v * latvec;
                    }
                    else if(Coordinate=="Cartesian")
                    {
                        atoms[it].tau[ia] = v ;// in unit lat0
                        //std::cout << " T=" << it << " I=" << ia << " tau=" << atoms[it].tau[ia].x << " " << 
                        //atoms[it].tau[ia].y << " " << atoms[it].tau[ia].z << std::endl;
                    }
                    else if(Coordinate=="Cartesian_angstrom")
                    {
                        atoms[it].tau[ia] = v / 0.529177 / lat0;
                    }    
                    else if(Coordinate=="Cartesian_angstrom_center_xy")
                    {
                        // calculate lattice center 
                        latcenter.x = (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
                        latcenter.y = (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
                        latcenter.z = 0.0;
                        atoms[it].tau[ia] = v / 0.529177 / lat0 + latcenter; 
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_xz")
                    {
                        // calculate lattice center 
                        latcenter.x = (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
                        latcenter.y = 0.0; 
                        latcenter.z = (latvec.e13 + latvec.e23 + latvec.e33)/2.0;    
                        atoms[it].tau[ia] = v / 0.529177 / lat0 + latcenter; 
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_yz")
                    {
                        // calculate lattice center 
                        latcenter.x = 0.0; 
                        latcenter.y = (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
                        latcenter.z = (latvec.e13 + latvec.e23 + latvec.e33)/2.0;    
                        atoms[it].tau[ia] = v / 0.529177 / lat0 + latcenter; 
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_xyz")
                    {
                        // calculate lattice center 
                        latcenter.x = (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
                        latcenter.y = (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
                        latcenter.z = (latvec.e13 + latvec.e23 + latvec.e33)/2.0;    
                        atoms[it].tau[ia] = v / 0.529177 / lat0 + latcenter; 
                    }
                    else if(Coordinate=="Cartesian_au")
                    {
                        atoms[it].tau[ia] = v / lat0;
                    }

                    if(Coordinate=="Cartesian" || 
                        Coordinate=="Cartesian_angstrom" || 
                        Coordinate=="Cartesian_angstrom_center_xy" || 
                        Coordinate=="Cartesian_angstrom_center_xz" || 
                        Coordinate=="Cartesian_angstrom_center_yz" || 
                        Coordinate=="Cartesian_angstrom_center_xyz" || 
                        Coordinate=="Cartesian_au")
                    {
                        double dx=0.0;
                        double dy=0.0;
                        double dz=0.0;
                        ModuleBase::Mathzone::Cartesian_to_Direct(atoms[it].tau[ia].x, atoms[it].tau[ia].y, atoms[it].tau[ia].z,
                        latvec.e11, latvec.e12, latvec.e13,
                        latvec.e21, latvec.e22, latvec.e23,
                        latvec.e31, latvec.e32, latvec.e33,
                        dx,dy,dz);
                    
                        atoms[it].taud[ia].x = dx;
                        atoms[it].taud[ia].y = dy;
                        atoms[it].taud[ia].z = dz;

                    }
                    
                    if(!PARAM.inp.fixed_atoms)
                    {
                        atoms[it].mbl[ia] = mv;
                    }
                    else
                    {
                        atoms[it].mbl[ia] = 0.0;
                        atoms[it].mbl[ia].print();
                    }
                    atoms[it].dis[ia].set(0, 0, 0);
                }//endj
            }    // end na
            // reset some useless parameters
            if (set_element_mag_zero)
            {
                magnet.start_magnetization[it] = 0.0;
            }
        } // end for ntype
        // Start Autoset magnetization
        // defaultly set a finite magnetization if magnetization is not specified
        int autoset_mag = 1;
        for (int it = 0;it < ntype; it++)
        {
            for (int ia = 0;ia < this->atoms[it].na; ia++)
            {
                if(std::abs(atoms[it].mag[ia]) > 1e-5)
                {
                    autoset_mag = 0;
                    break;
                }
            }
        }
        if (autoset_mag)
        {
            if(PARAM.inp.nspin==4)
            {
                for (int it = 0;it < ntype; it++)
                {
                    for (int ia = 0;ia < this->atoms[it].na; ia++)
                    {
                        atoms[it].m_loc_[ia].x = 1.0;
                        atoms[it].m_loc_[ia].y = 1.0;
                        atoms[it].m_loc_[ia].z = 1.0;
                        atoms[it].mag[ia] = sqrt(pow(atoms[it].m_loc_[ia].x,2)+pow(atoms[it].m_loc_[ia].y,2)+pow(atoms[it].m_loc_[ia].z,2));
                        ModuleBase::GlobalFunc::OUT(ofs_running,"Autoset magnetism for this atom", 1.0, 1.0, 1.0);
                    }
                }
            }
            else if(PARAM.inp.nspin==2)
            {
                for (int it = 0;it < ntype; it++)
                {
                    for (int ia = 0;ia < this->atoms[it].na; ia++)
                    {
                        atoms[it].mag[ia] = 1.0;
                        atoms[it].m_loc_[ia].x = atoms[it].mag[ia];
                        ModuleBase::GlobalFunc::OUT(ofs_running,"Autoset magnetism for this atom", 1.0);
                    }
                }
            }
        }
        // End Autoset magnetization
    }   // end scan_begin

//check if any atom can move in MD
    if(!this->if_atoms_can_move() && PARAM.inp.calculation=="md" && PARAM.inp.esolver_type!="tddft")
    {
        ModuleBase::WARNING("read_atoms", "no atom can move in MD!");
        return false;
    } 

    ofs_running << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running,"TOTAL ATOM NUMBER",nat);

    // mohan add 2010-06-30    
    //xiaohui modify 2015-03-15, cancel outputfile "STRU_READIN.xyz"
    //this->print_cell_xyz("STRU_READIN.xyz");
    this->check_dtau();

    if ( this->check_tau() )
    {

    }
    else
    {
        return false;
    }
    this->print_tau();
    //xiaohui modify 2015-03-15, cancel outputfile "STRU_READIN.xyz"
    //this->print_cell_xyz("STRU_READIN_ADJUST.xyz");

    return true;
}//end read_atom_positions

bool UnitCell::check_tau() const {
    ModuleBase::TITLE("UnitCell","check_tau");
    ModuleBase::timer::tick("UnitCell","check_tau");
    
    ModuleBase::Vector3<double> diff = 0.0;
    double norm = 0.0;
    double tolerence_bohr = 1.0e-3;

    //GlobalV::ofs_running << "\n Output nearest atom not considering periodic boundary condition" << std::endl;
    //GlobalV::ofs_running << " " << std::setw(5) << "TYPE" << std::setw(6) << "INDEX" 
    //<< std::setw(20) << "NEAREST(Bohr)" 
    //<< std::setw(20) << "NEAREST(Angstrom)" << std::endl; 
    for(int T1=0; T1< this->ntype; T1++)
    {
        for(int I1=0; I1< this->atoms[T1].na; I1++)
        {    
            double shortest_norm = 10000.0; // a large number
            //int nearest_atom_type = 0;
            //int nearest_atom_index = 0;
            for(int T2=0; T2<this->ntype; T2++)
            {
                for(int I2=0; I2<this->atoms[T2].na; I2++)
                {
                    if(T1==T2 && I1==I2)
                    {
                        shortest_norm = 0.0;
                        //nearest_atom_type = T1;
                        //nearest_atom_index = I2;
                        // self atom
                    }
                    else
                    {
                        diff = atoms[T1].tau[I1] - atoms[T2].tau[I2];
                        norm = diff.norm() * lat0;
                        if( shortest_norm > norm )
                        {
                            shortest_norm = norm;
                            //nearest_atom_type = T2;
                            //nearest_atom_index = I2;
                        }
                        if( norm < tolerence_bohr ) // unit is Bohr
                        {    
                            GlobalV::ofs_warning << " two atoms are too close!" << std::endl;
                            GlobalV::ofs_warning << " type:" << this->atoms[T1].label << " atom " << I1 + 1 << std::endl; 
                            GlobalV::ofs_warning << " type:" << this->atoms[T2].label << " atom " << I2 + 1 << std::endl; 
                            GlobalV::ofs_warning << " distance = " << norm << " Bohr" << std::endl;
                            return false;
                        }
                    }
                }
            }
            //GlobalV::ofs_running << " " << std::setw(5) << atoms[T1].label << std::setw(6) << I1+1 
            //<< std::setw(20) << shortest_norm  
            //<< std::setw(20) << shortest_norm * ModuleBase::BOHR_TO_A << std::endl;
        }
    }

    ModuleBase::timer::tick("UnitCell","check_tau");
    return true;
}

void UnitCell::print_stru_file(const std::string& fn, 
                               const int& nspin,
                               const bool& direct,
                               const bool& vel,
                               const bool& magmom,
                               const bool& orb,
                               const bool& dpks_desc,
                               const int& iproc) const
{
    ModuleBase::TITLE("UnitCell","print_stru_file");
    if (iproc != 0) {
        return; // old: if(GlobalV::MY_RANK != 0) return;
    }
    // ATOMIC_SPECIES
    std::string str = "ATOMIC_SPECIES\n";
    for(int it=0; it<ntype; it++){ str += FmtCore::format("%s %8.4f %s %s\n", atom_label[it], atom_mass[it], pseudo_fn[it], pseudo_type[it]); }
    // NUMERICAL_ORBITAL
    if(orb)
    {
        str += "\nNUMERICAL_ORBITAL\n";
        for(int it = 0; it < ntype; it++) { str += orbital_fn[it] + "\n"; }
    }
    // NUMERICAL_DESCRIPTOR
    if(dpks_desc) { str += "\nNUMERICAL_DESCRIPTOR\n" + descriptor_file + "\n"; }
    // LATTICE_CONSTANT
    str += "\nLATTICE_CONSTANT\n" + FmtCore::format("%-.10f\n", lat0);
    // LATTICE_VECTORS
    str += "\nLATTICE_VECTORS\n";
    str += FmtCore::format("%20.10f%20.10f%20.10f\n", latvec.e11, latvec.e12, latvec.e13);
    str += FmtCore::format("%20.10f%20.10f%20.10f\n", latvec.e21, latvec.e22, latvec.e23);
    str += FmtCore::format("%20.10f%20.10f%20.10f\n", latvec.e31, latvec.e32, latvec.e33);
    // ATOMIC_POSITIONS
    str += "\nATOMIC_POSITIONS\n";
    const std::string scale = direct? "Direct": "Cartesian";
    int nat_ = 0; // counter iat, for printing out Mulliken magmom who is indexed by iat
    str += scale + "\n";
    for(int it = 0; it < ntype; it++)
    {
        str += "\n" + atoms[it].label + " #label\n";
        str += FmtCore::format("%-8.4f #magnetism\n", magnet.start_magnetization[it]);
        str += FmtCore::format("%d #number of atoms\n", atoms[it].na);
        for(int ia = 0; ia < atoms[it].na; ia++)
        {
            // output position
            const double& x = direct? atoms[it].taud[ia].x: atoms[it].tau[ia].x;
            const double& y = direct? atoms[it].taud[ia].y: atoms[it].tau[ia].y;
            const double& z = direct? atoms[it].taud[ia].z: atoms[it].tau[ia].z;
            str += FmtCore::format("%20.10f%20.10f%20.10f", x, y, z);
            str += FmtCore::format(" m%2d%2d%2d", atoms[it].mbl[ia].x, atoms[it].mbl[ia].y, atoms[it].mbl[ia].z);
            if (vel) // output velocity
            {
                str += FmtCore::format(" v%20.10f%20.10f%20.10f", atoms[it].vel[ia].x, atoms[it].vel[ia].y, atoms[it].vel[ia].z);
            }
            if (nspin == 2 && magmom) // output magnetic information
            {
                str += FmtCore::format(" mag%8.4f", atom_mulliken[nat_][1]);
            }
            else if (nspin == 4 && magmom) // output magnetic information
            {
                str += FmtCore::format(" mag%8.4f%8.4f%8.4f", atom_mulliken[nat_][1], atom_mulliken[nat_][2], atom_mulliken[nat_][3]);
            }
            str += "\n";
            nat_++;
        }
    }
    std::ofstream ofs(fn.c_str());
    ofs << str;
    ofs.close();
    return;
}

void UnitCell::print_tau() const {
    ModuleBase::TITLE("UnitCell", "print_tau");
    // assert (direct || Coordinate == "Cartesian" || Coordinate == "Cartesian_angstrom"); // this line causes abort in unittest ReadAtomPositionsCACXY.
    // previously there are two if-statements, the first is `if(Coordinate == "Direct")` and the second is `if(Coordinate == "Cartesian" || Coordiante == "Cartesian_angstrom")`
    // however the Coordinate can also be value among Cartesian_angstrom_center_xy, Cartesian_angstrom_center_xz, Cartesian_angstrom_center_yz and Cartesian_angstrom_center_xyz
    // if Coordinate has value one of them, this print_tau will not print anything.
    std::regex pattern("Direct|Cartesian(_angstrom)?(_center_(xy|xz|yz|xyz))?");
    assert(std::regex_search(Coordinate, pattern));
    bool direct = (Coordinate == "Direct");
    std::string table;
    table += direct? "DIRECT COORDINATES\n": FmtCore::format("CARTESIAN COORDINATES ( UNIT = %20.12f Bohr ).\n", lat0);
    const std::string redundant_header = direct? "taud_": "tauc_";
    table += FmtCore::format("%8s%20s%20s%20s%8s%20s%20s%20s\n", "atom", "x", "y", "z", "mag", "vx", "vy", "vz");
    for(int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            const double& x = direct? atoms[it].taud[ia].x: atoms[it].tau[ia].x;
            const double& y = direct? atoms[it].taud[ia].y: atoms[it].tau[ia].y;
            const double& z = direct? atoms[it].taud[ia].z: atoms[it].tau[ia].z;
            table += FmtCore::format("%5s%-s%-5d%20.10f%20.10f%20.10f%8.4f%20.10f%20.10f%20.10f\n", // I dont know why there must be a redundant "tau[c|d]_" in the output. So ugly, it should be removed!
                                     redundant_header, atoms[it].label, ia+1, x, y, z, atoms[it].mag[ia], 
                                     atoms[it].vel[ia].x, atoms[it].vel[ia].y, atoms[it].vel[ia].z);
        }
    }
    table += "\n";
    GlobalV::ofs_running << table << std::endl;
    return;
}

/*
int UnitCell::find_type(const std::string &label)
{
    if(PARAM.inp.test_pseudo_cell) ModuleBase::TITLE("UnitCell","find_type");
    assert(ntype>0);
    for(int it=0;it<ntype;it++)
    {
        if(atoms[it].label == label)
        {
            return it;
        }
    }
    ModuleBase::WARNING_QUIT("UnitCell::find_type","Can not find the atom type!");
    return -1;
}
*/

void UnitCell::check_dtau() {
    for(int it=0; it<ntype; it++)
    {
        Atom* atom1 = &atoms[it];
        for(int ia=0; ia<atoms[it].na; ia++)
        {
            double dx2 = (atom1->taud[ia].x+10000) - int(atom1->taud[ia].x+10000);
            double dy2 = (atom1->taud[ia].y+10000) - int(atom1->taud[ia].y+10000);
            double dz2 = (atom1->taud[ia].z+10000) - int(atom1->taud[ia].z+10000);

            // mohan add 2011-04-07            
            while(dx2 >= 1) 
            {
                GlobalV::ofs_warning << " dx2 is >=1 " << std::endl;
                dx2 -= 1.0;
            }
            while(dy2 >= 1) 
            {
                GlobalV::ofs_warning << " dy2 is >=1 " << std::endl;
                dy2 -= 1.0;
            }
            while(dz2 >= 1) 
            {
                GlobalV::ofs_warning << " dz2 is >=1 " << std::endl;
                dz2 -= 1.0;
            }
            // mohan add 2011-04-07            
            while(dx2<0) 
            {
                GlobalV::ofs_warning << " dx2 is <0 " << std::endl;
                dx2 += 1.0;
            }
            while(dy2<0) 
            {
                GlobalV::ofs_warning << " dy2 is <0 " << std::endl;
                dy2 += 1.0;
            }
            while(dz2<0) 
            {
                GlobalV::ofs_warning << " dz2 is <0 " << std::endl;
                dz2 += 1.0;
            }

            atom1->taud[ia].x = dx2;
            atom1->taud[ia].y = dy2;
            atom1->taud[ia].z = dz2;

            double cx2=0.0;
            double cy2=0.0;
            double cz2=0.0;

            ModuleBase::Mathzone::Direct_to_Cartesian(
            atom1->taud[ia].x, atom1->taud[ia].y, atom1->taud[ia].z,
            latvec.e11, latvec.e12, latvec.e13,
            latvec.e21, latvec.e22, latvec.e23,
            latvec.e31, latvec.e32, latvec.e33,
            cx2, cy2, cz2);

            atom1->tau[ia].x = cx2;
            atom1->tau[ia].y = cy2;
            atom1->tau[ia].z = cz2;

    //        std::cout << std::setw(15) << dx2 << std::setw(15) << dy2 << std::setw(15) << dz2 
    //        << std::setw(15) << cx2 << std::setw(15) << cy2 << std::setw(15) << cz2
    //        << std::endl;
            
        }
    }
    return;
}

void UnitCell::read_orb_file(int it, std::string &orb_file, std::ofstream &ofs_running, Atom* atom)
{
    // the maximum L is 9 like cc-pV9Z, according to the basissetexchange https://www.basissetexchange.org/
    // there is no orbitals with L>9 presently
    const std::string spectrum = "SPDFGHIKLM";
    std::ifstream ifs(orb_file.c_str(), std::ios::in);  // pengfei 2014-10-13
    // mohan add return 2021-04-26
    if (!ifs)
    {
        std::cout << " Element index " << it+1 << std::endl;
        std::cout << " orbital file: " << orb_file << std::endl;
        ModuleBase::WARNING_QUIT("UnitCell::read_orb_file", "ABACUS Cannot find the ORBITAL file (basis sets)");
    }
    std::string word;
    atom->nw = 0;
    while (ifs.good())
    {
        ifs >> word;
        if (word == "Element")         // pengfei Li 16-2-29
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifs, atom->label_orb);
        }
        if (word == "Lmax")
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifs, atom->nwl);
            atom->l_nchi.resize(atom->nwl+1, 0);
        }
        // assert(atom->nwl<10); // cannot understand why restrict the maximum value of atom->nwl
        if (word == "Cutoff(a.u.)")         // pengfei Li 16-2-29
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifs, atom->Rcut);
        }
        if (FmtCore::endswith(word, "orbital-->"))
        {
            bool valid = false;
            for (int i = 0; i < spectrum.size(); i++)
            {
                if (word == spectrum.substr(i, 1) + "orbital-->")
                {
                    ModuleBase::GlobalFunc::READ_VALUE(ifs, atom->l_nchi[i]);
                    atom->nw += (2*i + 1) * atom->l_nchi[i];
                    std::stringstream ss;
                    ss << "L=" << i << ", number of zeta";
                    ModuleBase::GlobalFunc::OUT(ofs_running,ss.str(),atom->l_nchi[i]);
                    valid = true;
                    break;
                }
            }
            if (!valid)
            {
                ModuleBase::WARNING_QUIT("UnitCell::read_orb_file", 
                                         "ABACUS does not support numerical atomic orbital with L > 9, "
                                         "or an invalid orbital label is found in the ORBITAL file.");
            }
        }
    }
    ifs.close();
    if(!atom->nw)
    {
        ModuleBase::WARNING_QUIT("UnitCell::read_orb_file","get nw = 0");
    }
}
