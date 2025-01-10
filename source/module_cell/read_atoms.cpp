#include "unitcell.h"
#include "module_parameter/parameter.h"
#include "read_stru.h"
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
                if ((PARAM.inp.init_wfc.substr(0, 3) == "nao") || PARAM.inp.onsite_radius > 0.0)
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
            else if (na == 0)
            {
                std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                std::cout << " Warning: atom number is 0 for atom type: " << atoms[it].label << std::endl;
                std::cout << " If you are confident that this is not a mistake, please ignore this warning." << std::endl;
                std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                ofs_running << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                ofs_running << " Warning: atom number is 0 for atom type: " << atoms[it].label << std::endl;
                ofs_running << " If you are confident that this is not a mistake, please ignore this warning." << std::endl;
                ofs_running << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
            }
            else if (na > 0)
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

                    // ----------------------------------------------------------------------------
                    // recalcualte mag and m_loc_ from read in angle1, angle2 and mag or mx, my, mz
                    if(input_angle_mag)
                    {// angle1 or angle2 are given, calculate mx, my, mz from angle1 and angle2 and mag
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
                    {// mx, my, mz are given, calculate angle1 and angle2 from mx, my, mz
                        double mxy=sqrt(pow(atoms[it].m_loc_[ia].x,2)+pow(atoms[it].m_loc_[ia].y,2));
                        atoms[it].angle1[ia]=atan2(mxy,atoms[it].m_loc_[ia].z);
                        if(mxy>1e-8)
                        {
                            atoms[it].angle2[ia]=atan2(atoms[it].m_loc_[ia].y,atoms[it].m_loc_[ia].x);
                        }
                    }
                    else// only one mag is given, assume it is z
                    {
                        atoms[it].m_loc_[ia].x = 0;
                        atoms[it].m_loc_[ia].y = 0;
                        atoms[it].m_loc_[ia].z = atoms[it].mag[ia];
                    }

                    if(PARAM.inp.nspin==4)
                    {
                        if(!PARAM.inp.noncolin)
                        {
                            //collinear case with nspin = 4, only z component is used
                            atoms[it].m_loc_[ia].x = 0;
                            atoms[it].m_loc_[ia].y = 0;
                        }
                        //print only ia==0 && mag>0 to avoid too much output
                        //print when ia!=0 && mag[ia] != mag[0] to avoid too much output
                        //  'A || (!A && B)' is equivalent to 'A || B',so the following 
                        // code is equivalent to 'ia==0 || (...)'
                        if(ia==0 || (atoms[it].m_loc_[ia].x != atoms[it].m_loc_[0].x 
                                    || atoms[it].m_loc_[ia].y != atoms[it].m_loc_[0].y 
                                    || atoms[it].m_loc_[ia].z != atoms[it].m_loc_[0].z))
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
                    {// collinear case with nspin = 2, only z component is used
                        atoms[it].mag[ia] = atoms[it].m_loc_[ia].z;
                        //print only ia==0 && mag>0 to avoid too much output
                        //print when ia!=0 && mag[ia] != mag[0] to avoid too much output
                        if(ia==0 || (atoms[it].mag[ia] != atoms[it].mag[0]))
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
                    // end of calculating initial magnetization of each atom
                    // ----------------------------------------------------------------------------
            
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

    if (nat == 0)
    {
        ModuleBase::WARNING("read_atom_positions","no atom in the system!");
        return false;
    }

    // mohan add 2010-06-30    
    //xiaohui modify 2015-03-15, cancel outputfile "STRU_READIN.xyz"
    //this->print_cell_xyz("STRU_READIN.xyz");
    this->check_dtau();

    if (unitcell::check_tau(this->atoms, this->ntype, this->lat0))
    {
        this->print_tau();
    //xiaohui modify 2015-03-15, cancel outputfile "STRU_READIN.xyz"
    //this->print_cell_xyz("STRU_READIN_ADJUST.xyz");
        return true;
    }
    return false;

}//end read_atom_positions


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
