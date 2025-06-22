#include <cstring>        // Peize Lin fix bug about strcmp 2016-08-02
#include <cassert>
#include <regex>
#include <fstream>

#include "unitcell.h"
#include "module_parameter/parameter.h"
#include "source_cell/print_cell.h"
#include "source_cell/read_stru.h"
#include "module_elecstate/read_orb.h"
#include "source_base/timer.h"
#include "source_base/constants.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "source_base/formatter.h"
#include "source_base/mathzone.h"

#ifdef __LCAO
#include "source_basis/module_ao/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif

bool unitcell::read_atom_positions(UnitCell& ucell,
                         std::ifstream &ifpos, 
                         std::ofstream &ofs_running, 
                         std::ofstream &ofs_warning)
{
    ModuleBase::TITLE("UnitCell","read_atom_positions");

    std::string& Coordinate  = ucell.Coordinate;
    const int    ntype       = ucell.ntype;

    if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifpos, "ATOMIC_POSITIONS"))
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifpos, Coordinate);

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
        ucell.nat = 0;

        //======================================
        // calculate total number of ucell.atoms
        // and adjust the order of atom species
        //======================================
        for (int it = 0;it < ntype; it++)
        {
            ofs_running << "\n READING ATOM TYPE " << it+1 << std::endl;
            
            //=======================================
            // (1) read in atom label
            // start magnetization
            //=======================================
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, ucell.atoms[it].label);

            if(ucell.atoms[it].label != ucell.atom_label[it])
            {
                ofs_warning << " Label orders in ATOMIC_POSITIONS and ATOMIC_SPECIES sections do not match!" << std::endl;
                ofs_warning << " Label read from ATOMIC_POSITIONS is " << ucell.atoms[it].label << std::endl; 
                ofs_warning << " Label from ATOMIC_SPECIES is " << ucell.atom_label[it] << std::endl;
                return false;
            }
            ModuleBase::GlobalFunc::OUT(ofs_running, "Atom label", ucell.atoms[it].label);

            bool set_element_mag_zero = false;
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, ucell.magnet.start_mag[it]);

#ifndef __SYMMETRY
            //===========================================
            // (2) read in numerical orbital information
            // int ucell.atoms[it].nwl
            // int* ucell.atoms[it].l_nchi;
            //===========================================

            if ((PARAM.inp.basis_type == "lcao")||(PARAM.inp.basis_type == "lcao_in_pw"))
            {
                std::string orbital_file = PARAM.inp.orbital_dir + ucell.orbital_fn[it];
                elecstate::read_orb_file(it, orbital_file, ofs_running, &(ucell.atoms[it]));
            }
            else if(PARAM.inp.basis_type == "pw")
            {
                if ((PARAM.inp.init_wfc.substr(0, 3) == "nao") || PARAM.inp.onsite_radius > 0.0)
                {
                    std::string orbital_file = PARAM.inp.orbital_dir + ucell.orbital_fn[it];
                    elecstate::read_orb_file(it, orbital_file, ofs_running, &(ucell.atoms[it]));
                }
                else
                {
                    ucell.atoms[it].nw = 0;
                    ucell.atoms[it].nwl = 2;
                    if ( ucell.lmaxmax != 2 )
                    {
                        ucell.atoms[it].nwl = ucell.lmaxmax;
                    }
                    ucell.atoms[it].l_nchi.resize(ucell.atoms[it].nwl+1, 0);
                    for(int L=0; L<ucell.atoms[it].nwl+1; L++)
                    {
                        ucell.atoms[it].l_nchi[L] = 1;
                        // calculate the number of local basis(3D)
                        ucell.atoms[it].nw += (2*L + 1) * ucell.atoms[it].l_nchi[L];
                        std::stringstream ss;
                        ss << "L=" << L << ", number of zeta";
                        ModuleBase::GlobalFunc::OUT(ofs_running,ss.str(),ucell.atoms[it].l_nchi[L]);
                    }
                }
            } // end basis type
#endif

            //=========================
            // (3) read in atom number
            //=========================
            ModuleBase::GlobalFunc::READ_VALUE(ifpos, na);
            ucell.atoms[it].na = na;

            ModuleBase::GlobalFunc::OUT(ofs_running,"Number of atoms for this type",na);

            ucell.nat += na;

            /**
             * liuyu update 2023-05-11
             * In order to employ the DP model as esolver,
             * all atom types must be specified in the `STRU` in the order consistent with that of the DP model,
             * even if the number of ucell.atoms is zero!
             */
            if (na < 0)
            {
                ModuleBase::WARNING("read_atom_positions", " atom number < 0.");
                return false;
            }
            else if (na == 0)
            {
                std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                std::cout << " Warning: atom number is 0 for atom type: " << ucell.atoms[it].label << std::endl;
                std::cout << " If you are confident that this is not a mistake, please ignore this warning." << std::endl;
                std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                ofs_running << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
                ofs_running << " Warning: atom number is 0 for atom type: " << ucell.atoms[it].label << std::endl;
                ofs_running << " If you are confident that this is not a mistake, please ignore this warning." << std::endl;
                ofs_running << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
            }
            else if (na > 0)
            {
                ucell.atoms[it].tau.resize(na, ModuleBase::Vector3<double>(0,0,0));
                ucell.atoms[it].dis.resize(na, ModuleBase::Vector3<double>(0,0,0));
                ucell.atoms[it].taud.resize(na, ModuleBase::Vector3<double>(0,0,0));
                ucell.atoms[it].vel.resize(na, ModuleBase::Vector3<double>(0,0,0));
                ucell.atoms[it].mbl.resize(na, ModuleBase::Vector3<int>(0,0,0));
                ucell.atoms[it].mag.resize(na, 0);
                ucell.atoms[it].angle1.resize(na, 0);
                ucell.atoms[it].angle2.resize(na, 0);
                ucell.atoms[it].m_loc_.resize(na, ModuleBase::Vector3<double>(0,0,0));
                ucell.atoms[it].lambda.resize(na, ModuleBase::Vector3<double>(0,0,0));
                ucell.atoms[it].constrain.resize(na, ModuleBase::Vector3<int>(0,0,0));
                ucell.atoms[it].mass = ucell.atom_mass[it]; //mohan add 2011-11-07 
                for (int ia = 0;ia < na; ia++)
                {
                 // modify the reading of frozen ions and velocities  -- Yuanbo Li 2021/8/20
                    ifpos >> v.x >> v.y >> v.z;
                    mv.x = true ;
                    mv.y = true ;
                    mv.z = true ;
                    ucell.atoms[it].vel[ia].set(0,0,0);
                    ucell.atoms[it].mag[ia]=ucell.magnet.start_mag[it];
                    //if this line is used, default startmag_type would be 2
                    ucell.atoms[it].angle1[ia]=0;
                    ucell.atoms[it].angle2[ia]=0;
                    ucell.atoms[it].m_loc_[ia].set(0,0,0);
                    ucell.atoms[it].lambda[ia].set(0,0,0);
                    ucell.atoms[it].constrain[ia].set(0,0,0);

                    std::string tmpid;
                    tmpid = ifpos.get();

                    if( (int)tmpid[0] < 0 )
                    {
                        std::cout << "read_atom_positions, mismatch in atom number for atom type: " 
                                  << ucell.atoms[it].label << std::endl;
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
                                ifpos >> ucell.atoms[it].vel[ia].x >> ucell.atoms[it].vel[ia].y >> ucell.atoms[it].vel[ia].z;
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
                                ifpos >> ucell.atoms[it].m_loc_[ia].y>>ucell.atoms[it].m_loc_[ia].z;
                                ucell.atoms[it].m_loc_[ia].x=tmpamg;
                                ucell.atoms[it].mag[ia]=sqrt(pow(ucell.atoms[it].m_loc_[ia].x,2)
                                  +pow(ucell.atoms[it].m_loc_[ia].y,2)
                                  +pow(ucell.atoms[it].m_loc_[ia].z,2));
                                input_vec_mag=true;
                                
                            }
                            else
                            {
                                ifpos.putback(tmp);
                                ucell.atoms[it].mag[ia]=tmpamg;
                            }
                        }
                        else if ( tmpid == "angle1")
                        {
                                ifpos >> ucell.atoms[it].angle1[ia];
                                ucell.atoms[it].angle1[ia]=ucell.atoms[it].angle1[ia]/180 *ModuleBase::PI;
                                input_angle_mag=true;
                                set_element_mag_zero = true;
                        }
                        else if ( tmpid == "angle2")
                        {
                                ifpos >> ucell.atoms[it].angle2[ia];
                                ucell.atoms[it].angle2[ia]=ucell.atoms[it].angle2[ia]/180 *ModuleBase::PI;
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
                                ifpos >> ucell.atoms[it].lambda[ia].y>>ucell.atoms[it].lambda[ia].z;
                                ucell.atoms[it].lambda[ia].x=tmplam;
                            }
                            else
                            {
                                ifpos.putback(tmp);
                                ucell.atoms[it].lambda[ia].z=tmplam;
                            }
                            ucell.atoms[it].lambda[ia].x /= ModuleBase::Ry_to_eV;
                            ucell.atoms[it].lambda[ia].y /= ModuleBase::Ry_to_eV;
                            ucell.atoms[it].lambda[ia].z /= ModuleBase::Ry_to_eV;
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
                                ifpos >> ucell.atoms[it].constrain[ia].y>>ucell.atoms[it].constrain[ia].z;
                                ucell.atoms[it].constrain[ia].x=tmplam;
                            }
                            else
                            {
                                ifpos.putback(tmp);
                                ucell.atoms[it].constrain[ia].z=tmplam;
                            }
                        } 
                    }
                    // move to next line
                    while ( (tmpid != "\n") && (ifpos.good()) )
                    {
                            tmpid = ifpos.get();
                    }
                    std::string mags;

                    // ----------------------------------------------------------------------------
                    // recalcualte mag and m_loc_ from read in angle1, angle2 and mag or mx, my, mz
                    if(input_angle_mag)
                    {// angle1 or angle2 are given, calculate mx, my, mz from angle1 and angle2 and mag
                        ucell.atoms[it].m_loc_[ia].z = ucell.atoms[it].mag[ia] *
                            cos(ucell.atoms[it].angle1[ia]);
                        if(std::abs(sin(ucell.atoms[it].angle1[ia])) > 1e-10 )
                        {
                            ucell.atoms[it].m_loc_[ia].x = ucell.atoms[it].mag[ia] *
                                sin(ucell.atoms[it].angle1[ia]) * cos(ucell.atoms[it].angle2[ia]);
                            ucell.atoms[it].m_loc_[ia].y = ucell.atoms[it].mag[ia] *
                                sin(ucell.atoms[it].angle1[ia]) * sin(ucell.atoms[it].angle2[ia]);
                        }
                    }
                    else if (input_vec_mag)
                    {// mx, my, mz are given, calculate angle1 and angle2 from mx, my, mz
                        double mxy=sqrt(pow(ucell.atoms[it].m_loc_[ia].x,2)+pow(ucell.atoms[it].m_loc_[ia].y,2));
                        ucell.atoms[it].angle1[ia]=atan2(mxy,ucell.atoms[it].m_loc_[ia].z);
                        if(mxy>1e-8)
                        {
                            ucell.atoms[it].angle2[ia]=atan2(ucell.atoms[it].m_loc_[ia].y,ucell.atoms[it].m_loc_[ia].x);
                        }
                    }
                    else// only one mag is given, assume it is z
                    {
                        ucell.atoms[it].m_loc_[ia].x = 0;
                        ucell.atoms[it].m_loc_[ia].y = 0;
                        ucell.atoms[it].m_loc_[ia].z = ucell.atoms[it].mag[ia];
                    }

                    if(PARAM.inp.nspin==4)
                    {
                        if(!PARAM.inp.noncolin)
                        {
                            //collinear case with nspin = 4, only z component is used
                            ucell.atoms[it].m_loc_[ia].x = 0;
                            ucell.atoms[it].m_loc_[ia].y = 0;
                        }
                        //print only ia==0 && mag>0 to avoid too much output
                        //print when ia!=0 && mag[ia] != mag[0] to avoid too much output
                        //  'A || (!A && B)' is equivalent to 'A || B',so the following 
                        // code is equivalent to 'ia==0 || (...)'
                        if(ia==0 || (ucell.atoms[it].m_loc_[ia].x != ucell.atoms[it].m_loc_[0].x 
                                    || ucell.atoms[it].m_loc_[ia].y != ucell.atoms[it].m_loc_[0].y 
                                    || ucell.atoms[it].m_loc_[ia].z != ucell.atoms[it].m_loc_[0].z))
                        {
                            //use a stringstream to generate string: "concollinear magnetization of element it is:"
                            std::stringstream ss;
                            ss << "Magnetization for this type";
                            if(ia!=0) 
                            {
                                ss<<" (atom"<<ia+1<<")";
                            }
                            ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(),
                              ucell.atoms[it].m_loc_[ia].x, 
                              ucell.atoms[it].m_loc_[ia].y, 
                              ucell.atoms[it].m_loc_[ia].z);
                        }
                        ModuleBase::GlobalFunc::ZEROS(ucell.magnet.ux_ ,3);
                    }
                    else if(PARAM.inp.nspin==2)
                    {// collinear case with nspin = 2, only z component is used
                        ucell.atoms[it].mag[ia] = ucell.atoms[it].m_loc_[ia].z;
                        //print only ia==0 && mag>0 to avoid too much output
                        //print when ia!=0 && mag[ia] != mag[0] to avoid too much output
                        if(ia==0 || (ucell.atoms[it].mag[ia] != ucell.atoms[it].mag[0]))
                        {
                            //use a stringstream to generate string: "cocollinear magnetization of element it is:"
                            std::stringstream ss;
                            ss << "magnetization of element " << it+1;
                            if(ia!=0) 
                            {
                                ss<<" (atom"<<ia+1<<")";
                            }
                            ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(),ucell.atoms[it].mag[ia]);
                        }
                    }
                    // end of calculating initial magnetization of each atom
                    // ----------------------------------------------------------------------------
            
                    if(Coordinate=="Direct")
                    {
                        // change v from direct to cartesian,
                        // the unit is GlobalC::sf.ucell.lat0
                        ucell.atoms[it].taud[ia] = v;
                        ucell.atoms[it].tau[ia] = v * ucell.latvec;
                    }
                    else if(Coordinate=="Cartesian")
                    {
                        ucell.atoms[it].tau[ia] = v ;// in unit ucell.lat0
                    }
                    else if(Coordinate=="Cartesian_angstrom")
                    {
                        ucell.atoms[it].tau[ia] = v / 0.529177 / ucell.lat0;
                    }    
                    else if(Coordinate=="Cartesian_angstrom_center_xy")
                    {
                        // calculate lattice center 
                        ucell.latcenter.x = (ucell.latvec.e11 + ucell.latvec.e21 + ucell.latvec.e31)/2.0;
                        ucell.latcenter.y = (ucell.latvec.e12 + ucell.latvec.e22 + ucell.latvec.e32)/2.0;
                        ucell.latcenter.z = 0.0;
                        ucell.atoms[it].tau[ia] = v / 0.529177 / ucell.lat0 + ucell.latcenter; 
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_xz")
                    {
                        // calculate lattice center 
                        ucell.latcenter.x = (ucell.latvec.e11 + ucell.latvec.e21 + ucell.latvec.e31)/2.0;
                        ucell.latcenter.y = 0.0; 
                        ucell.latcenter.z = (ucell.latvec.e13 + ucell.latvec.e23 + ucell.latvec.e33)/2.0;    
                        ucell.atoms[it].tau[ia] = v / 0.529177 / ucell.lat0 + ucell.latcenter; 
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_yz")
                    {
                        // calculate lattice center 
                        ucell.latcenter.x = 0.0; 
                        ucell.latcenter.y = (ucell.latvec.e12 + ucell.latvec.e22 + ucell.latvec.e32)/2.0;
                        ucell.latcenter.z = (ucell.latvec.e13 + ucell.latvec.e23 + ucell.latvec.e33)/2.0;    
                        ucell.atoms[it].tau[ia] = v / 0.529177 / ucell.lat0 + ucell.latcenter; 
                    }
                    else if(Coordinate=="Cartesian_angstrom_center_xyz")
                    {
                        // calculate lattice center 
                        ucell.latcenter.x = (ucell.latvec.e11 + ucell.latvec.e21 + ucell.latvec.e31)/2.0;
                        ucell.latcenter.y = (ucell.latvec.e12 + ucell.latvec.e22 + ucell.latvec.e32)/2.0;
                        ucell.latcenter.z = (ucell.latvec.e13 + ucell.latvec.e23 + ucell.latvec.e33)/2.0;    
                        ucell.atoms[it].tau[ia] = v / 0.529177 / ucell.lat0 + ucell.latcenter; 
                    }
                    else if(Coordinate=="Cartesian_au")
                    {
                        ucell.atoms[it].tau[ia] = v / ucell.lat0;
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
						ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[it].tau[ia].x, 
								ucell.atoms[it].tau[ia].y, 
								ucell.atoms[it].tau[ia].z,
								ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13,
								ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23,
								ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33,
								dx,dy,dz);
                    
                        ucell.atoms[it].taud[ia].x = dx;
                        ucell.atoms[it].taud[ia].y = dy;
                        ucell.atoms[it].taud[ia].z = dz;

                    }
                    
                    if(!PARAM.inp.fixed_atoms)
                    {
                        ucell.atoms[it].mbl[ia] = mv;
                    }
                    else
                    {
                        ucell.atoms[it].mbl[ia] = 0.0;
                        ucell.atoms[it].mbl[ia].print();
                    }
                    ucell.atoms[it].dis[ia].set(0, 0, 0);
                }//endj
            }    // end na
            // reset some useless parameters
            if (set_element_mag_zero)
            {
                ucell.magnet.start_mag[it] = 0.0;
            }
        } // end for ntype

        // Start Autoset magnetization
        // defaultly set a finite magnetization if magnetization is not specified
        int autoset_mag = 1;
        for (int it = 0;it < ntype; it++)
        {
            for (int ia = 0;ia < ucell.atoms[it].na; ia++)
            {
                if(std::abs(ucell.atoms[it].mag[ia]) > 1e-5)
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
                    for (int ia = 0;ia < ucell.atoms[it].na; ia++)
                    {
                        ucell.atoms[it].m_loc_[ia].x = 1.0;
                        ucell.atoms[it].m_loc_[ia].y = 1.0;
                        ucell.atoms[it].m_loc_[ia].z = 1.0;
						ucell.atoms[it].mag[ia] = sqrt(pow(ucell.atoms[it].m_loc_[ia].x,2)
								+pow(ucell.atoms[it].m_loc_[ia].y,2)
								+pow(ucell.atoms[it].m_loc_[ia].z,2));
                        ModuleBase::GlobalFunc::OUT(ofs_running,"Autoset magnetism for this atom", 1.0, 1.0, 1.0);
                    }
                }
            }
            else if(PARAM.inp.nspin==2)
            {
                for (int it = 0;it < ntype; it++)
                {
                    for (int ia = 0;ia < ucell.atoms[it].na; ia++)
                    {
                        ucell.atoms[it].mag[ia] = 1.0;
                        ucell.atoms[it].m_loc_[ia].x = ucell.atoms[it].mag[ia];
                        ModuleBase::GlobalFunc::OUT(ofs_running,"Autoset magnetism for this atom", 1.0);
                    }
                }
            }
        }
        // End Autoset magnetization
    }   // end scan_begin

    //check if any atom can move in MD
    if(!ucell.if_atoms_can_move() && PARAM.inp.calculation=="md" && PARAM.inp.esolver_type!="tddft")
    {
        ModuleBase::WARNING("read_atoms", "no atoms can move in MD simulations!");
        return false;
    } 

    ofs_running << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running,"TOTAL ATOM NUMBER",ucell.nat);
    ofs_running << std::endl;

    if (ucell.nat == 0)
    {
        ModuleBase::WARNING("read_atom_positions","no atoms found in the system!");
        return false;
    }

    // mohan add 2010-06-30    
    unitcell::check_dtau(ucell.atoms,ucell.ntype, ucell.lat0, ucell.latvec);

    if (unitcell::check_tau(ucell.atoms, ucell.ntype, ucell.lat0))
    {
        unitcell::print_tau(ucell.atoms,ucell.Coordinate,ucell.ntype,ucell.lat0,ofs_running);
        return true;
    }
    return false;

}//end read_atom_positions
