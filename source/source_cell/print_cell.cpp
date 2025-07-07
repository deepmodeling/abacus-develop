#include <regex>
#include <cassert>

#include "print_cell.h"
#include "source_base/formatter.h"
#include "source_base/tool_title.h"
#include "source_base/global_variable.h"

namespace unitcell
{
    void print_tau(Atom* atoms,
                   const std::string& Coordinate,
                   const int ntype,
                   const double lat0,
                   std::ofstream &ofs)
    {
        ModuleBase::TITLE("UnitCell", "print_tau");
        // assert (direct || Coordinate == "Cartesian" || Coordinate == "Cartesian_angstrom"); // this line causes abort in unittest ReadAtomPositionsCACXY.
        // previously there are two if-statements, the first is `if(Coordinate == "Direct")` and the second is `if(Coordinate == "Cartesian" || Coordiante == "Cartesian_angstrom")`
        // however the Coordinate can also be value among Cartesian_angstrom_center_xy, Cartesian_angstrom_center_xz, Cartesian_angstrom_center_yz and Cartesian_angstrom_center_xyz

        // if Coordinate has value one of them, this print_tau will not print anything.
        std::regex pattern("Direct|Cartesian(_angstrom)?(_center_(xy|xz|yz|xyz))?");
        assert(std::regex_search(Coordinate, pattern));
        bool direct = (Coordinate == "Direct");

        //----------------------
        // print atom positions
        //----------------------
        std::string table;
        table += direct? " DIRECT COORDINATES\n": FmtCore::format(" CARTESIAN COORDINATES ( UNIT = %15.8f Bohr )\n", lat0);
        table += FmtCore::format("%5s%19s%19s%19s%8s\n", "atom", "x", "y", "z", "mag");
        for(int it = 0; it < ntype; it++)
        {
            for (int ia = 0; ia < atoms[it].na; ia++)
            {
                const double& x = direct? atoms[it].taud[ia].x: atoms[it].tau[ia].x;
                const double& y = direct? atoms[it].taud[ia].y: atoms[it].tau[ia].y;
                const double& z = direct? atoms[it].taud[ia].z: atoms[it].tau[ia].z;
                table += FmtCore::format("%5s%19.12f%19.12f%19.12f%8.4f\n", 
                                        atoms[it].label, 
                                        x, 
                                        y, 
                                        z, 
                                        atoms[it].mag[ia]); 
            }
        }
        table += "\n";
        ofs << table; 


        // print velocities
        ofs << " ATOMIC VELOCITIES" << std::endl;
        ofs << std::setprecision(12);
        ofs << std::setw(5) << "atom" 
            << std::setw(19) << "vx" 
            << std::setw(19) << "vy" 
            << std::setw(19) << "vz"
            << std::endl;
 
		for(int it = 0; it < ntype; it++)
		{
			for (int ia = 0; ia < atoms[it].na; ia++)
			{
                ofs << std::setw(5) << atoms[it].label;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].x;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].y;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].z;
                ofs << std::endl;
			}
		}
        ofs << std::endl;
        ofs << std::setprecision(6); // return to 6, as original


        return;
    }

    void print_stru_file(const UnitCell& ucell,
                         const Atom*     atoms,
                         const ModuleBase::Matrix3& latvec,
                         const std::string& fn, 
                         const int& nspin,
                         const bool& direct,
                         const bool& vel,
                         const bool& magmom,
                         const bool& orb,
                         const bool& dpks_desc,
                         const int& iproc)
    {
        ModuleBase::TITLE("UnitCell","print_stru_file");
        if (iproc != 0) 
        {
            return; // old: if(GlobalV::MY_RANK != 0) return;
        }
        // ATOMIC_SPECIES
        std::string str = "ATOMIC_SPECIES\n";
        for(int it=0; it<ucell.ntype; it++)
        { 
            str += FmtCore::format("%s %8.4f %s %s\n", 
                                    ucell.atom_label[it], 
                                    ucell.atom_mass[it], 
                                    ucell.pseudo_fn[it], 
                                    ucell.pseudo_type[it]); 
        }
        // NUMERICAL_ORBITAL
        if(orb)
        {
            str += "\nNUMERICAL_ORBITAL\n";
            for(int it = 0; it < ucell.ntype; it++) 
            { 
                str += ucell.orbital_fn[it] + "\n"; 
            }
        }
        // NUMERICAL_DESCRIPTOR
        if(dpks_desc) 
        { 
            str += "\nNUMERICAL_DESCRIPTOR\n" + ucell.descriptor_file + "\n"; 
        }
        // LATTICE_CONSTANT
        str += "\nLATTICE_CONSTANT\n" + FmtCore::format("%-.10f\n", ucell.lat0);
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
        for(int it = 0; it < ucell.ntype; it++)
        {
            str += "\n" + ucell.atoms[it].label + " #label\n";
            str += FmtCore::format("%-8.4f #magnetism\n", ucell.magnet.start_mag[it]);
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
                    str += FmtCore::format(" mag%8.4f", ucell.atom_mulliken[nat_][1]);
                }
                else if (nspin == 4 && magmom) // output magnetic information
                {
                    str += FmtCore::format(" mag%8.4f%8.4f%8.4f", 
                                            ucell.atom_mulliken[nat_][1], 
                                            ucell.atom_mulliken[nat_][2], 
                                            ucell.atom_mulliken[nat_][3]);
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
}
