#include <regex>
#include <cassert>

#include "print_cell.h"
#include "module_base/formatter.h"
#include "module_base/tool_title.h"
#include "module_base/global_variable.h"

namespace unitcell
{
    void print_tau(Atom* atoms,
                   const std::string& Coordinate,
                   const int ntype,
                   const double lat0)
    {
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
}