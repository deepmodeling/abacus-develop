#include <cstring>        // Peize Lin fix bug about strcmp 2016-08-02
#include <cassert>
#include <regex>
#include <fstream>

#include "unitcell.h"
#include "read_atoms_helper.h"
#include "source_io/module_parameter/parameter.h"
#include "print_cell.h"
#include "read_stru.h"
#include "source_estate/read_orb.h"
#include "source_base/timer.h"
#include "source_base/constants.h"
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
    const int ntype = ucell.ntype;
    const int nspin = PARAM.inp.nspin;
    assert (nspin==1 || nspin==2 || nspin==4);

    if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifpos, "ATOMIC_POSITIONS"))
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifpos, Coordinate);

        if (!unitcell::validate_coordinate_system(Coordinate, ofs_warning))
        {
            return false;
        }

        ucell.nat = 0;

        //======================================
        // calculate total number of ucell.atoms
        // and adjust the order of atom species
        //======================================
        for (int it = 0;it < ntype; it++)
        {
            ofs_running << "\n READING ATOM TYPE " << it+1 << std::endl;

            bool set_element_mag_zero = false;
            if (!unitcell::read_atom_type_header(it, ucell, ifpos, ofs_running,
                                       ofs_warning, set_element_mag_zero))
            {
                return false;
            }

            int na = ucell.atoms[it].na;
            ucell.nat += na;

            if (na > 0)
            {
                unitcell::allocate_atom_properties(ucell.atoms[it], na, ucell.atom_mass[it]);
                for (int ia = 0;ia < na; ia++)
                {
                 // modify the reading of frozen ions and velocities  -- Yuanbo Li 2021/8/20
                    ModuleBase::Vector3<double> v;
                    ModuleBase::Vector3<int> mv;
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

                    bool input_vec_mag=false;
                    bool input_angle_mag=false;

                    // Parse optional properties
                    if (!unitcell::parse_atom_properties(ifpos, ucell.atoms[it], ia, mv,
                                              input_vec_mag, input_angle_mag,
                                              set_element_mag_zero))
                    {
                        return false;
                    }

                    // Process magnetization
                    unitcell::process_magnetization(ucell.atoms[it], it, ia, nspin,
                                        input_vec_mag, input_angle_mag, ofs_running);

                    // Transform coordinates
                    unitcell::transform_atom_coordinates(ucell.atoms[it], ia, Coordinate,
                                             v, ucell.latvec, ucell.lat0, ucell.latcenter);

                    // Set movement flags
                    unitcell::set_atom_movement_flags(ucell.atoms[it], ia, mv);
                    ucell.atoms[it].dis[ia].set(0, 0, 0);
                }//endj
            }    // end na
            // reset some useless parameters
            if (set_element_mag_zero)
            {
                ucell.magnet.start_mag[it] = 0.0;
            }
        } // end for ntype

        // Auto-set magnetization if needed
        unitcell::autoset_magnetization(ucell, nspin, ofs_running);
    }   // end scan_begin

    // Final validation and output
    return unitcell::finalize_atom_positions(ucell, ofs_running, ofs_warning);

}//end read_atom_positions
