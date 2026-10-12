#include "read_atoms_helper.h"
#include "source_base/global_function.h"
#include "source_base/constants.h"
#include "source_base/mathzone.h"
#include "read_stru.h"
#include "print_cell.h"
#include "read_orb.h"
#include "cell_tools.h"
#include <cmath>
#include <iostream>
#include <sstream>

namespace {
    // Magic number constants for character code checks
    constexpr char DIGIT_START = '0';  // ASCII 48
    constexpr char DIGIT_END = '9';    // ASCII 57
    constexpr char LOWER_A = 'a';
    constexpr char LOWER_Z = 'z';
    constexpr char MINUS_SIGN = '-';

    // Below this value a magnetic moment component is treated as zero.
    constexpr double mag_threshold = 1e-5;
}

namespace unitcell {

bool validate_coordinate_system(const std::string& Coordinate,
                                std::ofstream& ofs_warning)
{
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
        return false;
    }
    return true;
}

void allocate_atom_properties(Atom& atom, int na)
{
    atom.tau.resize(na, ModuleBase::Vector3<double>(0,0,0));
    atom.dis.resize(na, ModuleBase::Vector3<double>(0,0,0));
    atom.taud.resize(na, ModuleBase::Vector3<double>(0,0,0));
    atom.boundary_shift.resize(na, ModuleBase::Vector3<int>(0,0,0));
    atom.vel.resize(na, ModuleBase::Vector3<double>(0,0,0));
    atom.mbl.resize(na, ModuleBase::Vector3<int>(0,0,0));
    atom.mag.resize(na, 0);
    atom.angle1.resize(na, 0);
    atom.angle2.resize(na, 0);
    atom.m_loc_.resize(na, ModuleBase::Vector3<double>(0,0,0));
    atom.lambda.resize(na, ModuleBase::Vector3<double>(0,0,0));
    atom.constrain.resize(na, ModuleBase::Vector3<int>(0,0,0));
}

void set_atom_movement_flags(Atom& atom, int ia,
                             const ModuleBase::Vector3<int>& mv,
                             const bool fixed_atoms)
{
    if(!fixed_atoms)
    {
        atom.mbl[ia] = mv;
    }
    else
    {
        atom.mbl[ia] = 0.0;
        atom.mbl[ia].print();
    }
}

bool is_magnetization_all_zero(const UnitCell& ucell)
{
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            if (std::abs(ucell.atoms[it].mag[ia]) > mag_threshold)
            {
                return false;
            }
        }
    }
    return true;
}

void autoset_magnetization(UnitCell& ucell, int nspin,
                           std::ofstream& ofs_running)
{
    // nspin=4 never autosets a moment. An all-zero moment in STRU is a
    // legitimate nonmagnetic start; warn_zero_magnetization_nspin4() reports
    // it instead, independently of the symmetry switch that gates this call.
    if (nspin != 2)
    {
        return;
    }

    if (!is_magnetization_all_zero(ucell))
    {
        return;
    }

    const int ntype = ucell.ntype;
    for (int it = 0; it < ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            ucell.atoms[it].mag[ia] = 1.0;
            ucell.atoms[it].m_loc_[ia].z = ucell.atoms[it].mag[ia];
            ModuleBase::GlobalFunc::OUT(ofs_running,"Autoset magnetism for this atom", 1.0);
        }
    }
}

void warn_zero_magnetization_nspin4(const UnitCell& ucell,
                                    std::ofstream& ofs_running)
{
    const double threshold_sq = mag_threshold * mag_threshold;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const ModuleBase::Vector3<double>& m = ucell.atoms[it].m_loc_[ia];
            const double norm_sq = m.x * m.x + m.y * m.y + m.z * m.z;
            if (norm_sq > threshold_sq)
            {
                return;
            }
        }
    }

    const std::string summary = "nspin=4 but no initial magnetization is set in STRU; "
        "all atoms start from zero magnetic moment.";
    std::cout << " Warning: " << summary << std::endl;
    ModuleBase::WARNING("read_atom_positions", summary);
    ofs_running << "\n WARNING: " << summary << std::endl;
    ofs_running << "          If a magnetic ground state is expected, set 'mag' explicitly "
                << "in STRU for the magnetic atoms." << std::endl;
    ofs_running << "          With noncolin=0 the x/y components of 'mag' are discarded, so a "
                << "moment given purely along x/y also leaves a zero starting moment."
                << std::endl;
}

void warn_xy_magnetization_ignored(const std::string& label,
                                   int n_discarded,
                                   int na,
                                   std::ofstream& ofs_running)
{
    if (n_discarded <= 0)
    {
        return;
    }

    std::stringstream ss;
    ss << "atom type " << label << ": " << n_discarded << " of " << na
       << " atoms give non-zero x/y magnetization in STRU, but nspin=4 with noncolin=0 "
       << "is a collinear calculation; only the z component is used. "
       << "x/y components are IGNORED. Set 'noncolin 1' to use the full vector.";
    const std::string summary = ss.str();
    std::cout << " Warning: " << summary << std::endl;
    ModuleBase::WARNING("read_atom_positions", summary);
    ofs_running << "\n WARNING: " << summary << std::endl;
}

bool finalize_atom_positions(UnitCell& ucell,
                             std::ofstream& ofs_running,
                             std::ofstream& ofs_warning,
                             const std::string& calculation,
                             const std::string& esolver_type)
{
    // Check if any atom can move in MD
    if(!unitcell::if_atoms_can_move(ucell.atoms, ucell.ntype) && calculation=="md" && esolver_type!="tddft")
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

    // Check atom positions
    unitcell::check_dtau(ucell.atoms, ucell.ntype, ucell.lat0, ucell.latvec);

    if (unitcell::check_tau(ucell.atoms, ucell.ntype, ucell.lat0))
    {
        print_tau(ucell.atoms, ucell.Coordinate, ucell.ntype, ucell.lat0, ofs_running);
        return true;
    }
    return false;
}

ModuleBase::Vector3<double> calculate_lattice_center(
    const ModuleBase::Matrix3& latvec,
    const std::string& center_mode)
{
    ModuleBase::Vector3<double> latcenter(0.0, 0.0, 0.0);

    if (center_mode == "xy" || center_mode == "xyz")
    {
        latcenter.x = (latvec.e11 + latvec.e21 + latvec.e31) / 2.0;
        latcenter.y = (latvec.e12 + latvec.e22 + latvec.e32) / 2.0;
    }

    if (center_mode == "xz" || center_mode == "xyz")
    {
        latcenter.x = (latvec.e11 + latvec.e21 + latvec.e31) / 2.0;
        latcenter.z = (latvec.e13 + latvec.e23 + latvec.e33) / 2.0;
    }

    if (center_mode == "yz" || center_mode == "xyz")
    {
        latcenter.y = (latvec.e12 + latvec.e22 + latvec.e32) / 2.0;
        latcenter.z = (latvec.e13 + latvec.e23 + latvec.e33) / 2.0;
    }

    return latcenter;
}

void transform_atom_coordinates(Atom& atom, int ia,
                               const std::string& Coordinate,
                               const ModuleBase::Vector3<double>& v,
                               const ModuleBase::Matrix3& latvec,
                               double lat0,
                               ModuleBase::Vector3<double>& latcenter)
{
    if(Coordinate=="Direct")
    {
        // change v from direct to cartesian,
        // the unit is ucell.lat0
        atom.taud[ia] = v;
        atom.tau[ia] = v * latvec;
    }
    else if(Coordinate=="Cartesian")
    {
        atom.tau[ia] = v;  // in unit ucell.lat0
    }
    else if(Coordinate=="Cartesian_angstrom")
    {
        atom.tau[ia] = v / ModuleBase::BOHR_TO_A / lat0;
    }
    else if(Coordinate=="Cartesian_angstrom_center_xy")
    {
        latcenter = calculate_lattice_center(latvec, "xy");
        atom.tau[ia] = v / ModuleBase::BOHR_TO_A / lat0 + latcenter;
    }
    else if(Coordinate=="Cartesian_angstrom_center_xz")
    {
        latcenter = calculate_lattice_center(latvec, "xz");
        atom.tau[ia] = v / ModuleBase::BOHR_TO_A / lat0 + latcenter;
    }
    else if(Coordinate=="Cartesian_angstrom_center_yz")
    {
        latcenter = calculate_lattice_center(latvec, "yz");
        atom.tau[ia] = v / ModuleBase::BOHR_TO_A / lat0 + latcenter;
    }
    else if(Coordinate=="Cartesian_angstrom_center_xyz")
    {
        latcenter = calculate_lattice_center(latvec, "xyz");
        atom.tau[ia] = v / ModuleBase::BOHR_TO_A / lat0 + latcenter;
    }
    else if(Coordinate=="Cartesian_au")
    {
        atom.tau[ia] = v / lat0;
    }

    // Convert to direct coordinates if using Cartesian
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
        ModuleBase::Mathzone::Cartesian_to_Direct(atom.tau[ia].x,
                atom.tau[ia].y,
                atom.tau[ia].z,
                latvec.e11, latvec.e12, latvec.e13,
                latvec.e21, latvec.e22, latvec.e23,
                latvec.e31, latvec.e32, latvec.e33,
                dx,dy,dz);

        atom.taud[ia].x = dx;
        atom.taud[ia].y = dy;
        atom.taud[ia].z = dz;
    }
}

bool process_magnetization(Atom& atom, int it, int ia,
                          int nspin, bool input_vec_mag,
                          bool input_angle_mag,
                          std::ofstream& ofs_running,
                          const bool noncolin)
{
    bool xy_discarded = false;

    // Recalculate mag and m_loc_ from read in angle1, angle2 and mag or mx, my, mz
    if(input_angle_mag)
    {
        // angle1 or angle2 are given, calculate mx, my, mz from angle1 and angle2 and mag
        atom.m_loc_[ia].z = atom.mag[ia] * cos(atom.angle1[ia]);
        if(std::abs(sin(atom.angle1[ia])) > 1e-10)
        {
            atom.m_loc_[ia].x = atom.mag[ia] *
                sin(atom.angle1[ia]) * cos(atom.angle2[ia]);
            atom.m_loc_[ia].y = atom.mag[ia] *
                sin(atom.angle1[ia]) * sin(atom.angle2[ia]);
        }
    }
    else if (input_vec_mag)
    {
        // mx, my, mz are given, calculate angle1 and angle2 from mx, my, mz
        double mxy=sqrt(pow(atom.m_loc_[ia].x,2)+pow(atom.m_loc_[ia].y,2));
        atom.angle1[ia]=atan2(mxy,atom.m_loc_[ia].z);
        if(mxy>1e-8)
        {
            atom.angle2[ia]=atan2(atom.m_loc_[ia].y,atom.m_loc_[ia].x);
        }
    }
    else
    {
        // only one mag is given, assume it is z
        atom.m_loc_[ia].x = 0;
        atom.m_loc_[ia].y = 0;
        atom.m_loc_[ia].z = atom.mag[ia];
    }

    if(nspin==4)
    {
        if(!noncolin)
        {
            // collinear case with nspin = 4, only z component is used
            // Note: a scalar "mag <value>" in STRU is stored purely in
            // m_loc_.z (x/y stay zero), so it passes the check below without
            // flagging a discard and still defines a non-zero z moment.
            // cal_ux() then turns that into the global XC quantization axis
            // ux_ = (0,0,1) with lsign_ = true. Any positive scalar suffices;
            // only the *direction* enters ux_, the magnitude (e.g. 1.7320508
            // vs 1.0) is irrelevant.
            //
            // The discard is only reported back to the caller, which
            // aggregates it per atom type: emitting a message here would
            // print one warning per atom for a cell where every atom tilts
            // the moment.
            if(std::abs(atom.m_loc_[ia].x) > mag_threshold
                || std::abs(atom.m_loc_[ia].y) > mag_threshold)
            {
                xy_discarded = true;
            }
            atom.m_loc_[ia].x = 0;
            atom.m_loc_[ia].y = 0;
        }
        // print only ia==0 && mag>0 to avoid too much output
        // print when ia!=0 && mag[ia] != mag[0] to avoid too much output
        if(ia==0 || (atom.m_loc_[ia].x != atom.m_loc_[0].x
                    || atom.m_loc_[ia].y != atom.m_loc_[0].y
                    || atom.m_loc_[ia].z != atom.m_loc_[0].z))
        {
            std::stringstream ss;
            ss << "Magnetization for this type";
            if(ia!=0)
            {
                ss<<" (atom"<<ia+1<<")";
            }
            ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(),
              atom.m_loc_[ia].x,
              atom.m_loc_[ia].y,
              atom.m_loc_[ia].z);
        }
        // Note: The original code had ZEROS(ucell.magnet.ux_, 3) here
        // but ucell is not available in this function scope
    }
    else if(nspin==2)
    {
        // collinear case with nspin = 2, only z component is used
        atom.mag[ia] = atom.m_loc_[ia].z;
        // print only ia==0 && mag>0 to avoid too much output
        // print when ia!=0 && mag[ia] != mag[0] to avoid too much output
        if(ia==0 || (atom.mag[ia] != atom.mag[0]))
        {
            std::stringstream ss;
            ss << "magnetization of element " << it+1;
            if(ia!=0)
            {
                ss<<" (atom"<<ia+1<<")";
            }
            ModuleBase::GlobalFunc::OUT(ofs_running, ss.str(),atom.mag[ia]);
        }
    }

    return xy_discarded;
}

bool parse_atom_properties(std::ifstream& ifpos,
                          Atom& atom, int ia,
                          ModuleBase::Vector3<int>& mv,
                          bool& input_vec_mag,
                          bool& input_angle_mag,
                          bool& set_element_mag_zero)
{
    std::string tmpid;
    tmpid = ifpos.get();

    if( (int)tmpid[0] < 0 )
    {
        std::cout << "read_atom_positions, mismatch in atom number for atom type: "
                  << atom.label << std::endl;
        exit(1);
    }

    // read if catch goodbit before "\n" and "#"
    while ( (tmpid != "\n") && (ifpos.good()) && (tmpid !="#") )
    {
        tmpid = ifpos.get();
        // old method of reading frozen ions
        char tmp = (char)tmpid[0];
        if ( tmp >= DIGIT_START && tmp <= DIGIT_END )
        {
            mv.x = std::stoi(tmpid);
            ifpos >> mv.y >> mv.z;
        }
        // new method of reading frozen ions and velocities
        if ( tmp >= LOWER_A && tmp <= LOWER_Z)
        {
            ifpos.putback(tmp);
            ifpos >> tmpid;
        }
        if ( tmpid == "m" )
        {
            ifpos >> mv.x >> mv.y >> mv.z;
        }
        else if ( tmpid == "v" ||tmpid == "vel" || tmpid == "velocity" )
        {
            ifpos >> atom.vel[ia].x >> atom.vel[ia].y >> atom.vel[ia].z;
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

            if((tmp >= DIGIT_START && tmp <= DIGIT_END) or tmp==MINUS_SIGN)
            {
                ifpos.putback(tmp);
                ifpos >> atom.m_loc_[ia].y>>atom.m_loc_[ia].z;
                atom.m_loc_[ia].x=tmpamg;
                atom.mag[ia]=sqrt(pow(atom.m_loc_[ia].x,2)
                  +pow(atom.m_loc_[ia].y,2)
                  +pow(atom.m_loc_[ia].z,2));
                input_vec_mag=true;

            }
            else
            {
                ifpos.putback(tmp);
                atom.mag[ia]=tmpamg;
            }
        }
        else if ( tmpid == "angle1")
        {
            ifpos >> atom.angle1[ia];
            atom.angle1[ia]=atom.angle1[ia]/180 *ModuleBase::PI;
            input_angle_mag=true;
            set_element_mag_zero = true;
        }
        else if ( tmpid == "angle2")
        {
            ifpos >> atom.angle2[ia];
            atom.angle2[ia]=atom.angle2[ia]/180 *ModuleBase::PI;
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
            if((tmp >= DIGIT_START && tmp <= DIGIT_END) or tmp==MINUS_SIGN)
            {
                ifpos.putback(tmp);
                ifpos >> atom.lambda[ia].y>>atom.lambda[ia].z;
                atom.lambda[ia].x=tmplam;
            }
            else
            {
                ifpos.putback(tmp);
                atom.lambda[ia].z=tmplam;
            }
            atom.lambda[ia].x /= ModuleBase::Ry_to_eV;
            atom.lambda[ia].y /= ModuleBase::Ry_to_eV;
            atom.lambda[ia].z /= ModuleBase::Ry_to_eV;
        }
        else if ( tmpid == "f" || tmpid == "force" || tmpid == "forces" )
        {
            // Consume the three force components emitted by print_stru_file.
            // Forces are intentionally not stored in Atom; reading them keeps
            // STRU round-trip parsing intact and prevents the legacy
            // movement-flag parser from misinterpreting force digits.
            double fx = 0.0;
            double fy = 0.0;
            double fz = 0.0;
            ifpos >> fx >> fy >> fz;
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
            if((tmp >= DIGIT_START && tmp <= DIGIT_END) or tmp==MINUS_SIGN)
            {
                ifpos.putback(tmp);
                ifpos >> atom.constrain[ia].y>>atom.constrain[ia].z;
                atom.constrain[ia].x=tmplam;
            }
            else
            {
                ifpos.putback(tmp);
                atom.constrain[ia].z=tmplam;
            }
        }
    }
    // move to next line
    while ( (tmpid != "\n") && (ifpos.good()) )
    {
        tmpid = ifpos.get();
    }

    return true;
}

bool read_atom_type_header(int it, UnitCell& ucell,
                          std::ifstream& ifpos,
                          std::ofstream& ofs_running,
                          std::ofstream& ofs_warning,
                          bool& set_element_mag_zero,
                          const std::string& basis_type,
                          const std::string& orbital_dir,
                          const std::string& init_wfc,
                          const double onsite_radius)
{
    //=======================================
    // (1) read in atom label
    // start magnetization
    //=======================================
    const std::string label_from_species = ucell.atoms[it].label;
    ModuleBase::GlobalFunc::READ_VALUE(ifpos, ucell.atoms[it].label);

    if(ucell.atoms[it].label != label_from_species)
    {
        ofs_warning << " Label orders in ATOMIC_POSITIONS and ATOMIC_SPECIES sections do not match!" << std::endl;
        ofs_warning << " Label read from ATOMIC_POSITIONS is " << ucell.atoms[it].label << std::endl;
        ofs_warning << " Label from ATOMIC_SPECIES is " << label_from_species << std::endl;
        return false;
    }
    ModuleBase::GlobalFunc::OUT(ofs_running, "Atom label", ucell.atoms[it].label);

    set_element_mag_zero = false;
    ModuleBase::GlobalFunc::READ_VALUE(ifpos, ucell.magnet.start_mag[it]);

#ifndef __SYMMETRY
    //===========================================
    // (2) read in numerical orbital information
    // int ucell.atoms[it].nwl
    // int* ucell.atoms[it].l_nchi;
    //===========================================

    if ((basis_type == "lcao")||(basis_type == "lcao_in_pw"))
    {
        std::string orbital_file = orbital_dir + ucell.orbital_fn[it];
        bool normal = unitcell::read_orb_file(it, orbital_file, ofs_running, &(ucell.atoms[it]));
        if(!normal)
        {
            return false;
        }
    }
    else if(basis_type == "pw")
    {
        if ((init_wfc.substr(0, 3) == "nao") || onsite_radius > 0.0)
        {
            std::string orbital_file = orbital_dir + ucell.orbital_fn[it];
            bool normal = unitcell::read_orb_file(it, orbital_file, ofs_running, &(ucell.atoms[it]));
            if(!normal)
            {
                return false;
            }
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
    int na = 0;
    ModuleBase::GlobalFunc::READ_VALUE(ifpos, na);
    ucell.atoms[it].na = na;

    ModuleBase::GlobalFunc::OUT(ofs_running,"Number of atoms for this type",na);

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

    return true;
}

} // namespace unitcell
