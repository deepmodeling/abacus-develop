#include "read_stru.h"

#include <cmath>
#include <sstream>

#include "source_base/tool_title.h"

namespace unitcell
{
namespace
{
// SCAN_LINE_BEGIN matches the keyword as any token in a line, including
// tokens after '#'; a lattice section header must be the first token of its
// line, otherwise a comment mention of the keyword is taken as the section.
bool scan_lattice_section(std::ifstream& ifs, const std::string& name)
{
    ifs.clear();
    ifs.seekg(0);
    std::string line;
    while (std::getline(ifs, line))
    {
        std::istringstream iss(line);
        std::string first_token;
        if (iss >> first_token && first_token == name)
        {
            return true;
        }
    }
    return false;
}
} // namespace

bool read_atom_species(std::ifstream& ifa,
                      std::ofstream& ofs_running,
                      UnitCell& ucell,
                      const std::string& basis_type,
                      const std::string& orbital_dir,
                      const std::string& init_wfc,
                      const double onsite_radius,
                      const bool deepks_setorb,
                      const bool rpa)
{
    ModuleBase::TITLE("UnitCell","read_atom_species");

    const int ntype = ucell.ntype;
    std::string word;

    //==========================================
    // read in information of each type of atom
    //==========================================
    if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa, "ATOMIC_SPECIES") )
    {    
        ModuleBase::GlobalFunc::OUT(ofs_running,"Number of elements",ntype);
        for (int i = 0;i < ntype;i++)
        {
            std::string one_line;
            std::string one_string;
            std::getline(ifa, one_line);
            std::stringstream ss;
            ss << one_line;
            ss >> ucell.atoms[i].label >> ucell.atoms[i].mass;
            if (ss.fail())
            {
                ModuleBase::WARNING_QUIT("read_atom_species",
                                         "invalid ATOMIC_SPECIES line: failed to read element label and mass");
            }
            ucell.pseudo_fn[i] = "auto";
            ucell.pseudo_type[i] = "auto";

            bool end = false;
            if (ss >> one_string)
            {
                if (one_string[0] != '#')
                {
                    ucell.pseudo_fn[i] = one_string;
                }
                else
                {
                    end = true;
                }
            }

            if (!end && ss >> one_string && one_string[0] != '#')
            {
                if (one_string == "auto" || one_string == "upf" 
                    || one_string == "vwr" || one_string == "upf201" || one_string == "blps")
                {
                    ucell.pseudo_type[i] = one_string;
                }
                else if (one_string == "1/r")
                {
                    ucell.atoms[i].coulomb_potential = true;
                }
                else
                {
                    GlobalV::ofs_warning << "unrecognized pseudopotential type: " 
                    << one_string << ", check your STRU file." << std::endl;
                    ModuleBase::WARNING_QUIT("read_atom_species", "unrecognized pseudopotential type.");
                }
            }

            // Peize Lin test for bsse 2021.04.07
            const std::string bsse_label = "empty";
            ucell.atoms[i].flag_empty_element = 
                (search( ucell.atoms[i].label.begin(), ucell.atoms[i].label.end(), 
                    bsse_label.begin(), bsse_label.end() ) != ucell.atoms[i].label.end())
                ? true : false;
        }
    }

    if((basis_type == "lcao")
      ||(basis_type == "lcao_in_pw")
      ||((basis_type == "pw")&&(init_wfc.substr(0, 3) == "nao"))
      || onsite_radius > 0.0)
    {
        if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa, "NUMERICAL_ORBITAL") )
        {
            for(int i=0; i<ntype; i++)
            {
                ifa >> ucell.orbital_fn[i];
            }
        }    
        // caoyu add 2021-03-16
        if(deepks_setorb)
        {
            if (ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa, "NUMERICAL_DESCRIPTOR")) {
                ifa >> ucell.descriptor_file;
            }
        }
        else
        {
            ucell.descriptor_file = orbital_dir + ucell.orbital_fn[0];
        }
    }
    // Peize Lin add 2016-09-23
    // Read the ABFS/JLE orbital filenames (used by LCAO EXX) into the UnitCell.
    // The EXX layer copies these into the global Exx_Info during its own setup, so
    // source_cell does not depend on the XC module. Absent sections are no-ops.
    if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa, "ABFS_ORBITAL") )
    {
        for(int i=0; i<ntype; i++)
        {
            std::string ofile;
            ifa >> ofile;
            ucell.abfs_orbital_files.push_back(ofile);
        }
    }
    if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa, "ABFS_JLES_ORBITAL") )
    {
        for(int i=0; i<ntype; i++)
        {
            std::string ofile;
            ifa >> ofile;
            ucell.jle_orbital_files.push_back(ofile);
        }
    }
    return true;
}

bool read_lattice_constant(std::ifstream& ifa,
                           std::ofstream& ofs_running,
                           Lattice& lat)
{
    //==========================
    // read in lattice constant
    //==========================
    double& lat0 = lat.lat0;
    double& lat0_angstrom =lat.lat0_angstrom;
    std::string& latName = lat.latName;
    ModuleBase::Matrix3& latvec = lat.latvec;
    if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa, "LATTICE_CONSTANT") )
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifa, lat0);
        if(lat0<=0.0)
        {
            ModuleBase::WARNING_QUIT("read_atom_species","Lattice constant <= 0.0");
        }
        lat0_angstrom = lat0 * ModuleBase::BOHR_TO_A;
        ModuleBase::GlobalFunc::OUT(ofs_running,"Lattice constant (Bohr)",lat0);
        ModuleBase::GlobalFunc::OUT(ofs_running,"Lattice constant (Angstrom)",lat0_angstrom);
        lat.tpiba  = ModuleBase::TWO_PI / lat0;
        lat.tpiba2 = lat.tpiba * lat.tpiba;
    }

    //===========================
    // Read in latticies vector
    //===========================

    if(latName=="user_defined_lattice")
    {
        // check the existence of keyword "LATTICE_PARAMETERS"
        if (scan_lattice_section(ifa, "LATTICE_PARAMETERS"))
        {
            ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
            "do not use LATTICE_PARAMETERS without explicit specification of lattice type");
        }

        // check the existence of keyword "LATTICE_VECTORS"
        if( !ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa, "LATTICE_VECTORS") )
        {
            ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
            "Please set LATTICE_VECTORS in the STRU file");
        }
        else if( ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa, "LATTICE_VECTORS") )
        {
            // Reading lattice vectors. notice
            // here that only one cpu read these
            // parameters.
            std::string latvec_line_rest;
            ifa >> latvec.e11 >> latvec.e12 >> latvec.e13;
            std::getline(ifa, latvec_line_rest); // swallow inline comments after each row
            ifa >> latvec.e21 >> latvec.e22 >> latvec.e23;
            std::getline(ifa, latvec_line_rest);
            ifa >> latvec.e31 >> latvec.e32 >> latvec.e33;
            if (ifa.fail())
            {
                ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                         "invalid LATTICE_VECTORS: failed to read 9 lattice vector components");
            }
            std::getline(ifa, latvec_line_rest);
        }
    }//supply lattice vectors
    else
    {
        if (ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(ifa,
                                               "LATTICE_VECTORS",
                                               true,
                                               false)) 
        {
            ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
            "do not use LATTICE_VECTORS along with explicit specification of lattice type");
        }
        if(latName=="sc")
        {//simple-cubic, ibrav = 1
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 1.0;
        }
        else if(latName=="fcc")
        {//face-centered cubic, ibrav = 2
            latvec.e11 =-0.5; latvec.e12 = 0.0; latvec.e13 = 0.5;
            latvec.e21 = 0.0; latvec.e22 = 0.5; latvec.e23 = 0.5;
            latvec.e31 =-0.5; latvec.e32 = 0.5; latvec.e33 = 0.0;
        }
        else if(latName=="bcc")
        {//body-centered cubic, ibrav = 3
            latvec.e11 = 0.5; latvec.e12 = 0.5; latvec.e13 = 0.5;
            latvec.e21 =-0.5; latvec.e22 = 0.5; latvec.e23 = 0.5;
            latvec.e31 =-0.5; latvec.e32 =-0.5; latvec.e33 = 0.5;
        }
        else if(latName=="hexagonal")
        {//hexagonal, ibrav = 4
            double e22 = sqrt(3.0) / 2.0;
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 =-0.5; latvec.e22 = e22; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> latvec.e33;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type hexagonal");
                }
            }
        }
        else if(latName=="trigonal")
        {//trigonal, ibrav = 5
            double t1 = 0.0;
            double t2 = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                double cosab=0.0;
                // READ_VALUE would hide an extraction failure here: a failed
                // read leaves cosab = 0, silently building the simple-cubic
                // (alpha = 90 degree) special case instead of the intended
                // trigonal cell.
                ifa >> cosab;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid cos(alpha) in LATTICE_PARAMETERS for lattice type trigonal");
                }
                t1 = sqrt(1.0 + 2.0*cosab);
                t2 = sqrt(1.0 - cosab);
            }
            double e11 = t2 / sqrt(2.0);
            double e12 = -t2 / sqrt(6.0);
            double e13 = t1 / sqrt(3.0);
            double e22 = sqrt(2.0) * t2 / sqrt(3.0);
            latvec.e11 = e11; latvec.e12 = e12; latvec.e13 = e13;
            latvec.e21 = 0.0; latvec.e22 = e22; latvec.e23 = e13;
            latvec.e31 =-e11; latvec.e32 = e12; latvec.e33 = e13;
        }
        else if(latName=="st")
        {//simple tetragonal, ibrav= 6
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> latvec.e33;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type st");
                }
            }
        }
        else if(latName=="bct")
        {//body-centered tetragonal, ibrav = 7
            double cba = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> cba;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type bct");
                }
                cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 =-0.5; latvec.e13 = cba;
            latvec.e21 = 0.5; latvec.e22 = 0.5; latvec.e23 = cba;
            latvec.e31 =-0.5; latvec.e32 =-0.5; latvec.e33 = cba;
        }
        else if(latName=="so")
        {//simple orthorhombic, ibrav = 8
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = 0.0; latvec.e22 = 0.0; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> latvec.e22 >> latvec.e33;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type so");
                }
            }
        }
        else if(latName=="baco")
        {//base-centered orthorhombic, ibrav = 9
            latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 =-0.5; latvec.e22 = 0.0; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> latvec.e12 >> latvec.e33;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type baco");
                }
                latvec.e12 = latvec.e12 / 2.0;
                latvec.e22 = latvec.e12;
            }
        }
        else if(latName=="fco")
        {//face-centered orthorhombic, ibrav = 10
            double bba = 0.0; double cba = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type fco");
                }
                bba = bba / 2.0; cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 = cba;
            latvec.e21 = 0.5; latvec.e22 = bba; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = bba; latvec.e33 = cba;
        }
        else if(latName=="bco")
        {//body-centered orthorhombic, ibrav = 11
            double bba = 0.0; double cba = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type bco");
                }
                bba = bba / 2.0; cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 = bba; latvec.e13 = cba;
            latvec.e21 =-0.5; latvec.e22 = bba; latvec.e23 = cba;
            latvec.e31 =-0.5; latvec.e32 =-bba; latvec.e33 = cba;
        }
        else if(latName=="sm")
        {//simple monoclinic, ibrav = 12
            double bba = 0.0; double cba = 0.0;
            double cosab = 0.0;
            double e21 = 0.0; double e22 = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba >> cosab;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type sm");
                }
                e21 = bba * cosab;
                e22 = bba * sqrt(1.0-cosab*cosab);
            }
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = e21; latvec.e22 = e22; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = cba;
        }
        else if(latName=="bacm")
        {//base-centered monoclinic, ibrav = 13
            double bba = 0.0; double cba = 0.0;
            double cosab = 0.0;
            double e21 = 0.0; double e22 = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba >> cosab;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type bacm");
                }
                e21 = bba * cosab;
                e22 = bba * sqrt(1.0-cosab*cosab);
                cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 =-cba;
            latvec.e21 = e21; latvec.e22 = e22; latvec.e23 = 0.0;
            latvec.e31 = 0.5; latvec.e32 = 0.0; latvec.e33 = cba;
        }
        else if(latName=="triclinic")
        {//triclinic, ibrav = 14
            double bba = 0.0; 
            double cba = 0.0;
            double cosab = 0.0; 
            double cosac = 0.0;
            double cosbc = 0.0; 
            double sinab = 0.0;
            double term = 0.0;
            if( scan_lattice_section(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba >> cosab >> cosac >> cosbc;
                if (ifa.fail())
                {
                    ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                             "invalid LATTICE_PARAMETERS for lattice type triclinic");
                }
                sinab = sqrt(1.0-cosab*cosab);
            }
            latvec.e11 = 1.0;         latvec.e12 = 0.0;         latvec.e13 = 0.0;
            latvec.e21 = bba * cosab; latvec.e22 = bba * sinab; latvec.e23 = 0.0;
            latvec.e31 = cba * cosac; latvec.e32 = cba * (cosbc - cosac*cosab) / sinab;
            term = 1.0 + 2.0 * cosab*cosac*cosbc - cosab*cosab - cosac*cosac - cosbc*cosbc;
            term = sqrt(term)/sinab;
            latvec.e33 = cba * term;
        }
        else
        {
            ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant","latname not supported!");
        }
        // Non-cubic named lattices cannot be constructed without the shape
        // parameters in LATTICE_PARAMETERS; without them the shape variables
        // stay zero and produce a zero-volume or NaN cell. The check comes
        // after the branches so an unsupported latName reports itself first.
        if(latName != "sc" && latName != "fcc" && latName != "bcc")
        {
            const bool has_params = scan_lattice_section(ifa, "LATTICE_PARAMETERS");
            if (!has_params)
            {
                const std::string msg = "LATTICE_PARAMETERS is required for lattice type " + latName;
                ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant", msg);
            }
        }
    }

    // Reject degenerate (zero-volume) or non-finite lattice vectors, whether
    // they came from LATTICE_VECTORS or from a named lattice.
    const double latvec_det = latvec.Det();
    if (!std::isfinite(latvec_det) || latvec_det == 0.0)
    {
        ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant",
                                 "the lattice vectors form a zero-volume or non-finite cell");
    }

    // lattice vectors in another form.
    lat.a1.x = latvec.e11;
    lat.a1.y = latvec.e12;
    lat.a1.z = latvec.e13;

    lat.a2.x = latvec.e21;
    lat.a2.y = latvec.e22;
    lat.a2.z = latvec.e23;

    lat.a3.x = latvec.e31;
    lat.a3.y = latvec.e32;
    lat.a3.z = latvec.e33;
    return true;
}

}
