#include "read_stru.h"

#include "module_parameter/parameter.h"
#include "module_base/tool_title.h"
#ifdef __EXX
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_ri/serialization_cereal.h"
#endif
namespace unitcell
{
bool read_atom_species(std::ifstream& ifa,
                      std::ofstream& ofs_running,
                      UnitCell& ucell)
{
    ModuleBase::TITLE("UnitCell","read_atom_species");
    const int ntype = ucell.ntype;
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
            std::string one_line;
            std::string one_string;
            std::getline(ifa, one_line);
            std::stringstream ss;
            ss << one_line;
            ss >> ucell.atom_label[i] >> ucell.atom_mass[i];
            ucell.pseudo_fn[i] = "auto";
            ucell.pseudo_type[i] = "auto";

            if(!PARAM.inp.use_paw)
            {
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
                    if (one_string == "auto" || one_string == "upf" || one_string == "vwr" || one_string == "upf201" || one_string == "blps")
                    {
                        ucell.pseudo_type[i] = one_string;
                    }
                    else if (one_string == "1/r")
                    {
                        ucell.atoms[i].coulomb_potential = true;
                    }
                    else
                    {
                        GlobalV::ofs_warning << "unrecongnized pseudopotential type: " << one_string << ", check your STRU file." << std::endl;
                        ModuleBase::WARNING_QUIT("read_atom_species", "unrecongnized pseudo type.");
                    }
                }

                // Peize Lin test for bsse 2021.04.07
                const std::string bsse_label = "empty";
                ucell.atoms[i].flag_empty_element = 
                    (search( ucell.atom_label[i].begin(), ucell.atom_label[i].end(), bsse_label.begin(), bsse_label.end() ) != ucell.atom_label[i].end())
                    ? true : false;
            }
        }
    }

    if((PARAM.inp.basis_type == "lcao")
      ||(PARAM.inp.basis_type == "lcao_in_pw")
      ||((PARAM.inp.basis_type == "pw")&&(PARAM.inp.init_wfc.substr(0, 3) == "nao"))
      || PARAM.inp.onsite_radius > 0.0)
    {
        if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "NUMERICAL_ORBITAL") )
        {
            for(int i=0; i<ntype; i++)
            {
                ifa >> ucell.orbital_fn[i];
            }
        }    
        // caoyu add 2021-03-16
        if(PARAM.globalv.deepks_setorb)
        {
            if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "NUMERICAL_DESCRIPTOR")) {
                ifa >> ucell.descriptor_file;
            }
        }
        else
        {
            ucell.descriptor_file = PARAM.inp.orbital_dir + ucell.orbital_fn[0];
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
    if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_CONSTANT") )
    {
        ModuleBase::GlobalFunc::READ_VALUE(ifa, lat0);
        if(lat0<=0.0)
        {
            ModuleBase::WARNING_QUIT("read_atom_species","lat0<=0.0");
        }
        lat0_angstrom = lat0 * 0.529177;
        ModuleBase::GlobalFunc::OUT(ofs_running,"lattice constant (Bohr)",lat0);
        ModuleBase::GlobalFunc::OUT(ofs_running,"lattice constant (Angstrom)",lat0_angstrom);
        lat.tpiba  = ModuleBase::TWO_PI / lat0;
        lat.tpiba2 = lat.tpiba * lat.tpiba;
    }

    //===========================
    // Read in latticies vector
    //===========================

    if(latName=="none")
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa,
                                               "LATTICE_PARAMETERS",
                                               true,
                                               false)) 
        {
            ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant","do not use LATTICE_PARAMETERS without explicit specification of lattice type");
        }
        if( !ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_VECTORS") )
        {
            ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant","Please set LATTICE_VECTORS in STRU file");
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
    else
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa,
                                               "LATTICE_VECTORS",
                                               true,
                                               false)) 
        {
            ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant","do not use LATTICE_VECTORS along with explicit specification of lattice type");
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
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
            }
        }
        else if(latName=="trigonal")
        {//trigonal, ibrav = 5
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
            latvec.e21 = 0.0; latvec.e22 = e22; latvec.e23 = e13;
            latvec.e31 =-e11; latvec.e32 = e12; latvec.e33 = e13;
        }
        else if(latName=="st")
        {//simple tetragonal, ibrav= 6
            latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
            }
        }
        else if(latName=="bct")
        {//body-centered tetragonal, ibrav = 7
            double cba = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cba);
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
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> latvec.e22;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
            }
        }
        else if(latName=="baco")
        {//base-centered orthorhombic, ibrav = 9
            latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 = 0.0;
            latvec.e21 =-0.5; latvec.e22 = 0.0; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = 0.0; latvec.e33 = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> latvec.e12;
                latvec.e12 = latvec.e12 / 2.0;
                latvec.e22 = latvec.e12;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
            }
        }
        else if(latName=="fco")
        {//face-centered orthorhombic, ibrav = 10
            double bba = 0.0; double cba = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cba);
                bba = bba / 2.0; cba = cba / 2.0;
            }
            latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 = cba;
            latvec.e21 = 0.5; latvec.e22 = bba; latvec.e23 = 0.0;
            latvec.e31 = 0.0; latvec.e32 = bba; latvec.e33 = cba;
        }
        else if(latName=="bco")
        {//body-centered orthorhombic, ibrav = 11
            double bba = 0.0; double cba = 0.0;
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cba);
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
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cosab);
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
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cosab);
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
            if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
            {
                ifa >> bba >> cba >> cosab >> cosac;
                ModuleBase::GlobalFunc::READ_VALUE(ifa, cosbc);
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
            std::cout << "latname is : " << latName << std::endl;
            ModuleBase::WARNING_QUIT("unitcell::read_lattice_constant","latname not supported!");
        }
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