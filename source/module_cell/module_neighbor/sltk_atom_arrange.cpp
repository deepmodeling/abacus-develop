#include "sltk_atom_arrange.h"

#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "sltk_grid.h"
#include "sltk_grid_driver.h"

// update the followig class in near future
#include "module_cell/unitcell.h"

atom_arrange::atom_arrange()
{
}

atom_arrange::~atom_arrange()
{
}

double atom_arrange::set_sr_NL(std::ofstream& ofs_in,
                               const std::string& output_level,
                               const double& rcutmax_Phi,
                               const double& rcutmax_Beta,
                               const bool gamma_only_local)
{
    ModuleBase::TITLE("atom_arrange", "set_sr_NL");
    // check in use_overlap_matrix,
    double sr = 0.0;
    if (gamma_only_local)
    {
        sr = 2 * rcutmax_Phi + 0.001;
    }
    else
    {
        sr = 2 * (rcutmax_Phi + rcutmax_Beta) + 0.001; // 0.001 is added to make safe.
                                                      // sr = 2 * longest_orb_rcut + 0.001;
    }

    if (output_level != "m") // xiaohui add 'output_level', 2015-09-16
    {
        ofs_in << "\n\n\n\n";
        ofs_in << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        ofs_in << " |                                                                    |" << std::endl;
        ofs_in << " | Search adjacent atoms:                                             |" << std::endl;
        ofs_in << " | Set the adjacent atoms for each atom and set the periodic boundary |" << std::endl;
        ofs_in << " | condition for the atoms on real space FFT grid. For k-dependent    |" << std::endl;
        ofs_in << " | algorithm, we also need to set the sparse H and S matrix element   |" << std::endl;
        ofs_in << " | for each atom.                                                     |" << std::endl;
        ofs_in << " |                                                                    |" << std::endl;
        ofs_in << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        ofs_in << "\n\n\n\n";

        ofs_in << "\n SETUP SEARCHING RADIUS FOR PROGRAM TO SEARCH ADJACENT ATOMS" << std::endl;
        ofs_in << std::setprecision(3);
        ModuleBase::GlobalFunc::OUT(ofs_in, "longest orb rcut (Bohr)", rcutmax_Phi);
        ModuleBase::GlobalFunc::OUT(ofs_in, "longest nonlocal projector rcut (Bohr)", rcutmax_Beta);
        ModuleBase::GlobalFunc::OUT(ofs_in, "search radius (Bohr)", sr);
	}
    return sr;
}

void atom_arrange::search(const bool pbc_flag,
                          std::ofstream& ofs_in,
                          Grid_Driver& grid_d,
                          const UnitCell& ucell,
                          const double& search_radius_bohr,
                          const int& test_atom_in,
                          const bool test_only)
{
    ModuleBase::TITLE("atom_arrange", "search");
    ModuleBase::timer::tick("atom_arrange", "search");

    if (search_radius_bohr < 0.0)
    {
        ModuleBase::WARNING_QUIT("atom_arrange::search", " search_radius_bohr < 0,forbidden");
    }

    ModuleBase::GlobalFunc::OUT(ofs_in, "searching radius is (Bohr))", search_radius_bohr);
    ModuleBase::GlobalFunc::OUT(ofs_in, "searching radius unit is (Bohr))", ucell.lat0);

    assert(ucell.nat > 0);

    /*
    2024-12-04 Zhang Haochong
        The neighboring atom search module has been completely rewritten.
        The new algorithm places atoms into boxes with an edge length of twice the atomic radius. The neighboring 
    atom list stores the data using the atom's type and its index within that type.
        By setting pbc_flag = false, periodic boundary conditions can be forcibly disabled. In this case, the search 
    process will not expand the supercell, and the neighboring atoms will only consider those within the original unit cell.
    */
    const double radius_lat0unit = search_radius_bohr / ucell.lat0;

    // Atom_input at(ofs_in, ucell, pbc_flag, radius_lat0unit, test_atom_in);

    grid_d.init(ofs_in, ucell, radius_lat0unit, pbc_flag);

	// The screen output is very time-consuming. To avoid interfering with the timing, we will insert logging here earlier.
    ModuleBase::timer::tick("atom_arrange", "search");

    if (test_only)
    {
        std::cout << "radius_lat0unit = " << radius_lat0unit << std::endl;
        std::cout << "search_radius_bohr = " << search_radius_bohr << std::endl;

        ofs_in << " " << std::setw(5) << "Type" << std::setw(5) << "Atom" << std::setw(8) << "AdjNum" << std::endl;
        std::cout << std::setw(8) << "Labels" << std::setw(15) << "tau.x" << std::setw(15) << "tau.y" << std::setw(15)
                  << "tau.z" << std::setw(8) << "box.x" << std::setw(8) << "box.y" << std::setw(8) << "box.z"
                  << std::endl;
        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                grid_d.Find_atom(ucell, ucell.atoms[it].tau[ia], it, ia);

                ofs_in << " " << std::setw(5) << it << std::setw(5) << ia << std::setw(8) << grid_d.getAdjacentNum() + 1
                       << std::endl;
                std::cout << " adjacent atoms of " << ucell.atoms[it].label + std::to_string(ia + 1) << ":" << std::endl;
                std::cout << "getAdjacentNum: " << grid_d.getAdjacentNum() + 1 << std::endl;
                /*
                for (int ad = 0; ad < grid_d.getAdjacentNum() + 1; ad++)
                {
                    ModuleBase::Vector3<double> tau = grid_d.getAdjacentTau(ad);
                    ModuleBase::Vector3<int> box = grid_d.getBox(ad);
                    std::cout << std::setw(8) << ucell.atoms[it].label + std::to_string(ia + 1) << std::setw(15)
                              << tau.x << " " << std::setw(15) << tau.y << " " << std::setw(15) << tau.z << " "
                              << std::setw(8) << box.x << std::setw(8) << box.y << std::setw(8) << box.z << std::endl;
                }*/
            }
        }
        ofs_in << "search neighboring atoms done." << std::endl;
    }

    return;
}
