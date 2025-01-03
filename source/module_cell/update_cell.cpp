#include "update_cell.h"
#include "bcast_cell.h"
#include "module_base/global_function.h"
namespace unitcell
{
void remake_cell(Lattice& lat)
{
    ModuleBase::TITLE("Lattice", "rmake_cell");

    // The idea is as follows: for each type of lattice, first calculate
    // from current latvec the lattice parameters, then use the parameters
    // to reconstruct latvec
    std::string& latName = lat.latName;
    ModuleBase::Matrix3&  latvec = lat.latvec;

    if (latName == "none") {
        ModuleBase::WARNING_QUIT("UnitCell",
            "to use fixed_ibrav, latname must be provided");
    } else if (latName == "sc") // ibrav = 1
    {
        double celldm = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                  + pow(latvec.e13, 2));

        latvec.Zero();
        latvec.e11 = latvec.e22 = latvec.e33 = celldm;
    } else if (latName == "fcc") // ibrav = 2
    {
        double celldm = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                + pow(latvec.e13, 2)) / std::sqrt(2.0);

        latvec.e11 = -celldm;
        latvec.e12 = 0.0;
        latvec.e13 = celldm;
        latvec.e21 = 0.0;
        latvec.e22 = celldm;
        latvec.e23 = celldm;
        latvec.e31 = -celldm;
        latvec.e32 = celldm;
        latvec.e33 = 0.0;
    } else if (latName == "bcc") // ibrav = 3
    {
        double celldm = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                  + pow(latvec.e13, 2))
                        / std::sqrt(3.0);

        latvec.e11 = celldm;
        latvec.e12 = celldm;
        latvec.e13 = celldm;
        latvec.e21 = -celldm;
        latvec.e22 = celldm;
        latvec.e23 = celldm;
        latvec.e31 = -celldm;
        latvec.e32 = -celldm;
        latvec.e33 = celldm;
    } else if (latName == "hexagonal") // ibrav = 4
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));
        double e22 = sqrt(3.0) / 2.0;

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = -0.5 * celldm1;
        latvec.e22 = celldm1 * e22;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "trigonal") // ibrav = 5
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22
                           + latvec.e13 * latvec.e23);
        double cos12 = celldm12 / celldm1 / celldm2;

        if (cos12 <= -0.5 || cos12 >= 1.0) {
            ModuleBase::WARNING_QUIT("unitcell", "wrong cos12!");
        }
        double t1 = sqrt(1.0 + 2.0 * cos12);
        double t2 = sqrt(1.0 - cos12);

        double e11 = celldm1 * t2 / sqrt(2.0);
        double e12 = -celldm1 * t2 / sqrt(6.0);
        double e13 = celldm1 * t1 / sqrt(3.0);
        double e22 = celldm1 * sqrt(2.0) * t2 / sqrt(3.0);

        latvec.e11 = e11;
        latvec.e12 = e12;
        latvec.e13 = e13;
        latvec.e21 = 0.0;
        latvec.e22 = e22;
        latvec.e23 = e13;
        latvec.e31 = -e11;
        latvec.e32 = e12;
        latvec.e33 = e13;
    } else if (latName == "st") // ibrav = 6
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));
        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = 0.0;
        latvec.e22 = celldm1;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "bct") // ibrav = 7
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::abs(latvec.e13);

        latvec.e11 = celldm1;
        latvec.e12 = -celldm1;
        latvec.e13 = celldm2;
        latvec.e21 = celldm1;
        latvec.e22 = celldm1;
        latvec.e23 = celldm2;
        latvec.e31 = -celldm1;
        latvec.e32 = -celldm1;
        latvec.e33 = celldm2;
    } else if (latName == "so") // ibrav = 8
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = 0.0;
        latvec.e22 = celldm2;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "baco") // ibrav = 9
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::abs(latvec.e22);
        double celldm3 = std::abs(latvec.e33);

        latvec.e11 = celldm1;
        latvec.e12 = celldm2;
        latvec.e13 = 0.0;
        latvec.e21 = -celldm1;
        latvec.e22 = celldm2;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "fco") // ibrav = 10
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::abs(latvec.e22);
        double celldm3 = std::abs(latvec.e33);

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = celldm3;
        latvec.e21 = celldm1;
        latvec.e22 = celldm2;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = celldm2;
        latvec.e33 = celldm3;
    } else if (latName == "bco") // ibrav = 11
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::abs(latvec.e12);
        double celldm3 = std::abs(latvec.e13);

        latvec.e11 = celldm1;
        latvec.e12 = celldm2;
        latvec.e13 = celldm3;
        latvec.e21 = -celldm1;
        latvec.e22 = celldm2;
        latvec.e23 = celldm3;
        latvec.e31 = -celldm1;
        latvec.e32 = -celldm2;
        latvec.e33 = celldm3;
    } else if (latName == "sm") // ibrav = 12
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));
        double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22
                           + latvec.e13 * latvec.e23);
        double cos12 = celldm12 / celldm1 / celldm2;

        double e21 = celldm2 * cos12;
        double e22 = celldm2 * std::sqrt(1.0 - cos12 * cos12);

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = e21;
        latvec.e22 = e22;
        latvec.e23 = 0.0;
        latvec.e31 = 0.0;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "bacm") // ibrav = 13
    {
        double celldm1 = std::abs(latvec.e11);
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm3 = std::abs(latvec.e13);

        double cos12 = latvec.e21 / celldm2;
        if (cos12 >= 1.0) {
            ModuleBase::WARNING_QUIT("unitcell", "wrong cos12!");
        }

        double e21 = celldm2 * cos12;
        double e22 = celldm2 * std::sqrt(1.0 - cos12 * cos12);

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = -celldm3;
        latvec.e21 = e21;
        latvec.e22 = e22;
        latvec.e23 = 0.0;
        latvec.e31 = celldm1;
        latvec.e32 = 0.0;
        latvec.e33 = celldm3;
    } else if (latName == "triclinic") // ibrav = 14
    {
        double celldm1 = std::sqrt(pow(latvec.e11, 2) + pow(latvec.e12, 2)
                                   + pow(latvec.e13, 2));
        double celldm2 = std::sqrt(pow(latvec.e21, 2) + pow(latvec.e22, 2)
                                   + pow(latvec.e23, 2));
        double celldm3 = std::sqrt(pow(latvec.e31, 2) + pow(latvec.e32, 2)
                                   + pow(latvec.e33, 2));
        double celldm12 = (latvec.e11 * latvec.e21 + latvec.e12 * latvec.e22
                           + latvec.e13 * latvec.e23);
        double cos12 = celldm12 / celldm1 / celldm2;
        double celldm13 = (latvec.e11 * latvec.e31 + latvec.e12 * latvec.e32
                           + latvec.e13 * latvec.e33);
        double cos13 = celldm13 / celldm1 / celldm3;
        double celldm23 = (latvec.e21 * latvec.e31 + latvec.e22 * latvec.e32
                           + latvec.e23 * latvec.e33);
        double cos23 = celldm23 / celldm2 / celldm3;

        double sin12 = std::sqrt(1.0 - cos12 * cos12);
        if (cos12 >= 1.0) {
            ModuleBase::WARNING_QUIT("unitcell", "wrong cos12!");
        }

        latvec.e11 = celldm1;
        latvec.e12 = 0.0;
        latvec.e13 = 0.0;
        latvec.e21 = celldm2 * cos12;
        latvec.e22 = celldm2 * sin12;
        latvec.e23 = 0.0;
        latvec.e31 = celldm3 * cos13;
        latvec.e32 = celldm3 * (cos23 - cos13 * cos12) / sin12;
        double term = 1.0 + 2.0 * cos12 * cos13 * cos23 - cos12 * cos12
                      - cos13 * cos13 - cos23 * cos23;
        term = sqrt(term) / sin12;
        latvec.e33 = celldm3 * term;
    } else {
        std::cout << "latname is : " << latName << std::endl;
        ModuleBase::WARNING_QUIT("UnitCell::read_atom_species",
                                 "latname not supported!");
    }
}

// LiuXh add a new function here,
// 20180515
void setup_cell_after_vc(UnitCell& ucell, std::ofstream& log) {
    ModuleBase::TITLE("UnitCell", "setup_cell_after_vc");
    assert(ucell.lat0 > 0.0);
    ucell.omega = std::abs(ucell.latvec.Det()) * 
                           pow(ucell.lat0, 3);
    if (ucell.omega <= 0) 
    {
        ModuleBase::WARNING_QUIT("setup_cell_after_vc", "omega <= 0 .");
    } else {
        log << std::endl;
        ModuleBase::GlobalFunc::OUT(log, "Volume (Bohr^3)", ucell.omega);
        ModuleBase::GlobalFunc::OUT(log, "Volume (A^3)",
                                    ucell.omega * pow(ModuleBase::BOHR_TO_A, 3));
    }

    ucell.lat0_angstrom = ucell.lat0 * 0.529177;
    ucell.tpiba = ModuleBase::TWO_PI / ucell.lat0;
    ucell.tpiba2 = ucell.tpiba * ucell.tpiba;

    // lattice vectors in another form.
    ucell.a1.x = ucell.latvec.e11;
    ucell.a1.y = ucell.latvec.e12;
    ucell.a1.z = ucell.latvec.e13;

    ucell.a2.x = ucell.latvec.e21;
    ucell.a2.y = ucell.latvec.e22;
    ucell.a2.z = ucell.latvec.e23;

    ucell.a3.x = ucell.latvec.e31;
    ucell.a3.y = ucell.latvec.e32;
    ucell.a3.z = ucell.latvec.e33;

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec has the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    ucell.GT = ucell.latvec.Inverse();
    ucell.G = ucell.GT.Transpose();
    ucell.GGT = ucell.G * ucell.GT;
    ucell.invGGT = ucell.GGT.Inverse();

    for (int it = 0; it < ucell.ntype; it++) {
        Atom* atom = &ucell.atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            atom->tau[ia] = atom->taud[ia] * ucell.latvec;
        }
    }

#ifdef __MPI
    ucell.bcast_unitcell();
#endif

    log << std::endl;
    output::printM3(log,
                    "Lattice vectors: (Cartesian coordinate: in unit of a_0)",
                    ucell.latvec);
    output::printM3(log,
                   "Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",
                    ucell.G);

    return;
}

void update_pos_tau(const Lattice& lat,
                    const double* pos,
                    const int ntype,
                    const int nat,
                    Atom* atoms) 
{
    int iat = 0;
    for (int it = 0; it < ntype; it++) {
        Atom* atom = &atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            for (int ik = 0; ik < 3; ++ik) {
                if (atom->mbl[ia][ik]) 
                {
                    atom->dis[ia][ik] = pos[3 * iat + ik] / lat.lat0 - atom->tau[ia][ik];
                    atom->tau[ia][ik] = pos[3 * iat + ik] / lat.lat0;
                }
            }
            // the direct coordinates also need to be updated.
            atom->dis[ia] = atom->dis[ia] * lat.GT;
            atom->taud[ia] = atom->tau[ia] * lat.GT;
            iat++;
        }
    }
    assert(iat == nat);
    periodic_boundary_adjustment(atoms,lat.latvec,ntype);
    bcast_atoms_tau(atoms, ntype);
}

void periodic_boundary_adjustment(Atom* atoms,
                                  const ModuleBase::Matrix3& latvec,
                                  const int ntype) 
{
    //----------------------------------------------
    // because of the periodic boundary condition
    // we need to adjust the atom positions,
    // first adjust direct coordinates,
    // then update them into cartesian coordinates,
    //----------------------------------------------
    for (int it = 0; it < ntype; it++) {
        Atom* atom = &atoms[it];
        for (int ia = 0; ia < atom->na; ia++) {
            // mohan update 2011-03-21
            for (int ik = 0; ik < 3; ik++) 
            {
                if (atom->taud[ia][ik] < 0) 
                {
                    atom->taud[ia][ik] += 1.0;
                }
                if (atom->taud[ia][ik] >= 1.0) 
                {
                    atom->taud[ia][ik] -= 1.0;
                }
            }
            if (atom->taud[ia].x < 0 
                || atom->taud[ia].y < 0
                || atom->taud[ia].z < 0 
                || atom->taud[ia].x >= 1.0
                || atom->taud[ia].y >= 1.0 
                || atom->taud[ia].z >= 1.0) 
            {
                GlobalV::ofs_warning << " it=" << it + 1 << " ia=" << ia + 1 << std::endl;
                GlobalV::ofs_warning << "d=" << atom->taud[ia].x << " "
                                     << atom->taud[ia].y << " "
                                     << atom->taud[ia].z << std::endl;
                ModuleBase::WARNING_QUIT("Ions_Move_Basic::move_ions",
                    "the movement of atom is larger than the length of cell.");
            }

            atom->tau[ia] = atom->taud[ia] * latvec;
        }
    }
    return;
}

} // namespace unitcell