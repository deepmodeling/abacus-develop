#include "symmetry_test.h"

void SymmetryTest::construct_ucell(stru_ &stru)
{
    std::vector<atomtype_> coord = stru.all_type;
    ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
    ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
    ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
    ucell.latvec.e11=ucell.a1.x;
    ucell.latvec.e12=ucell.a1.y;
    ucell.latvec.e13=ucell.a1.z;
    ucell.latvec.e21=ucell.a2.x;
    ucell.latvec.e22=ucell.a2.y;
    ucell.latvec.e23=ucell.a2.z;
    ucell.latvec.e31=ucell.a3.x;
    ucell.latvec.e32=ucell.a3.y;
    ucell.latvec.e33=ucell.a3.z;
    ucell.GT = ucell.latvec.Inverse();
    ucell.G = ucell.GT.Transpose();
    ucell.ntype = stru.all_type.size();
    ucell.atoms = new Atom[ucell.ntype];
    ucell.nat = 0;

    for (int i = 0; i < coord.size(); i++)
    {
        ucell.atoms[i].label = coord[i].atomname;
        ucell.atoms[i].na = coord[i].coordinate.size();
        ucell.atoms[i].tau = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
        ucell.atoms[i].taud = new ModuleBase::Vector3<double>[ucell.atoms[i].na];
        for (int j = 0; j < ucell.atoms[i].na; j++)
        {
            std::vector<double> this_atom = coord[i].coordinate[j];
            if (stru.coordtype == "C")
            {
                ucell.atoms[i].tau[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
                ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[i].tau[j].x,
                                                        ucell.atoms[i].tau[j].y,
                                                        ucell.atoms[i].tau[j].z,
                                                        ucell.a1.x,
                                                        ucell.a1.y,
                                                        ucell.a1.z,
                                                        ucell.a2.x,
                                                        ucell.a2.y,
                                                        ucell.a2.z,
                                                        ucell.a3.x,
                                                        ucell.a3.y,
                                                        ucell.a3.z,
                                                        ucell.atoms[i].taud[j].x,
                                                        ucell.atoms[i].taud[j].y,
                                                        ucell.atoms[i].taud[j].z);
            }
            else
            {
                ucell.atoms[i].taud[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
                ModuleBase::Mathzone::Direct_to_Cartesian(ucell.atoms[i].taud[j].x,
                                                        ucell.atoms[i].taud[j].y,
                                                        ucell.atoms[i].taud[j].z,
                                                        ucell.a1.x,
                                                        ucell.a1.y,
                                                        ucell.a1.z,
                                                        ucell.a2.x,
                                                        ucell.a2.y,
                                                        ucell.a2.z,
                                                        ucell.a3.x,
                                                        ucell.a3.y,
                                                        ucell.a3.z,
                                                        ucell.atoms[i].tau[j].x,
                                                        ucell.atoms[i].tau[j].y,
                                                        ucell.atoms[i].tau[j].z);
            }

        }
        ucell.nat += ucell.atoms[i].na;
    }
}

void SymmetryTest::ClearUcell()
{
    for (int i = 0; i < ucell.ntype; i++)
    {
        delete[] ucell.atoms[i].tau;
        delete[] ucell.atoms[i].taud;
    }
    delete[] ucell.atoms;
}