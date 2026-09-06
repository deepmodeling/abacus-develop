#include "dfpt_stru_fixture.h"

#include "source_base/mathzone.h"

DFPTStruTestFixture::DFPTStruTestFixture()
{
    stru_lib.push_back(stru_{1,
                             "O_h",
                             "m-3m",
                             "Pm-3m",
                             std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
                             std::vector<atomtype_>{atomtype_{"C",
                                                              std::vector<std::vector<double>>{
                                                                  {0., 0., 0.},
                                                              }}}});
}

void DFPTStruTestFixture::construct_ucell(stru_& stru)
{
    std::vector<atomtype_> coord = stru.all_type;
    ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
    ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
    ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
    ucell.ntype = stru.all_type.size();
    ucell.atoms = new Atom[ucell.ntype];
    ucell.nat = 0;
    ucell.latvec.e11 = ucell.a1.x;
    ucell.latvec.e12 = ucell.a1.y;
    ucell.latvec.e13 = ucell.a1.z;
    ucell.latvec.e21 = ucell.a2.x;
    ucell.latvec.e22 = ucell.a2.y;
    ucell.latvec.e23 = ucell.a2.z;
    ucell.latvec.e31 = ucell.a3.x;
    ucell.latvec.e32 = ucell.a3.y;
    ucell.latvec.e33 = ucell.a3.z;
    ucell.GT = ucell.latvec.Inverse();
    ucell.G = ucell.GT.Transpose();
    ucell.lat0 = 1.8897261254578281;
    for (size_t i = 0; i < coord.size(); i++)
    {
        ucell.atoms[i].label = coord[i].atomname;
        ucell.atoms[i].na = coord[i].coordinate.size();
        ucell.atoms[i].tau.resize(ucell.atoms[i].na);
        ucell.atoms[i].taud.resize(ucell.atoms[i].na);
        for (int j = 0; j < ucell.atoms[i].na; ++j)
        {
            std::vector<double> this_atom = coord[i].coordinate[j];
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
        ucell.nat += ucell.atoms[i].na;
    }
}

void DFPTStruTestFixture::ClearUcell()
{
    delete[] ucell.atoms;
}
