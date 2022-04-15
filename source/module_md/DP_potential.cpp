#include <unistd.h>
#include "DP_potential.h"
#ifdef __DPMD
#include "deepmd/DeepPot.h"
#endif

DP_potential::DP_potential(){}

DP_potential::~DP_potential(){}

void DP_potential::DP_pot(const UnitCell_pseudo &ucell_c, double &potential, ModuleBase::Vector3<double> *force, ModuleBase::matrix &stress)
{
#ifdef __DPMD
    ModuleBase::TITLE("DP_potential", "DP_pot");
    ModuleBase::timer::tick("DP_potential", "DP_pot");

    if(access("graph.pb", 0) == -1)
    {
        ModuleBase::WARNING_QUIT("DP_pot", "Can not find greph.pb !");
    }

    deepmd::DeepPot dp ("graph.pb");

    std::vector<double> cell(9);
    std::vector<int> atype;
    std::vector<double> coord;

    DP_init(ucell_c, cell, coord, atype);

    std::vector<double> f, v;

    dp.compute (potential, f, v, coord, atype, cell);

    potential/=ModuleBase::Hartree_to_eV;

    double fact_f = ModuleBase::Hartree_to_eV*ModuleBase::ANGSTROM_AU;
    double fact_v = ModuleBase::Hartree_to_eV*pow(ModuleBase::ANGSTROM_AU, 3);

    for(int i=0; i<ucell_c.nat;  ++i)
    {
        force[i].x = f[3*i]/fact_f;
        force[i].y = f[3*i+1]/fact_f;
        force[i].z = f[3*i+2]/fact_f;
    }

    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            stress(i, j) = v[3*i+j]/fact_v;
        }
    }

    ModuleBase::timer::tick("DP_potential", "DP_pot");
#else
    ModuleBase::WARNING_QUIT("DP_pot", "Please recompile with -D__DPMD !");
#endif
}

void DP_potential::DP_init(const UnitCell_pseudo &ucell_c, 
                std::vector<double> &cell, 
                std::vector<double> &coord, 
                std::vector<int> &atype)
{
    atype.resize(ucell_c.nat);
    coord.resize(3*ucell_c.nat);

    cell[0] = ucell_c.latvec.e11*ucell_c.lat0_angstrom;
    cell[1] = ucell_c.latvec.e12*ucell_c.lat0_angstrom;
    cell[2] = ucell_c.latvec.e13*ucell_c.lat0_angstrom;
    cell[3] = ucell_c.latvec.e21*ucell_c.lat0_angstrom;
    cell[4] = ucell_c.latvec.e22*ucell_c.lat0_angstrom;
    cell[5] = ucell_c.latvec.e23*ucell_c.lat0_angstrom;
    cell[6] = ucell_c.latvec.e31*ucell_c.lat0_angstrom;
    cell[7] = ucell_c.latvec.e32*ucell_c.lat0_angstrom;
    cell[8] = ucell_c.latvec.e33*ucell_c.lat0_angstrom;

    int iat = 0;
    for(int it=0; it<ucell_c.ntype; ++it)
    {
        for(int ia=0; ia<ucell_c.atoms[it].na; ++ia)
        {
            atype[iat] = it;
            coord[3*iat] = ucell_c.atoms[it].tau[ia].x * ucell_c.lat0_angstrom;
            coord[3*iat+1] = ucell_c.atoms[it].tau[ia].y * ucell_c.lat0_angstrom;
            coord[3*iat+2] = ucell_c.atoms[it].tau[ia].z * ucell_c.lat0_angstrom;
            iat++;
        }
    }
    assert(ucell_c.nat == iat);
}