#include <unistd.h>
#include "esolver_dp.h"
#ifdef __DPMD
#include "deepmd/DeepPot.h"
#endif


namespace ModuleESolver
{

void ESolver_DP::Init(Input &inp, UnitCell_pseudo &ucell)
{
    dp_potential = 0;
    dp_force.create(ucell.nat,3);
    dp_virial.create(3,3);

    cell.resize(9);
    atype.resize(ucell.nat);
    coord.resize(3*ucell.nat);

    cell[0] = ucell.latvec.e11*ucell.lat0_angstrom;
    cell[1] = ucell.latvec.e12*ucell.lat0_angstrom;
    cell[2] = ucell.latvec.e13*ucell.lat0_angstrom;
    cell[3] = ucell.latvec.e21*ucell.lat0_angstrom;
    cell[4] = ucell.latvec.e22*ucell.lat0_angstrom;
    cell[5] = ucell.latvec.e23*ucell.lat0_angstrom;
    cell[6] = ucell.latvec.e31*ucell.lat0_angstrom;
    cell[7] = ucell.latvec.e32*ucell.lat0_angstrom;
    cell[8] = ucell.latvec.e33*ucell.lat0_angstrom;

    int iat = 0;
    for(int it=0; it<ucell.ntype; ++it)
    {
        for(int ia=0; ia<ucell.atoms[it].na; ++ia)
        {
            atype[iat] = it;
            coord[3*iat] = ucell.atoms[it].tau[ia].x * ucell.lat0_angstrom;
            coord[3*iat+1] = ucell.atoms[it].tau[ia].y * ucell.lat0_angstrom;
            coord[3*iat+2] = ucell.atoms[it].tau[ia].z * ucell.lat0_angstrom;
            iat++;
        }
    }
    assert(ucell.nat == iat);
}

void ESolver_DP::Run(int istep, UnitCell_pseudo &ucell)
{
#ifdef __DPMD
    if(access("graph.pb", 0) == -1)
    {
        ModuleBase::WARNING_QUIT("DP_pot", "Can not find graph.pb !");
    }

    deepmd::DeepPot dp ("graph.pb");

    std::vector<double> f, v;

    dp.compute (dp_potential, f, v, coord, atype, cell);

    dp_potential /= ModuleBase::Hartree_to_eV;

    double fact_f = ModuleBase::Hartree_to_eV*ModuleBase::ANGSTROM_AU;
    double fact_v = ModuleBase::Hartree_to_eV*pow(ModuleBase::ANGSTROM_AU, 3);

    for(int i=0; i<ucell.nat;  ++i)
    {
        dp_force(i,0) = f[3*i]/fact_f;
        dp_force(i,1) = f[3*i+1]/fact_f;
        dp_force(i,2) = f[3*i+2]/fact_f;
    }

    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            dp_virial(i, j) = v[3*i+j]/fact_v;
        }
    }
#else
    ModuleBase::WARNING_QUIT("DP_pot", "Please recompile with -D__DPMD !");
#endif
}

void ESolver_DP::cal_Energy(energy &en)
{

}

void ESolver_DP::cal_Force(ModuleBase::matrix &force)
{
    force = dp_force;
}

void ESolver_DP::cal_Stress(ModuleBase::matrix &stress)
{
    stress = dp_virial;
}

}