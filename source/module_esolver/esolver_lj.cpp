#include "esolver_lj.h"

#include "../module_neighbor/sltk_atom_arrange.h"
#include "../module_neighbor/sltk_grid_driver.h"

namespace ModuleESolver
{

    void ESolver_LJ::Init(Input& inp, UnitCell_pseudo& ucell)
    {
        lj_potential = 0;
        lj_force.create(ucell.nat, 3);
        lj_virial.create(3, 3);

        lj_rcut = inp.mdp.lj_rcut;
        lj_epsilon = inp.mdp.lj_epsilon;
        lj_sigma = inp.mdp.lj_sigma;

        atom_arrange::search(
            GlobalV::SEARCH_PBC,
            GlobalV::ofs_running,
            grid_neigh,
            ucell,
            GlobalV::SEARCH_RADIUS,
            GlobalV::test_atom_input);
    }

    void ESolver_LJ::Run(const int istep, UnitCell_pseudo& ucell)
    {
        double distance;
        int index = 0;

        ModuleBase::Vector3<double> tau1, tau2, dtau;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            Atom* atom1 = &ucell.atoms[it];
            for (int ia = 0; ia < atom1->na; ++ia)
            {
                tau1 = atom1->tau[ia];
                grid_neigh.Find_atom(ucell, tau1, it, ia);
                for (int ad = 0; ad < grid_neigh.getAdjacentNum(); ++ad)
                {
                    tau2 = grid_neigh.getAdjacentTau(ad);
                    dtau = (tau1 - tau2) * ucell.lat0;
                    distance = dtau.norm();
                    if (distance <= lj_rcut)
                    {
                        lj_potential += LJ_energy(distance); // - LJ_energy(lj_rcut);
                        ModuleBase::Vector3<double> f_ij = LJ_force(distance, dtau);
                        lj_force(index, 0) += f_ij.x;
                        lj_force(index, 1) += f_ij.y;
                        lj_force(index, 2) += f_ij.z;
                        LJ_virial(f_ij, dtau);
                    }
                }
                index++;
            }
        }

        lj_potential /= 2.0;

        // Post treatment for virial
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                lj_virial(i, j) /= (2.0 * ucell.omega);
            }
        }
    }

    void ESolver_LJ::cal_Energy(energy& en)
    {

    }

    void ESolver_LJ::cal_Force(ModuleBase::matrix& force)
    {
        force = lj_force;
    }

    void ESolver_LJ::cal_Stress(ModuleBase::matrix& stress)
    {
        stress = lj_virial;
    }

    double ESolver_LJ::LJ_energy(const double d)
    {
        double temp = pow(lj_sigma / d, 6);
        return 4 * lj_epsilon * (temp - 1) * temp;
    }

    ModuleBase::Vector3<double> ESolver_LJ::LJ_force(const double d, const ModuleBase::Vector3<double> dr)
    {
        double temp = pow(lj_sigma / d, 6);
        double coff = 24 * lj_epsilon * (2 * temp - 1) * temp / pow(d, 2);
        return dr * coff;
    }

    void ESolver_LJ::LJ_virial(const ModuleBase::Vector3<double>& force,
        const ModuleBase::Vector3<double>& dtau)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                lj_virial(i, j) += dtau[i] * force[j];
            }
        }
    }

}
