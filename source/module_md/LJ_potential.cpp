#include "LJ_potential.h"
#include "../input.h"
#include "../module_base/timer.h"

LJ_potential::LJ_potential(){}

LJ_potential::~LJ_potential(){}

double LJ_potential::Lennard_Jones(const UnitCell_pseudo &ucell_c, 
                    Grid_Driver &grid_neigh, 
                    ModuleBase::Vector3<double> *force,
                    ModuleBase::matrix &stress)
{
    ModuleBase::TITLE("LJ_potential", "Lennard_Jones");
    ModuleBase::timer::tick("LJ_potential", "Lennard_Jones");

    double distance, potential = 0; //initialize
    int index = 0;
    double virial[6];
    ModuleBase::GlobalFunc::ZEROS(virial, 6);    // initialize

    ModuleBase::Vector3<double> tau1, tau2, dtau;
    for(int it=0; it<ucell_c.ntype; ++it)
    {
        Atom* atom1 = &ucell_c.atoms[it];
        for(int ia=0; ia<atom1->na; ++ia)
	    {	
		    force[index].set(0,0,0);   //initialize
            tau1 = atom1->tau[ia];
            grid_neigh.Find_atom(ucell_c, tau1, it, ia);
            for(int ad=0; ad<grid_neigh.getAdjacentNum(); ++ad)
            {
                tau2 = grid_neigh.getAdjacentTau(ad);
                dtau = (tau1 - tau2) * ucell_c.lat0;
                distance = dtau.norm();
                if(distance <= INPUT.mdp.lj_rcut)
			    {
				    potential += LJ_energy(distance); // - LJ_energy(INPUT.mdp.lj_rcut);
                    ModuleBase::Vector3<double> f_ij = LJ_force(distance, dtau);
				    force[index] += f_ij;
                    LJ_virial(virial, f_ij, dtau);
			    }
            }
            index++;
	    }
    }

    // Post treatment for virial
    stress(0, 0) = virial[0]/(2.0*ucell_c.omega);
    stress(1, 1) = virial[1]/(2.0*ucell_c.omega);
    stress(2, 2) = virial[2]/(2.0*ucell_c.omega);
    stress(0, 1) = stress(1, 0) = virial[3]/(2.0*ucell_c.omega);
    stress(0, 2) = stress(2, 0) = virial[4]/(2.0*ucell_c.omega);
    stress(1, 2) = stress(2, 1) = virial[5]/(2.0*ucell_c.omega);

    ModuleBase::timer::tick("LJ_potential", "Lennard_Jones");
	return potential/2.0;
}

#include "../module_base/mathzone.h"
double LJ_potential::Lennard_Jones(const UnitCell_pseudo &ucell_c, 
                    CMD_neighbor &cmd_neigh,
                    ModuleBase::Vector3<double> *force,
                    ModuleBase::matrix &stress)
{
    ModuleBase::TITLE("LJ_potential", "Lennard_Jones");
    ModuleBase::timer::tick("LJ_potential", "Lennard_Jones");

    double potential = 0; // initialize
    double virial[6];
    ModuleBase::GlobalFunc::ZEROS(virial, 6);    // initialize

	for(int i=0; i<ucell_c.nat; i++)
	{	
		force[i].set(0,0,0);   //initialize

        int T1 = ucell_c.iat2it[i];
        int I1 = ucell_c.iat2ia[i];
        ModuleBase::Vector3<double> taud1 = ucell_c.atoms[T1].taud[I1];

		for(int j=0; j<cmd_neigh.nlist[i]; j++)
		{
            int T2 = ucell_c.iat2it[cmd_neigh.list[i][j]];
            int I2 = ucell_c.iat2ia[cmd_neigh.list[i][j]];
            ModuleBase::Vector3<double> taud2 = ucell_c.atoms[T2].taud[I2];

			ModuleBase::Vector3<double> tempd = cmd_neigh.cell_periodic(taud1,taud2);
            ModuleBase::Vector3<double> temp;
            ModuleBase::Mathzone::Direct_to_Cartesian(
			tempd.x, tempd.y, tempd.z,
			ucell_c.latvec.e11, ucell_c.latvec.e12, ucell_c.latvec.e13,
			ucell_c.latvec.e21, ucell_c.latvec.e22, ucell_c.latvec.e23,
			ucell_c.latvec.e31, ucell_c.latvec.e32, ucell_c.latvec.e33,
			temp.x, temp.y, temp.z);

            double distance = temp.norm()*ucell_c.lat0;

			if(distance <= INPUT.mdp.lj_rcut)
			{
				potential += LJ_energy(distance); // - LJ_energy(INPUT.mdp.lj_rcut);
                ModuleBase::Vector3<double> f_ij = LJ_force(distance, temp*ucell_c.lat0);
				force[i] = force[i] + f_ij;
                LJ_virial(virial, f_ij, temp*ucell_c.lat0);
			}
		}
	}

    // Post treatment for virial
    stress(0, 0) = virial[0]/(2.0*ucell_c.omega);
    stress(1, 1) = virial[1]/(2.0*ucell_c.omega);
    stress(2, 2) = virial[2]/(2.0*ucell_c.omega);
    stress(0, 1) = stress(1, 0) = virial[3]/(2.0*ucell_c.omega);
    stress(0, 2) = stress(2, 0) = virial[4]/(2.0*ucell_c.omega);
    stress(1, 2) = stress(2, 1) = virial[5]/(2.0*ucell_c.omega);

    ModuleBase::timer::tick("LJ_potential", "Lennard_Jones");
	return potential/2.0;
}

double LJ_potential::LJ_energy(const double d)
{
	return 4*INPUT.mdp.lj_epsilon*( pow(INPUT.mdp.lj_sigma/d, 12) - pow(INPUT.mdp.lj_sigma/d, 6) );
}

ModuleBase::Vector3<double> LJ_potential::LJ_force(const double d, const ModuleBase::Vector3<double> dr)
{
	double coff = 4*INPUT.mdp.lj_epsilon*( 12*pow(INPUT.mdp.lj_sigma/d, 12) - 6*pow(INPUT.mdp.lj_sigma/d, 6) )/pow(d,2);
	return dr*coff;
}

void LJ_potential::LJ_virial(double *virial, const ModuleBase::Vector3<double> &force, const ModuleBase::Vector3<double> &dtau)
{
    virial[0] += dtau.x * force.x;
    virial[1] += dtau.y * force.y;
    virial[2] += dtau.z * force.z;
    virial[3] += dtau.x * force.y;
    virial[4] += dtau.x * force.z;
    virial[5] += dtau.y * force.z;
}