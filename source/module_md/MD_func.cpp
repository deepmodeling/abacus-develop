#include "MD_func.h"
#include "../module_base/global_variable.h"
#include "../module_base/timer.h"


double MD_func::gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) 
    {
        do 
        {
            double U1 = rand()/double(RAND_MAX);
            double U2 = rand()/double(RAND_MAX);
             
            V1 = 2.0 * U1 - 1.0;
            V2 = 2.0 * U2 - 1.0;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2.0 * log(S) / S);
    } 
    else
    {
        X = V2 * sqrt(-2.0 * log(S) / S);
    }

    phase = 1 - phase;
 
    return X;
}

double MD_func::GetAtomKE(
		const int &numIon,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass)
{
	double ke = 0;

	for(int ion=0; ion<numIon; ++ion)
	{
		ke += 0.5 * allmass[ion] * vel[ion].norm2();
	}

	return ke;
}

void MD_func::compute_stress(
		const UnitCell &unit_in,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass, 
        const ModuleBase::matrix &virial,
		ModuleBase::matrix &stress)
{
//--------------------------------------------------------------------------------------------
// DESCRIPTION:
//   This function calculates the contribution of classical kinetic energy of atoms to stress.
//--------------------------------------------------------------------------------------------

    if(GlobalV::CAL_STRESS)
    {
        ModuleBase::matrix t_vector;

        temp_vector(unit_in.nat, vel, allmass, t_vector);

        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                stress(i, j) = virial(i, j) + t_vector(i, j) / unit_in.omega;
            }
        }
    }
}

// Read Velocity from STRU liuyu 2021-09-24
void MD_func::ReadVel(
	const UnitCell &unit_in, 
	ModuleBase::Vector3<double>* vel)
{
	int iat=0;
    for(int it=0; it<unit_in.ntype; ++it)
    {
        for(int ia=0; ia<unit_in.atoms[it].na; ++ia)
        {
            vel[iat] = unit_in.atoms[it].vel[ia];
			if(unit_in.atoms[it].mbl[ia].x==0) vel[iat].x = 0;
			if(unit_in.atoms[it].mbl[ia].y==0) vel[iat].y = 0;
			if(unit_in.atoms[it].mbl[ia].z==0) vel[iat].z = 0;
            ++iat;
        }
    }
    assert(iat==unit_in.nat);
}

// Initial velocity randomly
void MD_func::RandomVel(
	const int& numIon, 
	const double& temperature, 
	const double* allmass,
	const int& frozen_freedom,
	const ModuleBase::Vector3<int> frozen,
	const ModuleBase::Vector3<int>* ionmbl,
	ModuleBase::Vector3<double>* vel)
{
    if(!GlobalV::MY_RANK)
    {
        double tot_mass = 0;
        ModuleBase::Vector3<double> tot_momentum;
        for(int i=0; i<numIon; i++)
        {
            tot_mass += allmass[i];
            double sigma = sqrt(temperature / allmass[i]);
            for(int k=0; k<3; ++k)
            {
                if(ionmbl[i][k]==0)
                {
                    vel[i][k] = 0;
                }
                else
                {
                    vel[i][k] = gaussrand() * sigma;
                }

                if(frozen[k] == 0)
                {
                    tot_momentum[k] += allmass[i] * vel[i][k];
                }
            }
        }

        for(int k=0; k<3; ++k)
        {
            if(frozen[k] == 0)
            {
                for(int i=0; i<numIon; i++)
                {
                    vel[i][k] -= tot_momentum[k] / tot_mass;
                }
            }
        }

        double factor = 0.5*(3*numIon-frozen_freedom)*temperature/GetAtomKE(numIon, vel, allmass);
        for(int i=0; i<numIon; i++)
        {
            vel[i] = vel[i]*sqrt(factor);
        }
    }

#ifdef __MPI
    MPI_Bcast(vel,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
    return;
}

void MD_func::InitVel(
	const UnitCell &unit_in, 
	const double& temperature, 
	double* allmass,
	int& frozen_freedom,
	ModuleBase::Vector3<int>* ionmbl,
	ModuleBase::Vector3<double>* vel)
{
	ModuleBase::Vector3<int> frozen;
	getMassMbl(unit_in, allmass, frozen, ionmbl);
	frozen_freedom = frozen.x + frozen.y + frozen.z;
	if(frozen.x == 0) ++frozen_freedom;
	if(frozen.y == 0) ++frozen_freedom;
	if(frozen.z == 0) ++frozen_freedom;

	if(unit_in.init_vel)
    {
        ReadVel(unit_in, vel);
    }
    else
    {
        RandomVel(unit_in.nat, temperature, allmass, frozen_freedom, frozen, ionmbl, vel);
    }
	std::cout << "--------------------------------- INITVEL DONE ------------------------------------" << std::endl;
}

void MD_func::InitPos(
	const UnitCell &unit_in, 
	ModuleBase::Vector3<double>* pos)
{
	int ion=0;
	for(int it=0;it<unit_in.ntype;it++)
	{
		for(int i=0;i<unit_in.atoms[it].na;i++)
		{
			pos[ion] = unit_in.atoms[it].tau[i]*unit_in.lat0;
			ion++;
		}
	}
}

//calculate potential, force and virial
void MD_func::force_virial(
		ModuleESolver::ESolver *p_esolver,
		const int &istep,
		UnitCell &unit_in,
		double &potential,
		ModuleBase::Vector3<double> *force,
		ModuleBase::matrix &virial)
{
	ModuleBase::TITLE("MD_func", "force_stress");
    ModuleBase::timer::tick("MD_func", "force_stress");

    p_esolver->Run(istep, unit_in);

    p_esolver->cal_Energy(potential);

    ModuleBase::matrix force_temp(unit_in.nat, 3); 
    p_esolver->cal_Force(force_temp);

    if(GlobalV::CAL_STRESS)
    {
        p_esolver->cal_Stress(virial);
    }

    // convert Rydberg to Hartree
    potential *= 0.5;
    force_temp *= 0.5;
    virial *= 0.5;

    for(int i=0; i<unit_in.nat; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            force[i][j] = force_temp(i, j);
        }
    }

    ModuleBase::timer::tick("MD_func", "force_stress");
}

void MD_func::outStress(const ModuleBase::matrix &virial, const ModuleBase::matrix &stress)
{
	GlobalV::ofs_running<<"\noutput Pressure for check!"<<std::endl;
    double stress_scalar = 0.0, virial_scalar = 0.0;
    for(int i=0;i<3;i++)
    {
        stress_scalar += stress(i,i)/3;
		virial_scalar += virial(i,i)/3;
    }
    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
    GlobalV::ofs_running<<"Virtual Pressure is "<<stress_scalar*unit_transform<<" Kbar "<<std::endl;
    GlobalV::ofs_running<<"Virial Term is "<<virial_scalar*unit_transform<<" Kbar "<<std::endl;
    GlobalV::ofs_running<<"Kenetic Term is "<<(stress_scalar-virial_scalar)*unit_transform<<" Kbar "<<std::endl;

    GlobalV::ofs_running.unsetf(ios::fixed);
    GlobalV::ofs_running << std::setprecision(8) << std::endl;
    ModuleBase::GlobalFunc::NEW_PART("MD STRESS (KBAR)");
    for (int i=0; i<3; i++)
    {
        GlobalV::ofs_running << std::setw(15) << stress(i,0)*unit_transform 
            << std::setw(15) << stress(i,1)*unit_transform 
            << std::setw(15) << stress(i,2)*unit_transform << std::endl;

    }
    GlobalV::ofs_running << std::setiosflags(ios::left);
}

void MD_func::MDdump(const int &step, 
        const UnitCell &unit_in,
        const ModuleBase::matrix &virial, 
        const ModuleBase::Vector3<double> *force)
{
    if(GlobalV::MY_RANK) return;

    std::stringstream file;
    file << GlobalV::global_out_dir << "MD_dump";
    std::ofstream ofs;
    if(step == 0)
    {
        ofs.open(file.str(), ios::trunc);
    }
    else
    {
        ofs.open(file.str(), ios::app);
    }

    const double unit_virial = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
    const double unit_force = ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU;

    ofs << "MDSTEP:  " << step << std::endl;
    ofs << std::setprecision(12) << std::setiosflags(ios::fixed);

    ofs << "LATTICE_CONSTANT: " << unit_in.lat0 << std::endl;

    ofs << "LATTICE_VECTORS" << std::endl;
    ofs << "  " << unit_in.latvec.e11 << "  " << unit_in.latvec.e12 << "  " << unit_in.latvec.e13 << std::endl; 
    ofs << "  " << unit_in.latvec.e21 << "  " << unit_in.latvec.e22 << "  " << unit_in.latvec.e23 << std::endl;
    ofs << "  " << unit_in.latvec.e31 << "  " << unit_in.latvec.e32 << "  " << unit_in.latvec.e33 << std::endl;

    if(GlobalV::CAL_STRESS)
    {
        ofs << "VIRIAL (KBAR)" << std::endl;
        for(int i=0; i<3; ++i)
        {
            ofs << "  " << virial(i, 0) * unit_virial 
                << "  " << virial(i, 1) * unit_virial 
                << "  " << virial(i, 2) * unit_virial << std::endl;
        }
    }

    ofs << "INDEX    LABEL    POSITIONS    FORCE (eV/Angstrom)" << std::endl;
    int index = 0;
    for(int it=0; it<unit_in.ntype; ++it)
    {
        for(int ia=0; ia<unit_in.atoms[it].na; ++ia)
        {
            ofs << "  " << index
            << "  " << unit_in.atom_label[it]
            << "  " << unit_in.atoms[it].tau[ia].x
            << "  " << unit_in.atoms[it].tau[ia].y
            << "  " << unit_in.atoms[it].tau[ia].z
            << "  " << force[index].x * unit_force 
            << "  " << force[index].y * unit_force 
            << "  " << force[index].z * unit_force << std::endl;
            index++;
        }
    }

    ofs << std::endl;
    ofs << std::endl;
    ofs.close();
}

void MD_func::getMassMbl(const UnitCell &unit_in, 
			double* allmass, 
			ModuleBase::Vector3<int> &frozen,
			ModuleBase::Vector3<int>* ionmbl)
{
//some prepared information
//mass and degree of freedom
	int ion=0;
	frozen.set(0,0,0);
	for(int it=0;it<unit_in.ntype;it++){
		for(int i=0;i<unit_in.atoms[it].na;i++)
		{
			allmass[ion]=unit_in.atoms[it].mass/ModuleBase::AU_to_MASS;
			ionmbl[ion]=unit_in.atoms[it].mbl[i];
			if (ionmbl[ion].x==0) ++frozen.x;
			if (ionmbl[ion].y==0) ++frozen.y;
			if (ionmbl[ion].z==0) ++frozen.z;

			ion++;
		}
	}
}

double MD_func::target_temp(const int &istep, const double &tfirst, const double &tlast)
{
    double delta = (double)(istep) / GlobalV::MD_NSTEP;
    return tfirst + delta * (tlast - tfirst);
}

double MD_func::current_temp(double &kinetic,
            const int &natom, 
            const int &frozen_freedom, 
            const double *allmass,
            const ModuleBase::Vector3<double> *vel)
{
    kinetic = GetAtomKE(natom, vel, allmass);

    return 2 * kinetic / (3 * natom - frozen_freedom);
}

void MD_func::temp_vector(const int &natom, 
            const ModuleBase::Vector3<double> *vel, 
            const double *allmass, 
            ModuleBase::matrix &t_vector)
{
    t_vector.create(3, 3);

    for(int ion=0; ion<natom; ++ion)
    {
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                t_vector(i, j) += allmass[ion] * vel[ion][i] * vel[ion][j];
            }
        }
    }
}