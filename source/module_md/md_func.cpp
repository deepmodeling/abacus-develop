#include "md_func.h"

#include "module_base/global_variable.h"
#include "module_base/timer.h"

namespace MD_func
{

double gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if (phase == 0)
    {
        do
        {
            double U1 = static_cast<double>(std::rand()) / RAND_MAX;
            double U2 = static_cast<double>(std::rand()) / RAND_MAX;

            V1 = 2.0 * U1 - 1.0;
            V2 = 2.0 * U2 - 1.0;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1 || S == 0);

        X = V1 * sqrt(-2.0 * log(S) / S);
    }
    else
    {
        X = V2 * sqrt(-2.0 * log(S) / S);
    }

    phase = 1 - phase;

    return X;
}

double kinetic_energy(const int& natom, const ModuleBase::Vector3<double>* vel, const double* allmass)
{
    double ke = 0;

    for (int ion = 0; ion < natom; ++ion)
    {
        ke += 0.5 * allmass[ion] * vel[ion].norm2();
    }

    return ke;
}

void compute_stress(const UnitCell& unit_in,
                    const ModuleBase::Vector3<double>* vel,
                    const double* allmass,
                    const bool& cal_stress,
                    const ModuleBase::matrix& virial,
                    ModuleBase::matrix& stress)
{
    if (cal_stress)
    {
        ModuleBase::matrix t_vector;

        temp_vector(unit_in.nat, vel, allmass, t_vector);

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                stress(i, j) = virial(i, j) + t_vector(i, j) / unit_in.omega;
            }
        }
    }
}

void read_vel(const UnitCell& unit_in, ModuleBase::Vector3<double>* vel)
{
    int iat = 0;
    for (int it = 0; it < unit_in.ntype; ++it)
    {
        for (int ia = 0; ia < unit_in.atoms[it].na; ++ia)
        {
            vel[iat] = unit_in.atoms[it].vel[ia];
            if (unit_in.atoms[it].mbl[ia].x == 0)
                vel[iat].x = 0;
            if (unit_in.atoms[it].mbl[ia].y == 0)
                vel[iat].y = 0;
            if (unit_in.atoms[it].mbl[ia].z == 0)
                vel[iat].z = 0;
            ++iat;
        }
    }
    assert(iat == unit_in.nat);
}

void rand_vel(const int& natom,
              const double& temperature,
              const double* allmass,
              const int& frozen_freedom,
              const ModuleBase::Vector3<int> frozen,
              const ModuleBase::Vector3<int>* ionmbl,
              const int& my_rank,
              ModuleBase::Vector3<double>* vel)
{
    if (!my_rank)
    {
        double tot_mass = 0;
        ModuleBase::Vector3<double> tot_momentum;
        for (int i = 0; i < natom; i++)
        {
            tot_mass += allmass[i];
            double sigma = sqrt(temperature / allmass[i]);
            for (int k = 0; k < 3; ++k)
            {
                if (ionmbl[i][k] == 0)
                {
                    vel[i][k] = 0;
                }
                else
                {
                    vel[i][k] = gaussrand() * sigma;
                }

                if (frozen[k] == 0)
                {
                    tot_momentum[k] += allmass[i] * vel[i][k];
                }
            }
        }

        for (int k = 0; k < 3; ++k)
        {
            if (frozen[k] == 0)
            {
                for (int i = 0; i < natom; i++)
                {
                    vel[i][k] -= tot_momentum[k] / tot_mass;
                }
            }
        }

        double factor;
        if (3 * natom == frozen_freedom || temperature == 0)
        {
            factor = 0;
        }
        else
        {
            factor = 0.5 * (3 * natom - frozen_freedom) * temperature / kinetic_energy(natom, vel, allmass);
        }

        for (int i = 0; i < natom; i++)
        {
            vel[i] = vel[i] * sqrt(factor);
        }
    }

#ifdef __MPI
    MPI_Bcast(vel, natom * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    return;
}

void init_vel(const UnitCell& unit_in,
              const int& my_rank,
              const bool& restart,
              double& temperature,
              double* allmass,
              int& frozen_freedom,
              ModuleBase::Vector3<int>* ionmbl,
              ModuleBase::Vector3<double>* vel)
{
    std::cout << " ----------------------------------- INIT VEL ---------------------------------------" << std::endl;
    ModuleBase::Vector3<int> frozen;
    get_mass_mbl(unit_in, allmass, frozen, ionmbl);
    frozen_freedom = frozen.x + frozen.y + frozen.z;
    if (frozen.x == 0)
        ++frozen_freedom;
    if (frozen.y == 0)
        ++frozen_freedom;
    if (frozen.z == 0)
        ++frozen_freedom;

    if (unit_in.init_vel)
    {
        std::cout << " READ VEL FROM STRU" << std::endl;
        read_vel(unit_in, vel);
        double kinetic = 0.0;
        double t_current = MD_func::current_temp(kinetic, unit_in.nat, frozen_freedom, allmass, vel);
        if (restart)
        {
            std::cout << " RESTART MD, CURRENT TEMPERATURE IS " << t_current * ModuleBase::Hartree_to_K << " K"
                      << std::endl;
        }
        else if (temperature < 0)
        {
            std::cout << " UNSET INITIAL TEMPERATURE, AUTOSET TO " << t_current * ModuleBase::Hartree_to_K << " K"
                      << std::endl;
            temperature = t_current;
        }
        else if (std::fabs(temperature - t_current) > 1e-8)
        {
            std::cout << " INITIAL TEMPERATURE IN INPUT  = " << temperature * ModuleBase::Hartree_to_K << " K"
                      << std::endl;
            std::cout << " READING TEMPERATURE FROM STRU = " << t_current * ModuleBase::Hartree_to_K << " K"
                      << std::endl;
            std::cout << " INCONSISTENCE, PLEASE CHECK" << std::endl;
            std::cout << " ------------------------------------- DONE -----------------------------------------"
                      << std::endl;
            exit(0);
        }
    }
    else
    {
        std::cout << " RANDOM VEL ACCORDING TO INITIAL TEMPERATURE: " << temperature * ModuleBase::Hartree_to_K << " K"
                  << std::endl;
        rand_vel(unit_in.nat, temperature, allmass, frozen_freedom, frozen, ionmbl, my_rank, vel);
    }
    std::cout << " ------------------------------------- DONE -----------------------------------------" << std::endl;
}

void force_virial(ModuleESolver::ESolver* p_esolver,
                  const int& istep,
                  UnitCell& unit_in,
                  double& potential,
                  ModuleBase::Vector3<double>* force,
                  const bool& cal_stress,
                  ModuleBase::matrix& virial)
{
    ModuleBase::TITLE("MD_func", "force_virial");
    ModuleBase::timer::tick("MD_func", "force_virial");

    p_esolver->Run(istep, unit_in);

    potential = p_esolver->cal_Energy();

    ModuleBase::matrix force_temp(unit_in.nat, 3);
    p_esolver->cal_Force(force_temp);

    if (cal_stress)
    {
        p_esolver->cal_Stress(virial);
    }

    /// convert Rydberg to Hartree
    potential *= 0.5;
    force_temp *= 0.5;
    virial *= 0.5;

    for (int i = 0; i < unit_in.nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            force[i][j] = force_temp(i, j);
        }
    }

    ModuleBase::timer::tick("MD_func", "force_virial");
}

void print_stress(std::ofstream& ofs, const ModuleBase::matrix& virial, const ModuleBase::matrix& stress)
{
    double stress_scalar = 0.0, virial_scalar = 0.0;
    for (int i = 0; i < 3; i++)
    {
        stress_scalar += stress(i, i) / 3;
        virial_scalar += virial(i, i) / 3;
    }
    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    ofs << "Virtual Pressure is " << stress_scalar * unit_transform << " kbar " << std::endl;
    ofs << "Virial Term is " << virial_scalar * unit_transform << " kbar " << std::endl;
    ofs << "Kinetic Term is " << (stress_scalar - virial_scalar) * unit_transform << " kbar " << std::endl;

    ofs.unsetf(std::ios::fixed);
    ofs << std::setprecision(8) << std::endl;
    ModuleBase::GlobalFunc::NEW_PART("MD STRESS (kbar)");
    for (int i = 0; i < 3; i++)
    {
        ofs << std::setw(15) << stress(i, 0) * unit_transform << std::setw(15) << stress(i, 1) * unit_transform
            << std::setw(15) << stress(i, 2) * unit_transform << std::endl;
    }
    ofs << std::setiosflags(std::ios::left);
}

void dump_info(const int& step,
               const std::string& global_out_dir,
               const UnitCell& unit_in,
               const MD_para& mdp,
               const ModuleBase::matrix& virial,
               const ModuleBase::Vector3<double>* force,
               const ModuleBase::Vector3<double>* vel)
{
    if (mdp.my_rank)
        return;

    std::stringstream file;
    file << global_out_dir << "MD_dump";
    std::ofstream ofs;
    if (step == 0)
    {
        ofs.open(file.str(), std::ios::trunc);
    }
    else
    {
        ofs.open(file.str(), std::ios::app);
    }

    const double unit_pos = unit_in.lat0 / ModuleBase::ANGSTROM_AU;                                  ///< Angstrom
    const double unit_vel = 1.0 / ModuleBase::ANGSTROM_AU / ModuleBase::AU_to_FS;                    ///< Angstrom/fs
    const double unit_virial = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8; ///< kBar
    const double unit_force = ModuleBase::Hartree_to_eV * ModuleBase::ANGSTROM_AU;                   ///< eV/Angstrom

    ofs << "MDSTEP:  " << step << std::endl;
    ofs << std::setprecision(12) << std::setiosflags(std::ios::fixed);

    ofs << "LATTICE_CONSTANT: " << unit_in.lat0_angstrom << " Angstrom" << std::endl;

    ofs << "LATTICE_VECTORS" << std::endl;
    ofs << "  " << unit_in.latvec.e11 << "  " << unit_in.latvec.e12 << "  " << unit_in.latvec.e13 << std::endl;
    ofs << "  " << unit_in.latvec.e21 << "  " << unit_in.latvec.e22 << "  " << unit_in.latvec.e23 << std::endl;
    ofs << "  " << unit_in.latvec.e31 << "  " << unit_in.latvec.e32 << "  " << unit_in.latvec.e33 << std::endl;

    if (mdp.cal_stress && mdp.dump_virial)
    {
        ofs << "VIRIAL (kbar)" << std::endl;
        for (int i = 0; i < 3; ++i)
        {
            ofs << "  " << virial(i, 0) * unit_virial << "  " << virial(i, 1) * unit_virial << "  "
                << virial(i, 2) * unit_virial << std::endl;
        }
    }

    ofs << "INDEX    LABEL    POSITION (Angstrom)";
    if (mdp.dump_force)
    {
        ofs << "    FORCE (eV/Angstrom)";
    }
    if (mdp.dump_vel)
    {
        ofs << "    VELOCITY (Angstrom/fs)";
    }
    ofs << std::endl;

    int index = 0;
    for (int it = 0; it < unit_in.ntype; ++it)
    {
        for (int ia = 0; ia < unit_in.atoms[it].na; ++ia)
        {
            ofs << "  " << index << "  " << unit_in.atom_label[it] << "  " << unit_in.atoms[it].tau[ia].x * unit_pos
                << "  " << unit_in.atoms[it].tau[ia].y * unit_pos << "  " << unit_in.atoms[it].tau[ia].z * unit_pos;

            if (mdp.dump_force)
            {
                ofs << "  " << force[index].x * unit_force << "  " << force[index].y * unit_force << "  "
                    << force[index].z * unit_force;
            }

            if (mdp.dump_vel)
            {
                ofs << "  " << vel[index].x * unit_vel << "  " << vel[index].y * unit_vel << "  "
                    << vel[index].z * unit_vel;
            }
            ofs << std::endl;
            index++;
        }
    }

    ofs << std::endl;
    ofs << std::endl;
    ofs.close();
}

void get_mass_mbl(const UnitCell& unit_in,
                  double* allmass,
                  ModuleBase::Vector3<int>& frozen,
                  ModuleBase::Vector3<int>* ionmbl)
{
    int ion = 0;
    frozen.set(0, 0, 0);
    for (int it = 0; it < unit_in.ntype; it++)
    {
        for (int i = 0; i < unit_in.atoms[it].na; i++)
        {
            allmass[ion] = unit_in.atoms[it].mass / ModuleBase::AU_to_MASS;
            ionmbl[ion] = unit_in.atoms[it].mbl[i];
            if (ionmbl[ion].x == 0)
                ++frozen.x;
            if (ionmbl[ion].y == 0)
                ++frozen.y;
            if (ionmbl[ion].z == 0)
                ++frozen.z;

            ion++;
        }
    }
}

double target_temp(const int& istep, const int& nstep, const double& tfirst, const double& tlast)
{
    double delta = static_cast<double>(istep) / nstep;
    return tfirst + delta * (tlast - tfirst);
}

double current_temp(double& kinetic,
                    const int& natom,
                    const int& frozen_freedom,
                    const double* allmass,
                    const ModuleBase::Vector3<double>* vel)
{
    if (3 * natom == frozen_freedom)
    {
        return 0;
    }
    else
    {
        kinetic = kinetic_energy(natom, vel, allmass);
        return 2 * kinetic / (3 * natom - frozen_freedom);
    }
}

void temp_vector(const int& natom,
                 const ModuleBase::Vector3<double>* vel,
                 const double* allmass,
                 ModuleBase::matrix& t_vector)
{
    t_vector.create(3, 3);

    for (int ion = 0; ion < natom; ++ion)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                t_vector(i, j) += allmass[ion] * vel[ion][i] * vel[ion][j];
            }
        }
    }
}

void current_md_info(const int& my_rank, const std::string& file_dir, int& md_step, double& temperature)
{
    bool ok = true;

    if (my_rank == 0)
    {
        std::stringstream ssc;
        ssc << file_dir << "Restart_md.dat";
        std::ifstream file(ssc.str().c_str());

        if (!file)
        {
            ok = false;
        }

        if (ok)
        {
            file >> md_step >> temperature;
            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    if (!ok)
    {
        ModuleBase::WARNING_QUIT("current_md_info", "no Restart_md.dat!");
    }

#ifdef __MPI
    MPI_Bcast(&md_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

} // namespace MD_func