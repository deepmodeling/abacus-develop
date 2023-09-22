#ifndef MD_FUNC_H
#define MD_FUNC_H

#include "module_esolver/esolver.h"

/**
 * @brief base functions in md
 *
 * Namespace MD_func contains several functions used in md.
 */
namespace MD_func
{

/**
 * @brief generate a Gaussian random number
 *
 * @return a Gaussian random number
 */
double gaussrand();

/**
 * @brief initialize the atomic velocities
 *
 * @param unit_in unitcell information
 * @param my_rank MPI rank of the processor
 * @param restart whether restart the md
 * @param temperature ion temperature
 * @param allmass atomic mass
 * @param frozen_freedom the fixed freedom
 * @param ionmbl determine whether the atomic freedom is fixed
 * @param vel the genarated atomic velocities
 */
void init_vel(const UnitCell& unit_in,
              const int& my_rank,
              const bool& restart,
              double& temperature,
              double* allmass,
              int& frozen_freedom,
              ModuleBase::Vector3<int>* ionmbl,
              ModuleBase::Vector3<double>* vel);

/**
 * @brief read in atomic velocities from STRU
 *
 * @param unit_in unitcell information
 * @param vel the read-in atomic velocities
 */
void read_vel(const UnitCell& unit_in, ModuleBase::Vector3<double>* vel);

/**
 * @brief generate atomic velocities that satisfy the Boltzmann distribution
 *
 * @param natom the number of atoms
 * @param temperature ion temperature
 * @param allmass atomic mass
 * @param frozen_freedom the fixed freedom
 * @param frozen the fixed freedom along three directions
 * @param ionmbl determine whether the atomic freedom is fixed
 * @param my_rank MPI rank of the processor
 * @param vel the genarated atomic velocities
 */
void rand_vel(const int& natom,
              const double& temperature,
              const double* allmass,
              const int& frozen_freedom,
              const ModuleBase::Vector3<int> frozen,
              const ModuleBase::Vector3<int>* ionmbl,
              const int& my_rank,
              ModuleBase::Vector3<double>* vel);

/**
 * @brief calculate energy, forces and virial tensor
 *
 * @param p_esolver enrergy solver
 * @param istep current md step
 * @param unit_in unitcell information
 * @param potential potential energy
 * @param force atomic forces
 * @param cal_stress whether calculate stress
 * @param virial lattice virial tensor
 */
void force_virial(ModuleESolver::ESolver* p_esolver,
                  const int& istep,
                  UnitCell& unit_in,
                  double& potential,
                  ModuleBase::Vector3<double>* force,
                  const bool& cal_stress,
                  ModuleBase::matrix& virial);
/**
 * @brief calculate the ionic kinetic energy
 *
 * @param natom the number of atoms
 * @param vel the atomic velocities
 * @param allmass atomic mass
 * @return the ionic kinetic energy
 */
double kinetic_energy(const int& natom, const ModuleBase::Vector3<double>* vel, const double* allmass);

/**
 * @brief calculate the total stress tensor
 *
 * @param unit_in unitcell information
 * @param vel the atomic velocities
 * @param allmass atomic mass
 * @param cal_stress whether calculate stress
 * @param virial lattice virial tensor
 * @param stress total stress tensor
 */
void compute_stress(const UnitCell& unit_in,
                    const ModuleBase::Vector3<double>* vel,
                    const double* allmass,
                    const bool& cal_stress,
                    const ModuleBase::matrix& virial,
                    ModuleBase::matrix& stress);

/**
 * @brief output the stress information
 *
 * @param ofs determine the output files
 * @param virial lattice virial tensor
 * @param stress total stress tensor
 */
void print_stress(std::ofstream& ofs, const ModuleBase::matrix& virial, const ModuleBase::matrix& stress);

/**
 * @brief dump the md information
 *
 * including md step, lattice constant, lattice vectors, lattice virial tensor, ion index, ion positions,
 * ion velocities, ion forces
 *
 * @param step current md step
 * @param global_out_dir directory of output files
 * @param unit_in unitcell information
 * @param mdp input parameters used in md
 * @param virial lattice virial tensor
 * @param force atomic forces
 * @param vel atomic velocities
 */
void dump_info(const int& step,
               const std::string& global_out_dir,
               const UnitCell& unit_in,
               const MD_para& mdp,
               const ModuleBase::matrix& virial,
               const ModuleBase::Vector3<double>* force,
               const ModuleBase::Vector3<double>* vel);

/**
 * @brief obtain the atomic mass and whether the freedom is fixed
 *
 * @param unit_in unitcell information
 * @param allmass atomic mass
 * @param frozen the fixed freedom along three directions
 * @param ionmbl determine whether the atomic freedom is fixed
 */
void get_mass_mbl(const UnitCell& unit_in,
                  double* allmass,
                  ModuleBase::Vector3<int>& frozen,
                  ModuleBase::Vector3<int>* ionmbl);

/**
 * @brief get the target temperature of the current md step
 *
 * @param istep the current md step
 * @param nstep the total md step
 * @param tfirst the initial temperature
 * @param tlast the final temperature
 * @return the target temperature
 */
double target_temp(const int& istep, const int& nstep, const double& tfirst, const double& tlast);

/**
 * @brief get the current temperature
 *
 * @param kinetic kinetic energy
 * @param natom the number of atoms
 * @param frozen_freedom the fixed freedom
 * @param allmass atomic mass
 * @param vel atomic velocities
 * @return the current temperature
 */
double current_temp(double& kinetic,
                    const int& natom,
                    const int& frozen_freedom,
                    const double* allmass,
                    const ModuleBase::Vector3<double>* vel);

/**
 * @brief get the temperature vectors
 *
 * @param natom the number of atoms
 * @param vel atomic velocities
 * @param allmass atomic mass
 * @param t_vector the temperature vectors
 */
void temp_vector(const int& natom,
                 const ModuleBase::Vector3<double>* vel,
                 const double* allmass,
                 ModuleBase::matrix& t_vector);

/**
 * @brief determine thr current md step
 *
 * determine thr current md step according to Restart_md.dat if md_restart is true
 *
 * @param my_rank MPI rank of the processor
 * @param file_dir the directory of read-in files
 * @return thr current md step
 */
double current_step(const int& my_rank, const std::string& file_dir);

} // namespace MD_func

#endif // MD_FUNC_H
