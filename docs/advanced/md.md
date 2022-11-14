# Molecular Dynamics

Molecular dynamics (MD) is a computer simulation method for analyzing the physical movements of atoms and molecules. The atoms and molecules are allowed to interact for a fixed period of time, giving a view of the dynamic "evolution" of the system. In the most common version, the trajectories of atoms and molecules are determined by numerically solving Newton's equations of motion for a system of interacting particles, where forces between the particles and their potential energies are calculated using first-principles calculations (first-principles molecular dynamics, FPMD), or interatomic potentials and molecular mechanics force fields (classical molecular dynamics, CMD).

By setting `calculation` to be `md`, ABACUS currently provides several different MD evolution methods, which is specified by keyword `md_type` in the `INPUT` file:

  - -1: FIRE method
  - 0: velocity Verlet algorithm (default: NVE ensemble)
  - 1: Nose-Hoover style non-Hamiltonian equations of motion
  - 2: NVT ensemble with Langevin thermostat
  - 4: MSST method

When `md_type` is set to 0, `md_thermostat` is used to specify the thermostat based on the velocity Verlet algorithm.

  - nve: NVE ensemble
  - anderson: NVT ensemble with Anderson thermostat
  - berendsen: NVT ensemble with Berendsen thermostat
  - rescaling: NVT ensemble with velocity Rescaling method 1
  - rescale_v: NVT ensemble with velocity Rescaling method 2

When `md_type` is set to 1, `md_pmode` is used to specify the NVT or NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.

  - none: NVT ensemble
  - iso: NPT ensemble with isotropic cetl fluctuations
  - aniso: NPT ensemble with anisotropic cetl fluctuations
  - tri: NPT ensemble with non-orthogonal (triclinic) simulation box

Furthermore, ABACUS also provides a [list of keywords](./input_files/input-main.md#molecular-dynamics) to control relevant parmeters used in MD simulations.

To employ CMD calculations, `esolver_type` should be set to be `lj` or `dp`.
If DP model is selected, the filename of DP model is specified by keyword `pot_file`.

[Examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/md/lcao_gammaonly_Sn64) of MD simulations are also provided.
There are eight INPUT files corresponding to eight different MD evolution methods in the directory.
For examlpe, `INPUT_0` shows how to employ the NVE simulation.

To run any of the fix cases, users may enter the directory, copy the corresponding input file to `INPUT`, and run ABACUS.

## FIRE
[FIRE](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201) (fast inertial relaxation engine) is a MD-based minimization algorithm. It is based on conventional molecular dynamics with additional velocity modifications and adaptive time steps. The MD trajectory will descend to an energy-minimum.

## NVE

NVE ensemble (i. e. microcanonical ensemble) is a statistical ensemble that represents the possible states of a mechanical system whose total energy is exactly specified. The system is assumed to be isolated in the sense that it cannot exchange energy or particles with its environment, so that the energy of the system does not change with time.

The primary macroscopic variables of the microcanonical ensemble are the total number of particles in the system (symbol: N), the system's volume (symbol: V), as well as the total energy in the system (symbol: E). Each of these is assumed to be constant in the ensemble.

Currently NVE ensemble in ABACUS is implemented based on the [velocity verlet algorithm](https://aip.scitation.org/doi/abs/10.1063/1.442716).


## Nose Hoover Chain

NVT ensemble (i. e. canonical ensemble)  is the statistical ensemble that represents the possible states of a mechanical system in thermal equilibrium with a heat bath at a fixed temperature. The system can exchange energy with the heat bath, so that the states of the system will differ in total energy.

The principal thermodynamic variable of the canonical ensemble, determining the probability distribution of states, is the absolute temperature (symbol: T). The ensemble typically also depends on mechanical variables such as the number of particles in the system (symbol: N) and the system's volume (symbol: V), each of which influence the nature of the system's internal states. An ensemble with these three parameters is sometimes called the NVT ensemble.


The isothermal–isobaric ensemble (constant temperature and constant pressure ensemble), also called NPT ensemble, is a statistical mechanical ensemble that maintains the number of particles N, constant temperature T, and constant pressure P. This ensemble plays an important role in chemistry as chemical reactions are usually carried out under constant pressure condition. The NPT ensemble is also useful for measuring the equation of state of model systems whose virial expansion for pressure cannot be evaluated, or systems near first-order phase transitions.

ABACUS perform time integration on [Nose-Hoover style non-Hamiltonian equations of motion](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.31.1695) which are designed to generate positions and velocities sampled from NVT and NPT ensemble.



## Langevin

[Langevin thermostat](https://en.wikipedia.org/wiki/Langevin_dynamics) can be used for molecular dynamics equations by assuming that the atoms being simulated are embedded in a sea of much smaller fictional particles. In many instances of solute-solvent systems, the behavior of the solute is desired, and the behavior of the solvent is non-interesting(e.g. proteins, DNA, nanoparticles in solution). In these cases, the solvent influences the dynamics of the solute(typically nanoparticles) via random collisions, and by imposing a frictional drag force on the motion of the nanoparticle in the solvent. The damping factor and the random force combine to give the correct NVT ensemble.

## Anderson

[Anderson thermostat](https://aip.scitation.org/doi/abs/10.1063/1.439486) couples the system to a heat bath that imposes the desired temperature to simulate the NVT ensemble. The coupling to a heat bath is represented by stochastic collision that act occasionally on randomly selected particles.

## Berendsen

Reset the temperature of a group of atoms by using a [Berendsen thermostat](https://aip.scitation.org/doi/10.1063/1.448118), which rescales their velocities every timestep. In this scheme, the system is weakly coupled to a heat bath with some temperature. Though the thermostat does not generate a correct canonical ensemble (especially for small systems), for large systems on the order of hundreds or thousands of atoms/molecules, the approximation yields roughly correct results for most calculated properties.

## Rescaling

Reset the temperature of a group of atoms by explicitly rescaling their velocities. Velocities are rescaled if the current and target temperature differ more than `md_tolerance` (Kelvin).

## Rescale_v

Reset the temperature of a group of atoms by explicitly rescaling their velocities. Every `md_nraise` steps the current temperature is rescaled to target temperature.

## MSST
ABACUS performs the [Multi-Scale Shock Technique (MSST) integration](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.90.235503) to update positions and velocities each timestep to mimic a compressive shock wave passing over the system. The MSST varies the cell volume and temperature in such a way as to restrain the system to the shock Hugoniot and the Rayleigh line. These restraints correspond to the macroscopic conservation laws dictated by a shock front.