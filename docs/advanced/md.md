# Molecular Dynamics

Molecular dynamics (MD) is a computer simulation method for analyzing the physical movements of atoms and molecules. The atoms and molecules are allowed to interact for a fixed period of time, giving a view of the dynamic "evolution" of the system. In the most common version, the trajectories of atoms and molecules are determined by numerically solving Newton's equations of motion for a system of interacting particles, where forces between the particles and their potential energies are calculated using first-principles calculations (first-principles molecular dynamics, FPMD), or interatomic potentials and molecular mechanics force fields (classical molecular dynamics, CMD).

By setting [calculation](./input_files/input-main.md#calculation) to be `md`, ABACUS currently provides several different MD evolution methods, which is specified by keyword [md_type](./input_files/input-main.md#md_type) in the `INPUT` file:

  - fire: a MD-based relaxation algorithm, see [details](#fire) here
  - nve: NVE ensemble with velocity Verlet algorithm
  - nvt: NVT ensemble
  - npt: Nose-Hoover style NPT ensemble
  - langevin: NVT ensemble with Langevin thermostat
  - msst: MSST method

When [md_type](./input_files/input-main.md#md_type) is set to nvt, [md_thermostat](./input_files/input-main.md#md_thermostat) is used to specify the temperature control method used in NVT ensemble.

  - nhc: Nose-Hoover chain
  - anderson: Anderson thermostat
  - berendsen: Berendsen thermostat
  - rescaling: velocity Rescaling method 1
  - rescale_v: velocity Rescaling method 2

When [md_type](./input_files/input-main.md#md_type) is set to npt, [md_pmode](./input_files/input-main.md#md_pmode) is used to specify the cell fluctuation mode in NPT ensemble based on the Nose-Hoover style non-Hamiltonian equations of motion.

  - iso: isotropic cell fluctuations
  - aniso: anisotropic cell fluctuations
  - tri: non-orthogonal (triclinic) simulation box

Furthermore, ABACUS also provides a [list of keywords](./input_files/input-main.md#molecular-dynamics) to control relevant parmeters used in MD simulations.

The MD output information will be written into the file `MD_dump`， in which the atomic forces, atomic velocities, and lattice virial are controlled by keyword [dump_force](./input_files/input-main.md#dump_force), [dump_vel](./input_files/input-main.md#dump_vel), and [dump_virial](./input_files/input-main.md#dump_virial), respectively.

[Examples](../../examples/md/lcao_gammaonly_Si8/) of MD simulations are also provided.
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

Reset the temperature of a group of atoms by explicitly rescaling their velocities. Velocities are rescaled if the current and target temperature differ more than [md_tolerance](./input_files/input-main.md#md_tolerance) (Kelvin).

## Rescale_v

Reset the temperature of a group of atoms by explicitly rescaling their velocities. Every [md_nraise](./input_files/input-main.md#md_nraise) steps the current temperature is rescaled to target temperature.

## MSST
ABACUS performs the [Multi-Scale Shock Technique (MSST) integration](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.90.235503) to update positions and velocities each timestep to mimic a compressive shock wave passing over the system. The MSST varies the cell volume and temperature in such a way as to restrain the system to the shock Hugoniot and the Rayleigh line. These restraints correspond to the macroscopic conservation laws dictated by a shock front.

## DPMD
Compiling ABACUS with [DeePMD-kit](https://github.com/deepmodeling/deepmd-kit), MD calculations based on machine learning DP model is enabled.

To employ DPMD calculations, [esolver_type](./input_files/input-main.md#esolver_type) should be set to `dp`.
And the filename of DP model is specified by keyword [pot_file](./input_files/input-main.md#pot_file).

First, we can find whether contains keyword `type_map` in the DP model through the shell command:
```bash
strings Al-SCAN.pb | grep type_map
```

```json
{"model": {"type_map": ["Al"], "descriptor": {"type": "se_e2_a", "sel": [150], "rcut_smth": 0.5, "rcut": 6.0, "neuron": [25, 50, 100], "resnet_dt": false, "axis_neuron": 16, "seed": 1, "activation_function": "tanh", "type_one_side": false, "precision": "default", "trainable": true, "exclude_types": [], "set_davg_zero": false}, "fitting_net": {"neuron": [240, 240, 240], "resnet_dt": true, "seed": 1, "type": "ener", "numb_fparam": 0, "numb_aparam": 0, "activation_function": "tanh", "precision": "default", "trainable": true, "rcond": 0.001, "atom_ener": []}, "data_stat_nbatch": 10, "data_stat_protect": 0.01}, "learning_rate": {"type": "exp", "decay_steps": 5000, "start_lr": 0.001, "stop_lr": 3.51e-08, "scale_by_worker": "linear"}, "loss": {"type": "ener", "start_pref_e": 0.02, "limit_pref_e": 1, "start_pref_f": 1000, "limit_pref_f": 1, "start_pref_v": 0, "limit_pref_v": 0, "start_pref_ae": 0.0, "limit_pref_ae": 0.0, "start_pref_pf": 0.0, "limit_pref_pf": 0.0, "enable_atom_ener_coeff": false}, "training": {"training_data": {"systems": ["../deepmd_data/"], "batch_size": "auto", "set_prefix": "set", "auto_prob": "prob_sys_size", "sys_probs": null}, "validation_data": {"systems": ["../deepmd_validation"], "batch_size": 1, "numb_btch": 3, "set_prefix": "set", "auto_prob": "prob_sys_size", "sys_probs": null}, "numb_steps": 1000000, "seed": 10, "disp_file": "lcurve.out", "disp_freq": 100, "save_freq": 1000, "save_ckpt": "model.ckpt", "disp_training": true, "time_training": true, "profiling": false, "profiling_file": "timeline.json", "enable_profiler": false, "tensorboard": false, "tensorboard_log_dir": "log", "tensorboard_freq": 1}}
```

If the keyword `type_map` is found, ABACUS will match the atom types between `STRU` and DP model.

Otherwise, all atom types must be specified in the `STRU` in the order consistent with that of the DP model, even if the number of atoms is zero!

For example, there is a Al-Cu-Mg ternary-alloy DP model, but the simulated cell is a Al-Cu binary alloy. Then the `STRU` should be written as follows:

```
ATOMIC_SPECIES
Al  26.982
Cu  63.546
Mg  24.305

LATTICE_CONSTANT
1.889727000000

LATTICE_VECTORS
4.0  0.0  0.0
0.0  4.0  0.0
0.0  0.0  4.0

ATOMIC_POSITIONS
Cartesian

Al
0
2
0.0  0.0  0.0
0.5  0.5  0.0

Cu
0
2
0.5  0.0  0.5
0.0  0.5  0.5

Mg
0
0
```