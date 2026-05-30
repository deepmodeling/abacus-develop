# Selective Dynamics for ABACUS+Phonopy

## Requirements

- [ase-abacus](https://gitlab.com/1041176461/ase-abacus)
- [Phonopy](https://github.com/phonopy/phonopy)

## Usage

### Setting

There is a setting file named `config.yaml`:
```yaml
origin_structure: 'STRU'     # the original structure
selected_indices: [18, 19]   # the atom index that you concerned about, start from 0
tasks_per_batch: 12          # how much jobs per batch
wait_time: 600               # s, sleep time between batches

setting_conf: |
  SYMMETRY = .FALSE. 
  DIM = 1 1 1 
  DISPLACEMENT_DISTANCE = 0.03

mesh.conf: |
  DIM = 1 1 1
  MESH = 31 31 31
  TMAX = 2000
  TSTEP = 2

input: |
  INPUT_PARAMETERS
  #Parameters (1.General)
  suffix                 phonon
  calculation            scf
  symmetry               1
  nspin                  1
  pseudo_dir             /fs2/home/chenkg/2_liuyu/3_abacus/PP_ORB/pseudo
  orbital_dir            /fs2/home/chenkg/2_liuyu/3_abacus/PP_ORB/efficiency
  kpoint_file            ../KPT
  
  #Parameters (2.Iteration)
  ecutwfc                100
  scf_thr                1e-8
  scf_nmax               100
  
  #Parameters (3.Basis)
  basis_type             lcao
  ks_solver              genelpa
  gamma_only             0
  
  #Parameters (4.Smearing)
  smearing_method        gaussian
  smearing_sigma         0.001
  
  #Parameters (5.Mixing)
  mixing_type            broyden
  mixing_beta            0.7
  
  cal_force              1
  cal_stress             1

kpt: |
  K_POINTS
  0
  Gamma
  5 5 1 0 0 0

job_script: |
  #!/bin/bash
  #SBATCH -p cp6
  #SBATCH -N 1
  #SBATCH -J abacus
  #SBATCH -n 28

  source /fs2/home/chenkg/2_liuyu/3_abacus/abacus_env.sh
  export OMP_NUM_THREADS=28

  mpirun -n 1 abacus
```

- origin_structure: The `STRU` filename, which contains both the fixed atoms and the free atoms.
- selected_indices: The indexs of the free atoms. Note that the index starts from 0.
- tasks_per_batch: How much jobs submitted per batch.
- wait_time: Sleep time between batches, the unit is second.
- setting_conf: The `setting.conf` file for Phonopy.
- mesh.conf: The `mesh.conf` file for Phonopy.
- input: The `INPUT` file for ABACUS.
- kpt: The `KPT` file for ABACUS.
- job_script: The script used to submit jobs.

### Submit jobs

Use the following command
```bash
python3 path_to_selective_dynamics.py --submit
```
to generate displaced structures and submit jobs.

### Postprocess

Use the following command
```bash
python3 path_to_selective_dynamics.py --post
```
to generate `FORCE_SETS` and results of phonon calculations.
