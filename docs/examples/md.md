# Molecular dynamics

[back to main page](../../README.md)

A typical INPUT file for MD simulation is given
below, which contains additional keywords that need to added.
```
INPUT_PARAMETERS
#Parameters (1.General)
pseudo_dir          ./
ntype               1
nbands              600
pseudo_type         upf
gamma_only          1
calculation         md
symmetry            0
out_level           m

#Parameters (2.SCF)
ecutwfc             60
scf_thr             1e-6
scf_nmax            100

#Parameters (3.Basis)
basis_type          lcao
ks_solver           genelpa

#Parameters (4.Smearing)
smearing_method     gaussian
smearing_sigma      0.02

#Parameters (5.Mixing)
mixing_type         pulay
mixing_beta         0.4
chg_extrap          second-order

#Parameters (6.MD)
md_nstep            10   // md steps
md_type             1    // choose ensemble
md_dt               1    // time step
md_tfirst           700  // the first target temperature
md_restart          0    // whether restart md
md_dumpfreq         10   // The period to dump MD information
```

These MD parameters means that ABACUS will use NVT ensemble with Nosé-hoover chain themostat; the time step is 1fs, and target temperature is 700K; start new MD without restart file.

Note: *Please turn off symmetry during MD simulations.*

- md_type : 
  - -1, FIRE : an ionic relaxation method using MD;
  - 0, NVE : microcanonical ensemble; 
  - 1, NHC : canonical ensemble (NVT) with Nosé-hoover chain themostat ; 
  - 2, LGV : ion dynamics is over-damped Langevin, modeling an interaction with a
background implicit solvent.; 
  - 3, ADS : NVT, control ionic temperature using Andersen thermostat; 
  - 4, MSST :  Multi-Scale Shock Technique (MSST) integration to mimic a compressive shock wave passing over the system.
- md_dt : time step in md simulation (fs)
- md_tfirst : target temperature in md simulation (K), you should set parameter md_tlast when you want to change temperature during md simulation.
- md_restart : 0, no need of restart ; 1, restart with restart file, you must repalce STRU file with STRU_MD before you run the restart task.
- md_dumpfreq : frequency for output consequence of md simulations.

The STRU file is:
```
ATOMIC_SPECIES
Sn 118.71  ./Sn.pz-bhs.UPF 

LATTICE_CONSTANT
23.45378

NUMERICAL_ORBITAL
./Sn_pz-bhs_8.0au_16Ry_2s2p1d

LATTICE_VECTORS
1.00     0.00     0.00
0.00     1.00     0.00
0.00     0.00     1.00

ATOMIC_POSITIONS
Direct 

Sn
0.0
64
0	0	0	1	1	1
0	0.25	0.25	1	1	1
0.25	0	0.25	1	1	1
0.25	0.25	0	1	1	1
0.375	0.125	0.375	1	1	1
0.125	0.125	0.125	1	1	1
0.125	0.375	0.375	1	1	1
0.375	0.375	0.125	1	1	1
0.5	0	0	1	1	1
0.5	0.25	0.25	1	1	1
0.75	0	0.25	1	1	1
0.75	0.25	0	1	1	1
0.875	0.125	0.375	1	1	1
0.625	0.125	0.125	1	1	1
0.625	0.375	0.375	1	1	1
0.875	0.375	0.125	1	1	1
0	0.5	0	1	1	1
0	0.75	0.25	1	1	1
0.25	0.5	0.25	1	1	1
0.25	0.75	0	1	1	1
0.375	0.625	0.375	1	1	1
0.125	0.625	0.125	1	1	1
0.125	0.875	0.375	1	1	1
0.375	0.875	0.125	1	1	1
0.5	0.5	0	1	1	1
0.5	0.75	0.25	1	1	1
0.75	0.5	0.25	1	1	1
0.75	0.75	0	1	1	1
0.875	0.625	0.375	1	1	1
0.625	0.625	0.125	1	1	1
0.625	0.875	0.375	1	1	1
0.875	0.875	0.125	1	1	1
0	0	0.5	1	1	1
0	0.25	0.75	1	1	1
0.25	0	0.75	1	1	1
0.25	0.25	0.5	1	1	1
0.375	0.125	0.875	1	1	1
0.125	0.125	0.625	1	1	1
0.125	0.375	0.875	1	1	1
0.375	0.375	0.625	1	1	1
0.5	0	0.5	1	1	1
0.5	0.25	0.75	1	1	1
0.75	0	0.75	1	1	1
0.75	0.25	0.5	1	1	1
0.875	0.125	0.875	1	1	1
0.625	0.125	0.625	1	1	1
0.625	0.375	0.875	1	1	1
0.875	0.375	0.625	1	1	1
0	0.5	0.5	1	1	1
0	0.75	0.75	1	1	1
0.25	0.5	0.75	1	1	1
0.25	0.75	0.5	1	1	1
0.375	0.625	0.875	1	1	1
0.125	0.625	0.625	1	1	1
0.125	0.875	0.875	1	1	1
0.375	0.875	0.625	1	1	1
0.5	0.5	0.5	1	1	1
0.5	0.75	0.75	1	1	1
0.75	0.5	0.75	1	1	1
0.75	0.75	0.5	1	1	1
0.875	0.625	0.875	1	1	1
0.625	0.625	0.625	1	1	1
0.625	0.875	0.875	1	1	1
0.875	0.875	0.625	1	1	1
```

The KPT file is:
```
K_POINTS
0
Gamma
1 1 1 0 0 0
```

Run the program, and see results in the output directory. The following files are about MD:
- STRU/STRU_MD_$num: optimized structures in direct coordinate
- Restart_md.dat: output the information of md for restart
- If you want to restart md, you must replace the STRU with STRU/STRU_MD_$num.

MD information can be found in file running_md.log.

```
------------------------------------------------------------------------------------------------
 Energy              Potential           Kinetic             Temperature         Pressure (KBAR)     
 -7.7824215          -7.782469           +4.7502231e-05      +10                 +124.19599          
 ------------------------------------------------------------------------------------------------
```
All energy values are in atomic unit (Ry), and temperature is Kelvin. 

Check these information to confirm whether temperature and conservation is steady.

The MD post-processing tools can be found as follows:
- Candela : https://github.com/MCresearch/Candela
- dpdata : https://github.com/deepmodeling/dpdata

[back to top](#molecular-dynamics)