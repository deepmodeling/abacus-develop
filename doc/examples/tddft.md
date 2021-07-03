# Real-time time dependent density functional theory

[back to main page](../../README.md)

Real-time time dependent density functional theory (rt-TDDFT) approaches, directly provide
time domain evolution of electronic wave functions together with ionic movements, presenting a versatile way for the real time tracking of ultrafast dynamics and phenomena either in perturbative or non-perturbative regime.

To run Real-time TDDFT (Time-Dependent DFT) simulations with ABACUS, the KPT file and the STRU file should be set in the same way as in Sec. 4.7. A typical INPUT file for TDDFT simulation is given below, which contains additional keywords that need to be added.

```
INPUT_PARAMETERS
#Parameters (General)
pseudo_dir ./
ntype 1
nbands 12
calculation md

#Parameter (Accuracy)
ecutwfc 50
niter 30
smearing gaussian
sigma 0.02
basis_type lcao
out_charge 1
gamma_only 0
nstep 2000
dr2 1.0e-6
md_mdtype 0
md_dt 0.01
md_tfirst 30
md_mdoutpath mdtest
md_rstmd 0
md_fixtemperature 20000
md_nresn 3
md_nyosh 3
md_qmass 0
md_tlast 30
md_dumpmdfred 1
md_domsd 0
md_domsdatom 0
md_outputstressperiod 0
tddft 1
val_elec_01 5
vext 1
vext_dire 2
```

Note: *The TDDFT simulation is based on molecular dynamics in ABACUS. Accomplished with these MD parameters, TDDFT needs some new key words. The TDDFT simulation can be calculated only in lcao bases.*
- tddft : 1, calculate TDDFT; 0, no calculate TDDFT.
- vext : 1, a laser material interaction is included; 0, no extern laser field.
- vext_dire : 1, 2, 3, the direction of the extern laser field added x, y, z.
The results of the TDDFT calculation are dependent on what you want. Dynamics, High-order harmonic generation and occupation could be achieved when you do TDDFT calculation. Users can also do some post-processing according to the basic TDDFT calculation.

[back to top](#real-time-time-dependent-density-functional-theory)