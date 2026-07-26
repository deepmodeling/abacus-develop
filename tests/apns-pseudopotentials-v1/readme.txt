Usage
-----
The pseudopotentials (PSPs) in this directory are
tested for its efficiency and precision on ABACUS-PW
code. The convergence thresholds for efficiency
tests are set to 
- 1e-3 eV/atom for Kohn-Sham energy,
- 0.1 kbar for lattice pressure
- 1e-2 eV for band structure
, the minimal-suggested kinetic energy cutoffs of 
planewave for expanding the wavefunction (ecutwfc) 
are, respectively:

    element  cutoff (Ry)
        Ag      90 
        Al      70 
        As      40 
        Au      80 
        B       90 
        Ba      40 
        Be      90 
        Bi      40 
        Br      50 
        C       90 
        Ca      70 
        Cd      60 
        Cl      60 
        Co      150 
        Cr      80 
        Cs      40 
        Cu      150 
        F       100 
        Fe      90 
        Ga      90 
        Ge      70 
        H       80 
        Hf      70 
        Hg      60 
        I       80 
        In      80 
        Ir      60 
        K       40 
        Kr      30 
        Li      80 
        Mg      60 
        Mn      100 
        Mo      100 
        N       90 
        Na      90 
        Nb      90 
        Ni      90 
        O       90 
        Os      100 
        P       40 
        Pb      40 
        Pd      90 
        Pt      100 
        Rb      30 
        Re      150 
        Rh      80 
        Ru      70 
        S       50 
        Sb      70 
        Sc      60 
        Se      40 
        Si      40 
        Sn      50 
        Sr      30 
        Ta      60 
        Tc      70 
        Te      40 
        Ti      70 
        Tl      40 
        V       90 
        W       100 
        Xe      90 
        Y       80 
        Zn      150 
        Zr      100 
        He      50 
        Ne      70 
        Ar      70 
        La      60 

It is also recommended to test the convergence 
of the practical system you are working on, before
any productive calculation.

Cite
----
You are using the pseudopotentials (PSPs) generated
by ONCVPSP code, you are kindly suggested to cite
the following paper:

Hamann D R. Optimized norm-conserving Vanderbilt 
pseudopotentials[J]. Physical Review B—Condensed 
Matter and Materials Physics, 2013, 88(8): 085117.

Additionally, upon the PSPs practically used, please
also cite the following paper:

SG15 family
-----------
If you use PSPs of the elements,

Au, Ba, Be, Bi, Cd, Co, Cr, Cu, Fe, Ge, 
He, Ir, K,  Kr, La, Mn, Mo, N,  Ni, Pb,
Pd, Pt, Rb, Rh, Ru, Sb, Sc, Se, Sn, Sr,
Ta, Tc, Ti, Tl, W,  Zn, Zr

, please cite the paper of SG15 PSPs:

Schlipf M, Gygi F. Optimization algorithm for the 
generation of ONCV pseudopotentials[J]. Computer 
Physics Communications, 2015, 196: 36-44.

Pseudo-DOJO family
-------------------
Otherwise, please cite the following paper:

Van Setten M J, Giantomassi M, Bousquet E, et al. 
The PseudoDojo: Training and grading a 85 element 
optimized norm-conserving pseudopotential table[J]. 
Computer Physics Communications, 2018, 226: 39-54.

Benchmarks
----------
These PSPs were benchmarked on ABACUS-PW code, 
against the all-electron data of systems:

Unaries:
- Body-centered cubic (bcc)
- Face-centered cubic (fcc)
- Diamond
Compounds:
- X2O
- XO
- X2O3
- XO2
- X2O5
- XO3

from the following references:

Bosoni E, Beal L, Bercx M, et al. How to verify 
the precision of density-functional-theory 
implementations via reproducible and universal 
workflows[J]. Nature Reviews Physics, 2024, 6(1): 45-58.

The benchmark data is open-sourced at AIS-Square:
https://aissquare.com/datasets/detail?pageType=datasets&name=ABACUS-stable-psp-orb-test-data&id=318