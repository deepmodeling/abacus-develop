1, modify file general_info

2, try this script for autotest:
./Autotest.sh > check.txt

3, if want to focus on No. xxx example, such as 101_PW_OU_PL:
cd 101_PW_OU_PL
./../run_debug.sh $parameter

# you can choose $parameter among "" (empty), debug or ref

# 'ref': generate result.ref file

NOTES:
****Purpose for this autotest package is:
1,  cover all functions and application scenarios
2,  autotest script to check if the version is correct.
    (reference result calculated by one core and saved in result.ref)
    (user can change NUMBEROFPROCESS value in general_info to test by multi cores)
3,  if there is any bug occurs, rename the bug example path to mark it please.

#######################
##################

"simple" scripts

---

0-99 Structure Tests (Bravis Lattices, symmetry analysis)

---

100-200 PW basic properties

100-110: pseudopotentials
(1) Semiconductor (Si), UPF1.0 (s,p,d,f)
(2) UPF2.0 (s,p,d,f)
(3) dojo (s,p,d)

110-120: smearing methods
(1) Metals
(2)

120-130: charge mixing
(1) Metals, Charge Mixing Methods,
(2) Surfaces, Charge Mixing Methods

130-140: diagonalization methods
(1) CG (number of bands/trial wave functions)
(2) Davidson

140-150: force tests

* test every components of forces
* small displacements .vs. large displacements

150-160: stress tests

160-170: ion relaxatoin

170-180: cell relaxation

180-190: md simulations

---

200-300 LCAO basic properties
300-400 LCAO basic Gamma Only properties

_PW    plain wave bases
_NO    numerical atoms orbitals bases
_NP    expand numerical atoms orbitals in plane-wave basis

_OU    old upf pseudopotential file
_15    SG15 pseudopotential file
_VW    vwr pseudopotential file
_DJ    DOJO pseudopotential file

_FD    smearing method Fermi-dirac
_FX    smearing method: Fixed occupations
_M2    smearing method: mp2
_MP    smearing method: Methfessel-Paxton (MP)
_MV    smearing method: Marzari-Vanderbilt
_SG    smearing method: Gaussian

_SY    turn on symmetry
_CG    cg diagonalization method
_DA    david diagonalization method

_S1    one spin channel
_S2    two spin channels
_S4    four spin channels

_GE    genelpa diagonalization method
_HP    hpseps  diagonalization method
_SC    scalapack diagonalization method

_RE    relax calculation
_CR    cell-relax calculation
_CF    calculate and output force
_CS    calculate and output stress
_MD    molecular dynamics
_TD    TDDFT: time dependent DFT

_OH    output Halmitonian matrix
_OB    output bands file
_OD    output DOS file
_OW    output wave functions
_OC    output charge density
_OK    output kinetic energy density

_GO    gamma_only method
_KP    all K-Points method

_FM    ferromagnetic nspin=2
_AF    anti-ferromagnetic nspin=2 anti initial magnetism

_PU    DFT plus U
_BS    bsse

_PL    mixing_type plain mixing
_KK    mixing_type kerker mixing
_PU    mixing_type pulay mixing
_PK    mixing_type pulay-kerker mixing
_BD    mixing_type broyden mixing

_SO    spin orbit coupling (SOC)
_NB    set NBANDS without default

_XX    EXX
_VD    VDW (both d2 or d3)

_MG    move ions method: cg
_MF    move ions method: FIRE
_MB    move ions method: bfgs

_1O    first-order charge extrapolation
_2O    second-order charge extrapolation

---

900-1000 Orbital free DFT properties

_OF    orbital free density functional theory (OFDFT)
_OP    optimization method used in OFDFT
_KE     kinetic energy functional used in OFDFT
_CO    convergence check
