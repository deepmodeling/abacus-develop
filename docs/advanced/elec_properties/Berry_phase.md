# Berry Phase Calculation

From version 2.0.0, ABACUS is capable of calculating macroscopic polarization of insulators by using the Berry phase method, known as the ["modern theory of polarization"](https://www.sciencedirect.com/science/article/abs/pii/S0022459612003234). To calculate the polarization, you need first to do a self-consistent calculation to get the converged charge density. Then, do a non-self-consistent calculation with berry_phase setting to 1. You need also to specify the direction of the polarization you want to calculate. An example is given in the directory [examples/berryphase/lcao_PbTiO3](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/berryphase/lcao_PbTiO3).

To run this example, first do a self-consistent calculation:
```
cp INPUT-scf INPUT
cp KPT-scf KPT
mpirun -np 4 abacus
```
Then run a non-self-consistent berry-phase calculation:
```
cp INPUT-nscf-c INPUT
cp KPT-nscf-c KPT
mpirun -np 4 abacus
```

In this example, we calculate the electric polarization along c axis for PbTiO~3~, and below are the INPUT file (nscf) and KPT file (nscf):

```
INPUT_PARAMETERS
pseudo_dir      ../../../tests/PP_ORB  //the path to locate the pesudopotential files
orbital_dir     ../../../tests/PP_ORB  //the path to locate the numerical orbital files
ntype         3
ecutwfc       50 // Ry
symmetry      0 // turn off symmetry
calculation   nscf // non-self-consistent calculation
basis_type    lcao // atomic basis
init_chg  file // read charge from files
berry_phase   1 // calculate Berry phase
gdir          3 // calculate polarization along c axis
```

Note: You need to turn off the symmetry when do Berry phase calculations. Currently, ABACUS support Berry phase calculation with nspin=1 and nspin=2. The Berry phase can be calculated in both pw and lcao bases.
- [berry_phase](../input_files/input-main.md#berry_phase) : 1, calculate berry phase; 0, no calculate berry phase.
- [gdir](../input_files/input-main.md#gdir) : 1, 2, 3, the lattice vector direction of the polarization you want to calculate.

The KPT file need to be modified according to gdir in the INPUT file. Generally, you need denser k points along this direction. For example, in the following KPT file, 4 k-points are taken along the a and b axes, and 8 k-points are taken along the c-axis. You should check the convergence of the k points when calculating the polarization.

```
K_POINTS
0
Gamma
4 4 8 0 0 0
```
The results of the berry phase calculation are written in the “running_nscf.log” in the OUT folder. You may search for these results by searching for keywords “POLARIZATION CALCULATION”.

The results are shown as follows:
```
 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 |                                                                    |
 | POLARIZATION CALCULATION:                                          |
 |                  Modern Theory of Polarization                     |
 | calculate the Macroscopic polarization of a crystalline insulator  |
 | by using Berry Phase method.                                       |
 |                                                                    |
 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 VALUES OF POLARIZATION

  The Ionic Phase:   -0.10600
 Electronic Phase:    0.92508

 The calculated polarization direction is in R3 direction

 P =    7.4095194  (mod   18.0922373)  (   0.0000000,   0.0000000,   7.4095194) (e/Omega).bohr

 P =    0.0155792  (mod    0.0380407)  (   0.0000000,   0.0000000,   0.0155792) e/bohr^2

 P =    0.8906925  (mod    2.1748536)  (   0.0000000,   0.0000000,   0.8906925) C/m^2
```

The electric polarization **P** is multivalued, which modulo a quantum e**R**/V~cell~. Note: the values in parentheses are the components of the **P** along the c axis in the x, y, z Cartesian coordinates when set gdir = 3 in INPUT file.