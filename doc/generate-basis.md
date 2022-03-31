# Generating atomic orbital bases

[back to main page](../README.md)

In ABACUS, the atomic orbital bases are generated using a scheme developed in the [paper](https://iopscience.iop.org/article/10.1088/0953-8984/22/44/445501). We provide a script named “generate_orbital.sh” under the directory tools/ to generate the atomic orbitals bases. In order to run this script, an ORBITAL_INPUT file is required.

An example of this ORBITAL_INPUT file can be found in $ABACUS/tools/SIAB/SimulatedAnnealing/example_N:
```
#1.exe_dir
#----------------------------------------------------------------------------
EXE_pw bin/ABACUS.fp_mpi.x
EXE_orbital bin/SIA_s.exe
#----------------------------------------------------------------------------
#( In this part, the direction of the two used exe is provided )
#2.electronic calculatation
#----------------------------------------------------------------------------
targets 07_N # element
ref_bands 5 # reference bands
nbands 8 # num of bands for calculate ( > reference bands)
Ecut 50 # cutoff energy (in Ry)
Rcut 6 # cutoff radius (in a.u.)
Pseudo_dir ./
Pseudo N.LDA.UPF
smearing_sigma 0.01 # energy range for gauss smearing (in Ry)
#----------------------------------------------------------------------------
#( In this part , some parameters of calculating are given )
#3.structure information
#----------------------------------------------------------------------------
Dis 1.0 1.1 1.5 2.0 3.0
#----------------------------------------------------------------------------
#( In this part , it gives us the bond length of the reference system( in
#angstrom) )
#4.orbital calculatation
#----------------------------------------------------------------------------
maxL 2 # the max angular momentum (L<=2)
Level 2 # num of levels to generate orbitals(<=5)
#(num) (the max ang) (num of S) (num of P) (num of D)
level1 1 1 1
level2 2 1 1 1
#----------------------------------------------------------------------------
#( In this part, some information of orbital is given )
#5.Metropolis parameters (in most cases do not need to change)
#----------------------------------------------------------------------------
Start_tem_S 1.0e-4 # start temperature for optimize Spillage(default 1.0e-4)
Start_tem_K 1.0e-2 # start temperature for optimize Kinetic (default 1.0e-2)
Step_S 20 # num of steps for optimize Spillage (default 20)
Step_K 15 # num of steps for optimize Kinetic (default 15)
Delta_kappa 0.01 # delta kappa (default 0.01)
#----------------------------------------------------------------------------
#(In this part, some parameters of Metropolis is given. In most cases, they
#do not need to be changed , only when you run into a situation , that the
#Kinnetic energy is larger than the maximum value allowed , you can enlarge
#the start temperature appropritely, or you can enlarge the delta_kappa, e.g.
#start_tem_k 1.0e-2 to 5.0e-1, delta_kappa 0.01 to 0.02. more steps can make 
#the orbitals better , too)
```

The ORBITAL_INPUT file contains 5 parts :
1. exe_dir

    The paths of two executable files:

    - EXE_pw: executable file of ABACUS
    - EXE_orbital: executable file of orbital generation

    We can get the exe file of orbital generation as below:
    ```
    cd tools/SIAB/1_Source/
    make s
    ```
2. electronic calculation

    Parameters for electronic calculations:

    - targets
    
        the element type of which orbitals are to be generated. Its value has form of ‘element.id_element’, for example 07_N.
    - ref_bands
    
        the number of reference bands for orbital generation. We usually take the number of occupied bands of the system. For the N element, we take its “dimers” as the reference systems, so the number of ref_bands should be 5 (valence electrons of this element)\*2(number of the atoms of the system)/2 (1 band contain 2 electrons)=5. While for Na element, we take trimer as reference systems, and the number of the ref_bands should be 1\*3/2=1.5 for 1.5 is not a integer, here we use 2 for its ref_bands. Most elements use dimer as reference systems, except for Li, Na, K, Ca, which use trimer instead.
    - nbands

        the number of bands to be calculated in electronic calculation. Here, we use the gaussian smearing for the electronic structure calculation, so the value of this parameter can not be smaller than the value of ref_bands
    - Ecut

        the Energy cutoff (in Ry)
    - Rcut

        the radius cutoff of atomic orbital (in a.u)
    - Pseudo_dir

        the path to the directory where the pseudopotential file is.
    - Pseudo

        the file name of pseudopotential
    - smearing_sigma

        the gaussian smearing (in Ry) for scf calculations. The default vaule is 0.01. In case that the scf iterations don’t converage (which could happen, e.g., for transition metal dimers), the user may increase this parameter, say, to 0.05.
3. structure information

    This part gives the bond lengths of the reference systems (dimer or trimer). Generally, the bond lengths are chosen to distribute on both sides of the equilibrium value. For example, for N dimer we use (in Å):
     - Dis 1.0 1.1 1.5 2.0 3.0

    It means we take 5 reference systems (dimer), and the bond lengths are 1.0 1.1 1.5 2.0 3.0 angstrom, respectively. Every element has reference systems with different bond lengths, which could be found in file $ABACUS/tools/SIAB/DIS.txt.
4. orbital generation

    The main parameters for orbital generation

    - maxL
        
       the max angular momentum for orbitals to be generated
    -  level
       
       number of levels to generate orbitals. In the main part of this section, level1, level2... provide the information of each layer, the max angular momentum and the number of s, p, d orbitals.
    
    For example, if we want to use 2 steps to generate DZP basis for N, we can set this part like this:
    ```
    maxL 2 # the max angular momentum (L<=2)
    level 2 # num of levels to generate orbitals(<=5)
    (num) (the max ang) (num of S) (num of P) (num of D)
    level1 1 1 1
    level2 2 1 1 1
    ```
    Because of the property of symmetry, taking dimer or trimer as reference systems can not generate the f orbitals very well. And test results show that f orbital has little effect on improving the results. It means we generate one s orbital and one p orbital in first step (level1), and generate one s, p, d orbital in the second step (level2).
5. Metropolis parameters

    The main parameters for Metropolis optimization.

    - Start_tem_S

        start temperature for spillage optimization
    - Start_tem_K

        start temperature for kinetic energy optimization
    - Step_S

        number of steps for spillage optimization
    - Step_K

        number of steps for kinetic energy optimization
    - Delta_kappa

        the accept rise of spillage when optimizing the kinetic energy
        
After preparing the ORBITAL_INPUT file, one just needs to run the script "$PATH_TO/generate_orbital.sh ORBITAL_INPUT" and wait for the results. The results will be written into several output files under the directory $element.id_element/$Rcut/, for example 07_N/6/.

Some output files listed here are useful.
- ORBITAL_RESULTS.txt

    this file shows some important information of the orbital, we can see the spillage of the orbital to judge whether the orbital is good enough for us to use.
- running_1.txt
    the details of generating orbitals
- ORBITAL_PLOTU.dat
    you can open it by using any drawing softwares to visualize the shape of the orbital
- ORBITAL_7U.dat
    the general type is ORBITAL_($element.id)U.dat . This is the orbital file we will use in the calculation. And you can rename it as anything you want. We usually use the covention “element_xc_rcut_ecut_XsY pZd”, e.g., N_lda_6.0au_50Ry_2s2p1d, which tells the key parameters for the basis set construction.

For some elements, you can download the reference ORBITAL_INPUT files and pseudopotentials from our [website](http://abacus.ustc.edu.cn/pseudo.html).

A file README is also given and you can decide the parameters with it as a reference.
In most cases, you just need to modify the parameters in Section 1, 2. Section 4 may be
partially modified if you need higher precision orbitals. The users are not encouraged to change
the settings in sections 5, unless you are very familiar with the code generating algorithms.
