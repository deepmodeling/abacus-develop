# Hefei-NAMD

[Hefei-NAMD](https://github.com/QijingZheng/Hefei-NAMD) Non-adiabatic molecular dynamics applies surface hopping to incorporate quantum mechanical effects into molecular dynamics simulations. Surface hopping partially incorporates the non-adiabatic effects by including excited adiabatic surfaces in the calculations, and allowing for ‘hops’ between these surfaces.  

Detailed instructions on installing and running Hefei-NAMD can be found on its official [website](http://staff.ustc.edu.cn/~zqj/posts/Hefei-NAMD-Training/). 

ABACUS provides results of molecular dynamics simulations for Hefei-NAMD to do non-adiabatic molecular dynamics simulations. 

The steps are as follows :
1. Add output parameters in INPUT when running MD using ABACUS .
```
out_wfc_lcao    1
out_mat_hs      1  
```
Then we obtain output files of hamiltonian matrix, overlap matrix, and wavefunction to do NAMD simulation. 

2. Clone Hefei-NAMD codes optimized for ABACUS from [website](https://github.com/vtzf/abacus-namd).

3. Modify parameters in `Args.py` including directory of ABACUS output files and NAMD parameters. We can see detailed explanation for all parameters in `Args.py`. 

4. Run `NAC.py` to prepare related files for NAMD simulations. 
```
sbatch sub_nac
```

5. Run `SurfHop.py` to perform NAMD simulations. 
```
sbatch sub_sh
```
And results are under directory namddir in `Args.py`.
