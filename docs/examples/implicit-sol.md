# Implicit solvation model

[back to main page](../../README.md)

Solid-liquid interfaces are ubiquitous in nature and frequently encountered and employed in materials simulation. The solvation effect should be taken into account in accurate first-principles calculations of such systems.  
Implicit solvation model is a well-developed method to deal with solvation effects, which has been widely used in finite and periodic systems. This approach treats the solvent as a continuous medium instead of individual “explicit” solvent molecules, which means that the solute embedded in an implicit solvent and the average over the solvent degrees of freedom becomes implicit in the properties of the solvent bath.

## Input
```
INPUT_PARAMETERS
imp_sol                 1
eb_k                    80
tau                     0.000010798
sigma_k                 0.6
nc_k                    0.00037
```
- imp_sol  

    If set to 1, an implicit solvation correction is considered. 0：vacuum calculation(default).
- eb_k  
    
    The relative permittivity of the bulk solvent, 80 for water. Used only if `imp_sol` == true.
- tau 

    The effective surface tension parameter, which describes the cavitation, the dispersion, and the repulsion interaction between the solute and the solvent that are not captured by the electrostatic terms.
    We use the values of `tau`, `sigma_k`, `nc_k` that were obtained by a fit of the model to experimental solvation energies for molecules in water. tau = 0.525 $meV/Å^{2}$ = 1.0798e-05 $Ry/Bohr^{2}$.
- sigma_k 
    
    We assume a diffuse cavity that is implicitly determined by the electronic structure of the solute. 
    `sigma_k` is the parameter that describes the width of the diffuse cavity. The specific value is sigma_k = 0.6.
- nc_k
    
    `nc_k` determines at what value of the electron density the dielectric cavity forms. 
    The specific value is nc_k = 0.0025 $Å^{-3}$ = 0.00037 $Bohr^{-3}$.

## Output
In this example, we calculate the implicit solvation correction for H2O.
The results of the energy calculation are written in the “running_nscf.log” in the OUT folder.
```
       Energy                       Rydberg                            eV
   E_KohnSham                -34.3200995971                -466.948910448
     E_Harris                -34.2973698556                -466.639656449
       E_band                -7.66026117767                -104.223200184
   E_one_elec                -56.9853883251                -775.325983964
    E_Hartree                +30.0541108968                +408.907156521
         E_xc                -8.32727420734                -113.298378028
      E_Ewald               +0.961180728747                +13.0775347188
      E_demet                            +0                            +0
      E_descf                            +0                            +0
     E_efield                            +0                            +0
        E_exx                            +0                            +0
     E_sol_el              -0.0250553663339               -0.340895747619
    E_sol_cav             +0.00232667606131               +0.031656051834
      E_Fermi               -0.499934383866                 -6.8019562467

```
- E_sol_el: Electrostatic contribution to the solvation energy.
- E_sol_cav: Cavitation and dispersion contributions to the solvation energy.
Both `E_sol_el` and `E_sol_cav` corrections are included in `E_KohnSham`. 



[back to top](#implicit-solvation-model)