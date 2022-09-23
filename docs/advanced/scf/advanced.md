# SCF in Complex Environments

## Implicit Solvation Model

Solid-liquid interfaces are ubiquitous in nature and frequently encountered and employed in materials simulation. The solvation effect should be taken into account in first-principles calculations of such systems so as to obtain accurate predictions.  
Implicit solvation model is a well-developed method to deal with solvation effects, which has been widely used in finite and periodic systems. This approach treats the solvent as a continuous medium instead of individual “explicit” solvent molecules, which means that the solute embedded in an implicit solvent, and the average over the solvent degrees of freedom becomes implicit in the properties of the solvent bath. Compared to the “explicit” method, such implicit solvation model can provide qualitatively correct results with much less computational cost, which is particularly suited for large and complex systems. 

Input parameters that control the implicit solvation model are listed as follows:

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

    If set to 1, an implicit solvation correction is considered; 0 is not applying the solvation model correction (default).
- eb_k  
    
    The relative permittivity of the bulk solvent, e.g., 80 for water. Used only if `imp_sol` == 1.
- tau 

    The effective surface tension parameter, which describes the cavitation, the dispersion, and the repulsion interaction between the solute and the solvent that are not captured by the electrostatic terms.
    We use the values of `tau`, `sigma_k`, `nc_k` that were obtained by a fit of the model to experimental solvation energies for molecules in water. tau = 0.525 $meV/Å^{2}$ = 1.0798e-05 $Ry/Bohr^{2}$.
- sigma_k 
    
    We assume a diffuse cavity that is implicitly determined by the electronic structure of the solute. 
    `sigma_k` is the parameter that describes the width of the diffuse cavity. The specific value is sigma_k = 0.6.
- nc_k
    
    `nc_k` determines at what value of the electron density the dielectric cavity forms. 
    The specific value is nc_k = 0.0025 $Å^{-3}$ = 0.00037 $Bohr^{-3}$.

## External Electric Field and Dipole Correction

### Electric field
A saw-like potential simulating an electric field
is added to the bare ionic potential. 
```
INPUT_PARAMETERS
efield_flag        1
dip_cor_flag       0
efield_dir         2
efield_pos_max     0.5
efield_pos_dec     0.1
efield_amp         0.001
```
- efield_flag : If set to true, a saw-like potential simulating an electric field
is added to the bare ionic potential. 
- dip_cor_flag : If `dip_cor_flag` == true and `efield_flag` == true,  a dipole correction is also
added to the bare ionic potential (see below for details). 
- efield_dir : The direction of the electric field or dipole correction is parallel to the reciprocal lattice vector, so the potential is constant in planes defined by FFT grid points, efield_dir = 0, 1 or 2. Used only if `efield_flag` == true.
- efield_pos_max : Position of the maximum of the saw-like potential along crystal axis `efield_dir`, within the  unit cell, 0 < `efield_pos_max` < 1. Used only if `efield_flag` == true.
- efield_pos_dec : Zone in the unit cell where the saw-like potential decreases, 0 < `efield_pos_dec` < 1. Used only if `efield_flag` == true.
- efield_amp : Amplitude of the electric field, in ***Hartree*** a.u.; 1 a.u. = 51.4220632*10^10 V/m. Used only if `efield_flag` == true. The saw-like potential increases with slope `efield_amp`  in the region from (`efield_pos_max`+`efield_pos_dec`-1) to (`efield_pos_max`), then decreases until (`efield_pos_max`+`efield_pos_dec`), in units of the crystal vector `efield_dir`. Important: the change of slope of this potential must be located in the empty region, or else unphysical forces will result.


### Dipole correction
A dipole correction is added to the bare ionic potential, which can compensate for the artificial dipole field within the context of a periodic supercell calculation. This correction must be used ONLY in a slab geometry, for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE. 
```
INPUT_PARAMETERS
efield_flag        1
dip_cor_flag       1
efield_dir         2
efield_pos_max     0.5
efield_pos_dec     0.1
efield_amp         0
```
Note: *efield_amp must be zero if no electric field is applied.*

### Electric field and Dipole correction
The external electric field and dipole correction can be added to the bare ionic potential. 
```
INPUT_PARAMETERS
efield_flag        1
dip_cor_flag       1
efield_dir         2
efield_pos_max     0.5
efield_pos_dec     0.1
efield_amp         0.001
```

## Compensating Charge

Modeling a constant-potential electronchemcial surface reaction requires adjustment of electron numbers in a simulation cell. At the mean time, we need to maintain the supercell's neutrality due to the periodic boundary condition. A distribution of compensating charge thus needs to be implemented in the vacuum region of surface models when extra electrons are added/extracted from the system.




## Van-der-Waals Correction
