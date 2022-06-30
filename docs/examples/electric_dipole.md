# Electric field and dipole correction

[back to main page](../../README.md)

## Electric field
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
added to the bare ionic potential. If you want no electric field, `efield_amp`  should be zero. Must be used ONLY in a slab geometry for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE.
- efield_dir : The direction of the electric field or dipole correction is parallel to the reciprocal lattice vector, so the potential is constant in planes defined by FFT grid points, efield_dir = 0, 1 or 2. Used only if `efield_flag` == true.
- efield_pos_max : Position of the maximum of the saw-like potential along crystal axis `efield_dir`, within the  unit cell, 0 < `efield_pos_max` < 1. Used only if `efield_flag` == true.
- efield_pos_dec : Zone in the unit cell where the saw-like potential decreases, 0 < `efield_pos_dec` < 1. Used only if `efield_flag` == true.
- efield_amp : Amplitude of the electric field, in ***Hartree*** a.u.; 1 a.u. = 51.4220632*10^10 V/m. Used only if `efield_flag` == true. The saw-like potential increases with slope `efield_amp`  in the region from (`efield_pos_max`+`efield_pos_dec`-1) to (`efield_pos_max`), then decreases until (`efield_pos_max`+`efield_pos_dec`), in units of the crystal vector `efield_dir`. Important: the change of slope of this potential must be located in the empty region, or else unphysical forces will result.


## Dipole correction
A dipole correction is added to the bare ionic potential.Must be used ONLY in a slab geometry, for surface calculations, with the discontinuity FALLING IN THE EMPTY SPACE.
```
INPUT_PARAMETERS
efield_flag        1
dip_cor_flag       1
efield_dir         2
efield_pos_max     0.5
efield_pos_dec     0.1
efield_amp         0
```
Note: *efield_amp must be zero so that there is no electric field.*

## Electric field and Dipole correction
Both external electric field and dipole correction are added to the bare ionic potential. 
```
INPUT_PARAMETERS
efield_flag        1
dip_cor_flag       1
efield_dir         2
efield_pos_max     0.5
efield_pos_dec     0.1
efield_amp         0.001
```


[back to top](#electric-field-and-dipole-correction)