# Initializing SCF
Good initializing would abate the number of iteration steps in SCF.
Charge density should be initialed for constructing the initial hamiltonian operator. 

In PW basis, wavefunction should be initialized for iterate diagonalization method.
In LCAO basis, wavefunction can be read to calculate initial charge density. The wavefunction itself does not have to be initialized.

## Charge Density
`init_chg` is used for choosing the method of charge density initialization.
 - `atomic` : initial charge density by atomic charge density from pseudopotential file under keyword `PP_RHOATOM` 
 - `file` : initial charge density from files produced by previous calculations with [`out_chg 1`](../elec_properties/charge.md).

## Wave function
`init_wfc` is used for choosing the method of wavefunction coefficient initialization.

When `basis_type=pw`, setting of `random` and `atomic` are supported.
Atomic wave function is read from pseudopotential file under keyword `PP_PSWFC`, if setting is `atomic` and number of band of atomic wavefunction less than `nbands` in INPUT file, the extra bands will be initialed by random.

When `basis_type=lcao`, we further support reading of initial wavefunction by setting `init_wfc` to `file`.
In LCAO code, wave function is used to initialize density matrix and real-space charge density.
For such purpose, a file containing wavefunction must be prepared. Such files can be generated from previous calculations with [`out_wfc_lcao 1`](../elec_properties/wfc.md).
