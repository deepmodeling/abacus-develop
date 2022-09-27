# Initializing SCF
Good initializing would abate the number of iteration steps in SCF.
Charge density should be initialed for constructing the initial hamiltonian operator. 

In PW base, wavefunction should be initialed for iterate diagonalization method.
In LCAO base, wavefunction can be read to calculate initial charge density. The wavefunction itself does not have to be initialized.

## Charge Density
`init_chg` is used for choosing the method of charge density distribution initialization.
 - `atomic` : initial charge density by atomic charge density from pseudopotential file under keyword `PP_RHOATOM` 
 - `file` : initial charge density by files which output by setting [`out_chg 1`](../elec_properties/charge.md)

## Wave function
`init_wfc` is used for choosing the method of wavefunction coefficient initialization.

When `basis_type=pw`, setting of `random` and `atomic` are supported.
Atomic wave function is read from pseudopotential file under keyword `PP_PSWFC`, if setting is `atomic` and number of band of atomic wavefunction less than `nbands` in INPUT file, the extra bands will be initialed by random.

When `basis_type=lcao`, setting of `file` are supported.
In LCAO code, wave function is used for init charge density matrix and charge density distribution.
Files which output by setting [`out_wfc_lcao 1`](../elec_properties/wfc.md) need preparation.
