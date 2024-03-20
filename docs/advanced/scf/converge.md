# Converging SCF

As in any non-linear systems, numerical instabilities during SCF iterations may lead to nonconvergence. In ABACUS, we offer the following options to facilitate SCF convergence.

## Charge Mixing

In ABACUS, KS-DFT is solved by self-consistent field (SCF) iteration method. By mixing the electron density with that obtained from previous steps, numerical instabilities can be ameliorated while also accelerating convergence. ABACUS offers several mixing schemes, and users may make a selection by adjusting the [mixing_type](../input_files/input-main.md#mixing_type) keyword in INPUT file.

For each of the mixing types, we also provide variables for controlling relevant parameters, including `mixing_ndim`, `mixing_type`, `mixing_beta`, `mixing_gg0`, `mixing_beta_mag`, `mixing_gg0_mag`, `mixing_gg0_min`, `mixing_angle`.

`mixing_ndim` is the mixing dimensions in DIIS (broyden or pulay) mixing. Gerenally, a larger `mixing_ndim` leads to a better convergence. the default choice `mixing_ndim=8` should work fine in most cases. For `mixing_type`, the default choice is `broyden`, which is slightly better than `Pulay` typically. Besides that, a large `mixing_beta` means a larger change in electron density for each SCF step. For well-behaved systems, a larger `mixing_beta` leads to faster convergence. However, for some difficult cases, a smaller `mixing_beta` is preferred to avoid numerical instabilities.

For non-spin-polarized calculations, the default choices usually achieve convergence. If convergence issue arises in metallic systems, you can try different value of Kerker preconditioning [mixing_gg0](../input_files/input-main.md#mixing_gg0) and [mixing_gg0_min](../input_files/input-main.md#mixing_gg0_min), and try to reduce `mixing_beta`, which is 0.8 defaultly for `nspin=1`.

For magnetic calculations, `mixing_beta_mag` and `mixing_gg0_mag` are activated. Considering collinear calculations, you can rely on the default value for most cases. If convergence issue arises, you can try to reduce `mixing_beta` and `mixing_beta_mag` together. For non-collinear calculations, tradtional broyden usually works, especially for a given magnetic configuration. If one is not interested in the energies of a given magnetic configuration but wants to determine the ground state by relaxing the magnetic moments’ directions, the standard Broyden mixing algorithm sometimes fails to find the correct magnetic configuration. If so, we can set [mixing_angle=1.0](../input_files/input-main.md#mixing_angle), which is a promising mixing method proposed by J. Phys. Soc. Jpn. 82 (2013) 114706.

An example showcasing different charge mixing methods can be found in our [repository](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/charge_mixing/pw_Al). Four INPUT files are provided, with description given in README.

## Smearing

Thermal smearing is an efficient tool for accelerating SCF convergence by allowing fractional occupation of molecular orbitals near the band edge. It is important for metallic systems.

In ABACUS, we provide a few smearing methods, which can be controlled using the keyword [smearing_method](../input_files/input-main.md#smearing_method). We also provide keyword `smearing_sigma` or `smearing_sigma_temp` to control the energy range of smearing. A larger value of smearing sigma leads to a more diffused occupation curve.

> Note : The two keywords `smearing_sigma` and `smearing_sigma_temp` should not be used concurrently.

We provide an example showing the importance of smearing in our [repository](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/smearing/lcao_fe). Two INPUT fiels rae provided, with description given in README.