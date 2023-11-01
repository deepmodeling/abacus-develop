# Converging SCF

As in any non-linear systems, numerical instabilities during SCF iterations may lead to nonconvergence. In ABACUS, we offer the following options to facilitate SCF convergence.

## Charge Mixing

By mixing the electron density with that obtained from previous steps, numerical instabilities can be ameliorated. ABACUS offers several mixing schemes, and users may make a selection by adjusting the [mixing_type](../input_files/input-main.md#mixing_type) keyword in INPUT file.

For each of the mixing types, we also provide variables for controlling relevant parameters, including `mixing_beta`, `mixing_ndim`, and `mixing_gg0`.

The default choice is `broyden`, which should work fine in most cases. If convergence issue arises in metallic systems, inclusion of Kerker preconditioning may be helpful, which can be achieved by setting [mixing_gg0](../input_files/input-main.md#mixing_gg0) to be a positive number. For the default broyden method, a choice of 1.5 might be a good start.

A large `mixing_beta` means a larger change in electron density for each SCF step. For well-behaved systems, a larger `mixing_beta` leads to faster convergence. However, for some difficult cases, a smaller `mixing_beta` is preferred to avoid numerical instabilities.

An example showcasing different charge mixing methods can be found in our [repository](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/charge_mixing/pw_Al). Four INPUT files are provided, with description given in README.

## Smearing

Thermal smearing is an efficient tool for accelerating SCF convergence by allowing fractional occupation of molecular orbitals near the band edge. It is important for metallic systems.

In ABACUS, we provide a few smearing methods, which can be controlled using the keyword [smearing_method](../input_files/input-main.md#smearing_method). We also provide keyword `smearing_sigma` or `smearing_sigma_temp` to control the energy range of smearing. A larger value of smearing sigma leads to a more diffused occupation curve.

> Note : The two keywords `smearing_sigma` and `smearing_sigma_temp` should not be used concurrently.

We provide an example showing the importance of smearing in our [repository](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/smearing/lcao_fe). Two INPUT fiels rae provided, with description given in README.