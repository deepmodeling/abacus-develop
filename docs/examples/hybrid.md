# Using hybrid functionals in ABACUS

[back to main page](../../README.md)

In ABACUS, we allow LCAO calculations with hybrid functionals. To use hybrid functionals, an input must be provided to the INPUT variable exx_hybrid_type. (Check the section [Exact Exchange](../input-main.md#exact-exchange) in the list of input variables for more information.)

Here we provide an example for Silicon crystal. In this example, HSE functional is used, and it is required that ABACUS is compiled with LIBXC in this case. We provide a short [introduction](../install.md#link-libxc) on compiling ABACUS with LIBXC.

The input files (INPUT, STRU, KPT, pseudopotential and atomic orbitals) are all provided in the directory ${ABACUS_DIR}/examples/hse-Si-example/.

In the directory we further provide two files: running_scf.log_ref and log_ref, which contains reference for running_scf.log and standard output from the program, respectively.

In the log_ref file, you will see the repetitive appearance of a piece of warning message:
```
 The angular momentum larger than 4 (g orbitals) may be error about eggbox.
 Check file ./module_orbital/ORB_atomic_lm.cpp line 272
```
This is normal and it will not affect the results of calculation.

The calculation uses 8\*8*8 k-point mesh. In the reference calculation, the top of valence band is located at the first k point (i.e. gamma point), at 6.11586 eV; the bottom of conduction band is located at the 217-th and 361-th k points, at 7.45955 eV. The calculated fundamental band gap is thus 1.34 eV.

[back to top](#Using-hybrid-functionals-in-ABACUS)