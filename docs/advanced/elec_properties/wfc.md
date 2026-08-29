# Extracting Wave Functions

ABACUS is able to output electron wave functions in both PW and LCAO basis calculations. One can find the examples in [examples/11_wfc](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/11_wfc).

## Wave Function in G-Space

For `basis_type=pw` and `esolver_type=ksdft`, [`out_wfc_pw`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-pw) controls the output of plane-wave Kohn-Sham coefficients:

* `0`: Do not write wave-function coefficients.
* `1`: Write text files with the `.txt` suffix.
* `2`: Write binary files with the `.dat` suffix.

The files are stored in `OUT.${suffix}/`. Their pattern is `wfk{k}[s{spin}][g{geometry step}][e{electronic iteration}]_pw.txt` for `out_wfc_pw=1` and `wfk{k}[s{spin}][g{geometry step}][e{electronic iteration}]_pw.dat` for `out_wfc_pw=2`. The `s*` label is omitted for `nspin=1`, is `s1` or `s2` for `nspin=2`, and is `s4` for `nspin=4`. All PW files include the `k*` label, including Gamma-only calculations.

With `out_freq_ion=0`, files are written only when the electronic calculation converges or reaches `scf_nmax`, and the names contain neither `g*` nor `e*`. During structural relaxation or molecular dynamics, each later ionic step overwrites the same files. With `out_freq_ion>0`, output is restricted to the ionic steps selected by `out_freq_ion` and occurs at multiples of `out_freq_elec`, at convergence, or at `scf_nmax`; both `g*` and `e*` are included in the file names. A static `calculation=scf` or `calculation=nscf` run also receives `g1e*` indices when `out_freq_ion>0`.

The normal [`init_wfc=file`](../scf/initialization.md#wave-function) path reads only unindexed binary `wf*_pw.dat` files from `read_file_dir`. Generate directly reusable files with `out_wfc_pw=2` and normally `out_freq_ion=0`. Text `wf*_pw.txt` files and files containing `g*` or `e*` indices are not matched automatically.

For `basis_type=lcao`, [`out_wfc_lcao`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-lcao) controls the output of numerical-atomic-orbital coefficients: `1` writes text files with the `.txt` suffix and `2` writes binary files with the `.dat` suffix. Multi-k calculations generate `wfs{spin}k{k-point}_nao.txt` or `.dat`, while Gamma-only calculations generate `wfs{spin}_nao.txt` or `.dat`.

To reuse LCAO coefficients, set `init_wfc=file`, point `read_file_dir` to the coefficient files, and set `read_wfc_lcao=1` for text or `read_wfc_lcao=2` for binary input. The selected extension is required; there is no automatic format fallback. Only independent files without geometry-step indices are supported, not files accumulated with `out_app_flag`.

## Wave Function in Real Space

One can also choose to output real-space wave functions with the keyword [`out_wfc_norm`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-norm) or [`out_wfc_re_im`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-re-im).

Notice: When the [`basis_type`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#basis-type) is `lcao`, only `get_wf` [`calculation`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#calculation) is effective. An example is [examples/11_wfc/lcao_ienvelope_Si2](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/11_wfc/lcao_ienvelope_Si2).
