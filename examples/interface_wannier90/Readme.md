ABACUS to wannier90 step:

1. You need a wannier90 input file `diamond.win`. then,you should run `wannier90 -pp diamond.win` and you will get `diamond.nnkp` file that is the necessary file for ABACUS.
2. Prepare ABACUS's `INPUT`, `STRU`, `KPT`, pseudopotential file, `diamond.nnkp` file. Especially pay attention to kpt file in nscf calculation.
3. First, you should do scf calculation, just like you always do. This step doesn't need `diamond.nnkp` file. Second, you should do nscf calculation, the kpt file is similar to "begin kpoints ..." in the `diamond.win` file.
4. After you run ABACUS nscf calculation, you will get `diamond.amn`, `diamond.mmn` and `UNK*` file in OUT.* folder.
5. Copy `diamond.amn`, `diamond.mmn`, `UNK*` file to wannier folder, then you can run `wannier90 diamond.win`, and you will get MLWF. There is a important things, you must set `wvfn_formatted = .true.` in `diamond.win`, otherwise the wannier90 code cannot read ABACUS file.
