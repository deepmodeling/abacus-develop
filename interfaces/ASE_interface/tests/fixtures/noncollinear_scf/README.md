# Noncollinear SCF handoff fixtures

Unmodified ABACUS output from the local two-engine acceptance on 2026-10-18:
`abacuslite-dual-engine-build-20261018/report-runtime-{develop,lts}/raw/magnetic/`.
The engine-side `magnetic.py` case used two Fe atoms, `nspin=4`, SOC,
noncollinear moments, a 2x2x2 k-grid and one SCF iteration. These are interface
fixtures, not converged scientific reference calculations.

- `latest`: develop `73cd094361f77574f9b28e825913e00f1d2c742a`
  (`v3.11.0-beta10`), run `abacus_ase_nc_dtm5ksbx`.
- `legacy`: LTS `725bf475a7159db4c1890a1f55160393e14b6b13`
  (`v3.10.1`), run `abacus_ase_nc_r9hqihpl`.

`STRU` records the input; the running log and (latest only) `eig_occ.txt`
are original engine outputs. The public-boundary tests replay these outputs
at the launcher boundary, retaining real input writing and result parsing.
They do not patch ASE's global magnetic-property definitions.
