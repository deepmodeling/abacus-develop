# Native SCC-DFTB in ABACUS

The optional `dftbnative` solver evaluates a periodic SCC-DFTB2 or SCC-DFTB3 model directly in ABACUS. It reads Slater–Koster parameter files and does not launch or link to DFTB+. `basis_type=dftb` selects the SKF minimal basis; UPF pseudopotentials and ABACUS `.orb` files are not used by this solver. The species labels and geometry still come from `STRU`.

At present, use `calculation=scf`; the SCC solve and optional frozen-potential band path are performed in that run. `nscf`, relaxation, and molecular-dynamics driver modes are rejected or unsupported because native DFTB restart and force/stress interfaces are not implemented.

## Selecting the solver

Build ABACUS with `-DENABLE_DFTB_NATIVE=ON`, then use:

```text
INPUT_PARAMETERS
suffix dftb_native
calculation scf
basis_type dftb
esolver_type dftbnative
nspin 1
```

The input file named by `dftb_native_input` (default `dftb_native.in`) selects SKF data and DFTB-specific controls. For example:

```text
skf_dir /path/to/3ob-3-1/skfiles
hubbard_deriv C -0.1492
hubbard_deriv N -0.1535
temperature_kelvin 20
scc_tolerance 1e-6
max_scc_iterations 200
mixing_parameter 0.2
mixing_method broyden
mixing_history 6
broyden_inverse_jacobi_weight 0.01
broyden_minimal_weight 1.0
broyden_maximal_weight 1.0e5
broyden_weight_factor 1.0e-2
output_precision 12
third_order yes
band_path_file dftb_band_path.in
```

`hubbard_deriv` is required for every species when `third_order yes` is selected. Supported mixers are damped linear, charge-residual Pulay/DIIS, and Johnson modified-Broyden mixing. Broyden uses secant differences of the SCC charge residual and charge updates, with a regularized small history-space solve; its default controls match DFTB+ (`InverseJacobiWeight=0.01`, `MinimalWeight=1`, `MaximalWeight=1e5`, `WeightFactor=1e-2`). `mixing_history` bounds the stored Pulay residuals or Broyden secant vectors and must be between 2 and 20. If an accelerated solve is singular, ill-conditioned, or non-finite, the iteration falls back to linear damping. Each update is projected to preserve total charge. `output_precision` controls significant digits in native DFTB output files, including `dftb.log`; ABACUS driver logs retain the standard `out_ndigits` formatting.

For Broyden, let `r_n = q_out - q_n` be the SCC charge residual. The mixer stores normalized residual differences and corresponding charge-update differences, solves a regularized weighted secant system in the small history space, and applies the resulting multisecant correction to `q_n + beta r_n`. This is distinct from Pulay residual minimization, although both use a bounded history to accelerate the fixed-point iteration.

## K-point and band inputs

The SCC integration KPT parser accepts ABACUS `K_POINTS`, `KPOINTS`, or `K` headers, Gamma/Monkhorst–Pack automatic meshes, and weighted explicit lists in `Direct`/`D` or `Cartesian`/`C` coordinates. Explicit weights are normalized internally to sum to one, with spin degeneracy carried by occupations (0–2 electrons per state). Cartesian coordinates are converted through the cell vectors.

ABACUS `Line`, `Line_Direct`, and `Line_Cartesian` KPT modes describe a path and do not define Brillouin-zone integration weights. They are therefore rejected for SCC. Supply a mesh or weighted integration list in `KPT` and a separate `band_path_file` for the frozen-SCC band solve. The band-path file starts with `<number_of_segments> <intervals_per_segment>`, followed by connected segments of the form `Gamma 0 0 0 M 0.5 0 0`.

## Method and current scope

The solver builds non-orthogonal Bloch Hamiltonian and overlap matrices from legacy-format two-centre SKF tables with an `s+p` basis (four orbitals per atom), then solves the generalized Hermitian eigenproblem. SCC charges are Mulliken electron populations relative to the SKF neutral valence. Their electrostatic potential is formed from the periodic gamma matrix; DFTB3 adds the implemented charge-dependent third-order potential. Finite-temperature Fermi filling, SKF repulsive energies, 3D Ewald electrostatics, and a post-SCC frozen-potential band-path eigensolve are included.

This makes the implemented DFTB2/DFTB3 charge functional self-consistent within the supported model and parameter format; it is not a full implementation of all DFTB variants or DFTB+ capabilities. Current restrictions are: periodic three-dimensional cells, legacy SKF format, an `s+p` basis, spin-degenerate (`nspin=1`) calculations, monopole SCC charges, and dense diagonalization. Spin-polarized DFTB, shell/multipole charges, modern SKF variants, dispersion add-ons, DOS/PDOS, restart files, analytic forces/stress, structural relaxation, and molecular dynamics are not implemented. Do not request a geometry optimization: force and stress callbacks stop with an explicit unsupported-feature error. Slab electrostatics currently use 3D periodic Ewald sums, so vacuum and image-interaction convergence must be checked for 2D systems.

## Outputs

All files are written under `OUT.{suffix}/`:

- `running_scf.log`: ABACUS driver log and timing summary. DFTB-specific setup, SCC iterations, energies, populations, and output-file details are recorded in `dftb.log`.
- `dftb.log`: DFTB+-style run report with the selected SKF file pairs, species parameters, cell and atom coordinates, weighted SCC k points, SCC iteration table (`E_electronic`, energy change, charge error), final energy decomposition, and atom-resolved Mulliken population/excess charges. The file is flushed each iteration, so it retains the iteration history if SCC fails to converge.
- `eig_occ.txt`: final SCC k-point eigenvalues and occupations.
- `mulliken.txt`: final atom populations and electron-excess charges.
- `band.txt`: eigenvalues along the configured frozen-SCC path, when `band_path_file` is set.

The final energy is reported as the band free energy corrected for SCC-potential double counting, plus the quadratic SCC term, third-order term, and SKF repulsive energy. A converged SCC solution and an accurate band path do not by themselves establish that the parameter set is transferable to other materials; compare energies, charges, and bands against a DFTB+ reference with the same SKF set, Hubbard derivatives, cell, k-point mesh, temperature, and electrostatics.
