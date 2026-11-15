# LTS SCF cell-rounding fixture

Taken from the publicly available fixture at QuantumMisaka/abacuslite commit
39ab1d5595bc8456397a1d211946116cf6b15b0c,
`tests/fixtures/lts_scf_cell_rounding`. It records ABACUS v3.10.1 output;
this backport does not represent a new engine run.

The running log rounds lattice components to 10, 11 and 12 Angstrom, whereas
STRU retains 10.0000165221, 11.000002299 and 11.9999880759 Angstrom.
The maximum difference is 1.65221e-5 Angstrom: a 1e-5 identity tolerance
rejects this valid frame. The regression checks that 1e-4 accepts it without
reading abacus.json. The original JSON file is intentionally omitted.
