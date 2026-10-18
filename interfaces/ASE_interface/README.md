# ABACUS ASE interface

[`abacuslite`](https://github.com/QuantumMisaka/abacuslite) is the
independently installed ASE calculator plugin for ABACUS. Install it from
PyPI:

```bash
python -m pip install abacuslite
```

The plugin owns the Python calculator, ABACUS input/output readers, examples
and user documentation. This directory contains only the engine-side
compatibility harness used by ABACUS CI; it does not ship a second copy of the
plugin.

`abacuslite` is not vendored and is not a git submodule. Downloading the ABACUS
source archive therefore keeps its engine tree complete; the plugin is resolved
as a normal Python dependency.

## What this harness checks

The harness builds one explicit ABACUS executable and runs the installed
`abacuslite` package against it. Required checks cover SCF, ASE relaxation and
MD, native MD, band structure, magnetic moments and fixed-cell socket handoff.
The develop lane uses the latest output format; the same runner accepts a
legacy-format executable for LTS compatibility.

Every case runs in a separate Python process and preserves its raw ABACUS
output. The runner verifies that imports come from the installed distribution,
that runtime files match the supplied wheel, and that no required case was
skipped or omitted.

## Run locally

Build ABACUS, install a wheel in the Python environment that will run the
harness, and pass both explicitly:

```bash
python -m pip install /path/to/abacuslite-*.whl

python interfaces/ASE_interface/tests/run_engine_compatibility.py \
  --engine /path/to/abacus \
  --mpi-command 'mpirun -np 2' \
  --assets /path/to/abacus-develop/tests/PP_ORB \
  --format latest \
  --wheel /path/to/abacuslite-*.whl \
  --report-dir /tmp/abacus-ase-engine-tests \
  --engine-ref "$(git rev-parse HEAD)"
```

Use `--format legacy` with the supported 3.10.x executable. The LTS engine has
no socket driver; that is a declared capability difference, not a skipped
required case.

For the offline harness contracts:

```bash
cd interfaces/ASE_interface/tests
python -m unittest test_external_provenance test_harness_facts \
  test_native_md_fixture test_portable_runner -v
```

## CI and cross-repository changes

The ASE workflow installs the published, tested `abacuslite` wheel (currently
1.1.0) by default. A manual `workflow_dispatch` run can build a reviewed
`abacuslite` ref for a coordinated change. The `abacuslite` repository
separately builds pinned ABACUS develop/LTS snapshots and runs the same engine
harness on a schedule.

ABACUS owns the engine, output formats and socket wire protocol. `abacuslite`
owns the Python adapter and parser. Parser or calculator fixes belong in
`abacuslite`; engine format changes belong here and require a linked adapter
change.

The socket metadata schema is `abacus.socket.properties.v1`. Additive fields
may be added within v1. A breaking change requires a new schema version and a
transition period in both repositories. See
[the socket contract](https://github.com/QuantumMisaka/abacuslite/blob/main/docs/socket-contract.md)
and [ABACUS ASE documentation](../../docs/advanced/interface/ase.md).

## License

The applicable license terms are provided in the repository
[LICENSE](../../LICENSE).
