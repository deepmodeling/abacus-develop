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

## Public interface and parser ownership

The plugin's [public API contract](https://github.com/QuantumMisaka/abacuslite/blob/main/docs/public-api.md)
defines supported imports, result units/shapes, missing values and compatibility
rules. That is the canonical specification; this directory does not duplicate
it or add a forwarding parser.

Use `Abacus`, `AbacusProfile`, `AbacusSocketIO` and the documented explicit I/O
functions. Tests use explicit curated pseudopotential/orbital filenames and
the socket's public `close()` lifecycle. They do not inspect template reader
flags, use asset-discovery heuristics, manipulate server process handles or
modify ASE's global magnetic-property definitions.

Report parser/calculator defects to abacuslite, engine/output defects to ABACUS,
and launcher/workflow policy defects to the consuming project. Link both sides
when a change affects an output format or socket protocol. The
[maintenance options](https://github.com/QuantumMisaka/abacuslite/blob/main/docs/maintenance-options.md)
are a community discussion proposal; the repository hosting and maintenance
appointments remain undecided.

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

The v2 report records the engine's public version and the requested oracle
format, without claiming to inspect an internal reader selection. The current
lanes accept ABACUS 3.10/LTS and 3.11/develop respectively. Actual parser
compatibility is demonstrated by the real calculations and independent raw-file
assertions. Socket close is included in the whole-case process-group timeout;
timeout/failure remains a failing case with retained raw output.

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

This public-interface follow-up needs the paired abacuslite candidate
`1.2.2.dev0` (or its subsequent reviewed release) for noncollinear result
projection without an ASE registry patch. The published 1.1.0 workflow default
does not include that fix. Validate the candidate pair with `plugin_ref` after
the plugin commit is public; before integrating this follow-up, publish the
reviewed plugin and update the workflow's stable version/hash together.
Do not report the old stable-default CI as acceptance of this candidate.

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
