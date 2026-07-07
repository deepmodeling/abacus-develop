# Optional ATST-Tools Workflow Examples

These examples show how ATST-Tools can use abacuslite as the ABACUS ASE
calculator backend. They are an optional workflow layer on top of abacuslite:
abacuslite itself does not depend on ATST-Tools.

Install the optional workflow package in an environment where this checkout is
available:

```bash
pip install atst-tools
```

Each case is intentionally small and uses pseudopotential and orbital files from
`tests/PP_ORB`. Run commands from the case directory so relative structure and
PP/ORB paths resolve as written.

| Case | Workflow | Purpose |
| --- | --- | --- |
| `relax_si` | `calculation.type: relax` | Basic geometry optimization through `atst run`. |
| `neb_si` | `calculation.type: neb` | Generated Si path with a short two-stage NEB setup. |
| `sella_h2` | `calculation.type: sella` | Standalone single-ended saddle search syntax. |
| `ccqn_h2` | `calculation.type: ccqn` | CCQN with an explicit reactive bond. |
| `md_si` | `calculation.type: md` | ASE-driven MD and ABACUS-native MD templates. |
| `abacus_helper` | `atst abacus prepare/collect` | INPUT/KPT/STRU preparation and conservative output collection. |

General workflow:

```bash
cd interfaces/ASE_interface/examples/atst_tools/relax_si
python make_inputs.py
atst config validate config.yaml --print-normalized
atst run --dry-run config.yaml
atst run config.yaml
```

For SAI GPU runs, keep the same YAML structure but change
`calculator.abacus.parameters.ks_solver` from `genelpa` to `cusolver` if the
loaded ABACUS build requires the GPU solver.

The helper case does not launch ABACUS:

```bash
cd interfaces/ASE_interface/examples/atst_tools/abacus_helper
python make_inputs.py
atst abacus prepare config.yaml --structure inputs/init.extxyz --output-dir prepared_abacus --force
atst abacus collect prepared_abacus --output abacus_results.json
```

`collect` is normally useful after ABACUS has populated a run directory; using it
on `prepared_abacus` only exercises the conservative file-summary path.
