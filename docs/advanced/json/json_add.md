# ABACUS JSON Development Guide

## Overview

ABACUS uses [nlohmann-json](https://github.com/nlohmann/json) for its optional JSON output. The implementation lives in `source/source_io/module_json` and uses `Json::jsonValue`, an alias for `nlohmann::ordered_json`, to retain object-key insertion order.

`AbacusJson` provides access to the shared document and writes it to a file. Its declarations are in namespace `Json`:

```cpp
using jsonValue = nlohmann::ordered_json;

class AbacusJson
{
  public:
    static jsonValue& document();
    static void write_to_json(const std::string& filename);

  private:
    static jsonValue doc;
};
```

Keep the document root an object. Its state remains shared within each process; this change does not introduce independent output contexts or make concurrent writes safe. The mutable accessor is for the schema generators and tests in `module_json`. Other modules should continue to pass data to functions such as `add_output_energy()` instead of directly editing the document.

The old path-component type and generic set/append interface have been removed. Use native object assignment, shallow `update()`, and array `push_back()` inside the schema generators; do not introduce another generic path wrapper.

`abacusjson.h` includes only `nlohmann/json_fwd.hpp`. Source files that construct or manipulate JSON values must include `nlohmann/json.hpp` under `__JSON`. The existing CMake option `ENABLE_JSON` controls this feature. Callers using only the higher-level declarations in `init_info.h` or `output_info.h` do not need the backend header.

## Constructing metadata

`gen_general_info()` owns the whole `general_info` section and assigns it as a complete object:

```cpp
AbacusJson::document()["general_info"] = {
    {"version", version},
    {"commit", commit},
    {"device", param.inp.device},
    {"mpi_num", mpi_num},
    {"omp_num", omp_num},
    {"pseudo_dir", param.inp.pseudo_dir},
    {"orbital_dir", param.inp.orbital_dir},
    {"stru_file", param.globalv.global_in_stru},
    {"kpt_file", param.inp.kpoint_file},
    {"start_time", start_time_str},
    {"end_time", end_time_str}};
```

The `init` section is shared by `gen_stru()`, `gen_init()`, and `add_nkstot()`. The first two construct the fields they own in a local object, then apply a **shallow** update:

```cpp
// Inside init_info.cpp; init_section() is local to this source file.
init_section().update(info);
```

The local helper creates a missing `init` object but rejects an existing non-object, including `null`. The update preserves fields supplied by the other generators and replaces each supplied value as a whole. In particular, per-species maps and coordinate arrays must not retain stale entries or accumulate on repeated generation. Do not assign a newly generated object to the entire `init` section, and do not enable recursive object merging here.

`add_nkstot()` only sets its own field:

```cpp
init_section()["nkstot"] = nkstot;
```

## Output-record lifecycle

The workflow starts each record with `init_output_array_obj()` **before** the corresponding solver writes SCF or other result data. That function alone creates the `output` array and appends the initial record. It rejects an existing `output` value that is not an array; an explicit `null` is not treated as a missing field.

The existing workflow entry points own this initialization:

| Workflow | Record initialization |
| --- | --- |
| SCF/relaxation | `Relax_Driver::iter_info()` starts the record, except for the first `ks-lr` step described below. |
| `ks-lr` | `ESolver_LR::before_all_runners()` starts the record before its embedded KS calculation; the first relaxation-driver step reuses it. |
| UnitCell-backed MD | `Run_MD::md_line()` starts a record at the beginning of each MD iteration when `mdcell.has_backing_unitcell()` is true. |
| Socket/i-PI | `SocketHandlers::handle_posdata()` starts a record before running the solver for the received `POSDATA` frame. |

Do not move record creation into individual field writers, create a second record for the same step, or reset the whole document to start a new step.

The result writers use `current_output()`, a helper local to `output_info.cpp`. It rejects a missing or non-array `output`, an empty array, or a final element that is not an object. It never creates a record as a side effect of writing a result.

For example, inside namespace `Json` in `output_info.cpp`:

```cpp
void add_output_energy(const double energy)
{
    current_output()["energy"] = energy;
}
```

Coordinate, force, stress, magnetic-moment, and cell arrays are built locally and assigned as complete arrays. Repeatedly updating the same record must replace these arrays rather than append rows.

SCF iterations are different: they form a history and must be appended. `add_output_scf_mag()` creates a missing `scf` array, rejects an existing non-array history, and appends one iteration object. Its implementation uses:

```cpp
jsonValue& output = current_output();
output["total_mag"] = total_mag;
output["absolute_mag"] = absolute_mag;
jsonValue& scf = *output.emplace("scf", jsonValue::array()).first;
if (!scf.is_array())
{
    throw std::invalid_argument("JSON SCF history must be an array");
}
scf.push_back({{"energy", energy}, {"ediff", ediff},
               {"drho", drho}, {"time", time}});
```

`ordered_json` may invalidate references to child values when new members are inserted into their parent object. Acquire the `scf` reference after inserting `total_mag` and `absolute_mag`, and do not retain a record reference across appending another `output` record. The same caution applies to references to root sections when new root keys are inserted.

## Serialization and tests

`document()` and `write_to_json()` do not perform MPI rank filtering. The existing `json_output()` wrapper writes `abacus.json` only on rank 0 in MPI builds; callers outside `module_json` should retain the existing integration wrappers.

`write_to_json()` preserves the existing four-space formatting and reports file-open and write/close failures. It serializes the document before opening the destination, so a serialization error does not first truncate the file. Non-finite numbers serialize as JSON `null`; decimal versus scientific float notation is not part of the schema contract.

The tests reset the shared document through `document()` in their fixture; no access-control macro or friend accessor is needed. Focus coverage on ABACUS behavior: generated fields and units, repeated metadata updates, record initialization and SCF accumulation, invalid section types, insertion order, escaping and non-finite values through the real writer, and file errors. Do not replace removed path-walker tests with tests of nlohmann-json's generic container API.

## Code structure

```text
source/source_io/module_json/
├── abacusjson.cpp/.h   # shared document and file output
├── general_info.cpp/.h # general_info section
├── init_info.cpp/.h    # comment and init sections
├── output_info.cpp/.h  # output records and lifecycle checks
├── para_json.cpp/.h    # integration-facing wrappers
└── test/              # focused unit tests
```

`init_section()` and `current_output()` are file-local helpers, not public interfaces for workflow callers.

## Guidelines for extending JSON output

Keep construction in the existing schema generator, pass its required data explicitly, and avoid adding `GlobalV`, `GlobalC`, or `PARAM` access. Preserve field names, value types, units, and order unless a schema change is intentional. Add focused tests for new fields and lifecycle behavior, and update the [JSON output reference](json_para.md) when the public schema changes.

Keep examples and new implementation code compatible with the C++11 baseline. Include complete domain-type definitions in the source or test file that needs them, keep public header dependencies minimal, and do not reintroduce access-control macros for testing.
