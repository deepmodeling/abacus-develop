# ABACUS JSON Development Guide

## Overview

ABACUS uses [nlohmann-json](https://github.com/nlohmann/json) as the backend for its optional JSON output. The JSON implementation is kept under `source/source_io/module_json`, with `AbacusJson` providing the small interface used to build and write `abacus.json`.

The public alias and mutation interfaces are:

```cpp
using jsonValue = nlohmann::ordered_json;

// Public static members of Json::AbacusJson:
static void set_json(const std::vector<jsonKeyNode>& keys, jsonValue value);
static void append_json(const std::vector<jsonKeyNode>& keys, jsonValue value);
```

`jsonValue` uses `nlohmann::ordered_json` so that object keys are written in insertion order. `jsonKeyNode` accepts either a string key or an integer array index, so paths can mix JSON objects and arrays.

`abacusjson.h` includes only `nlohmann/json_fwd.hpp`. A source file that constructs or operates on `jsonValue` must include `<nlohmann/json.hpp>` itself, inside the `__JSON` guard. Callers of the higher-level functions in `init_info.h` and `output_info.h` do not need the backend header.

## Adding values

### Add or replace an object member

Use `set_json()` to assign a value at a path:

```cpp
Json::AbacusJson::set_json({"general_info", "version"}, version);
```

Missing intermediate named nodes are created as objects. The final value is replaced regardless of its previous type, including when it is an array or an object. For example, setting a complete coordinate array replaces the old coordinates rather than adding another nested array:

```cpp
Json::AbacusJson::set_json({"init", "coordinate"}, coordinates);
```

Replacing a complete object also replaces all of its members; this is not a merge operation.

### Append to an array

Use `append_json()` to append one value to an array:

```cpp
Json::AbacusJson::append_json({"init", "label"}, label);
```

A missing final named member is created as an array. An existing destination must already be an array: appending to a scalar, an object, or `null` is an error rather than an implicit conversion.

For nested arrays, construct the value with `jsonValue::array()`:

```cpp
Json::jsonValue coordinate = Json::jsonValue::array({x, y, z});
Json::AbacusJson::append_json({"init", "coordinate"}, coordinate);
```

The coordinate is appended as **one row**; its elements are not flattened into the destination array. An empty path is a no-op for both `set_json()` and `append_json()`.

### Construct objects and arrays

Use the nlohmann-json initializer syntax through the `Json::jsonValue` alias. There is no need for backend-specific helper macros.

Object example:

```cpp
Json::jsonValue scf = {
    {"energy", energy},
    {"ediff", ediff},
    {"drho", drho},
    {"time", time},
};
```

Array example:

```cpp
Json::jsonValue row = Json::jsonValue::array({x, y, z});
```

Append a completed SCF record with:

```cpp
Json::AbacusJson::append_json({"output", -1, "scf"}, scf);
```

Construct complete sections or arrays locally before storing them where practical. `gen_general_info()` assigns its complete section once. `gen_stru()` constructs each structure field locally, and `gen_init()` does the same for calculation metadata. These two generators share `init` with `add_nkstot()`, so they replace only their own fields through a file-local helper; they must not replace the entire `init` object and discard fields written by another generator.

For a current output record, coordinates, magnetic moments, the cell, forces, and stress are replaced as complete arrays. Repeating the geometry update for the same record therefore does not accumulate extra rows. Only genuinely sequential data, such as `output` records and `scf` iteration records, use `append_json()`.

## Addressing array elements

Integer path components address existing array elements. Non-negative indices count from the beginning, while negative indices count from the end (`-1` is the last element). Indexed traversal never grows an array.

For example, given:

```json
{
    "Json": {
        "key6": {
            "key7": [
                {"a": 1, "new": 2},
                "vasp",
                "abacus"
            ]
        }
    }
}
```

replace `"vasp"` with `"cp2k"` using either its forward index:

```cpp
Json::AbacusJson::set_json({"Json", "key6", "key7", 1}, "cp2k");
```

or the corresponding negative index:

```cpp
Json::AbacusJson::set_json({"Json", "key6", "key7", -2}, "cp2k");
```

When the destination selected by an integer is itself an array, `append_json()` appends to that nested array; it does not replace the selected element. Out-of-range indices and mismatched object/array path components are errors.

The workflow must call `init_output_array_obj()` before filling the corresponding calculation/ionic-step record. `set_json()` and `append_json()` do not create an implicit current output record when traversing `{"output", -1, ...}`. Record initialization remains the responsibility of the existing driver/solver entry points, not the generic path interface.

## Migrating older JSON call sites

The former `add_json(keys, value, is_array)` interface has been removed. Choose the new operation by intent, not just by the old boolean:

- Use `set_json()` for scalar assignments, whole-container replacement, and replacement of an indexed element.
- Use `append_json()` for adding one element to a named or indexed array.

The old interface appended to an existing named array even when `is_array` was `false`, and it replaced an indexed element even when the flag was `true`. Neither implicit behavior is retained by the new operation names.

## Code structure

The JSON implementation is organized as follows:

```text
source/source_io/module_json/
├── abacusjson.cpp/.h   # set/append path handling and file output
├── json_node.h         # object-key / array-index path component
├── general_info.cpp/.h # general_info section
├── init_info.cpp/.h    # comment and init sections
├── output_info.cpp/.h  # output section
├── para_json.cpp/.h    # integration-facing wrappers
└── test/              # focused unit tests
```

JSON support is compiled under `__JSON`, which is enabled by the CMake option `ENABLE_JSON`.

## Guidelines for extending JSON output

When adding JSON output:

1. Keep JSON construction in `source/source_io/module_json` whenever practical, rather than spreading nlohmann-json details into unrelated modules.
2. Pass the data required for output explicitly through function parameters. Do not add new `GlobalV`, `GlobalC`, or `PARAM` accesses merely to obtain a value for JSON output.
3. Prefer existing domain objects or small scalar/reference parameters over introducing new cross-module dependencies.
4. Use `Json::jsonValue` for compound JSON values, `set_json()` for assignment, and `append_json()` for sequence growth.
5. Preserve the existing JSON schema unless the change intentionally modifies the public output format.
6. Add or update focused tests under `source/source_io/module_json/test` for new fields and for array/object behavior.

For example, `output_info` receives the required values as function arguments and adds them to the current output record:

```cpp
void add_output_scf_mag(const double total_mag,
                        const double absolute_mag,
                        const double energy,
                        const double ediff,
                        const double drho,
                        const double time)
{
    AbacusJson::set_json({"output", -1, "total_mag"}, total_mag);
    AbacusJson::set_json({"output", -1, "absolute_mag"}, absolute_mag);
    AbacusJson::append_json({"output", -1, "scf"},
                            {{"energy", energy},
                             {"ediff", ediff},
                             {"drho", drho},
                             {"time", time}});
}
```

This keeps the JSON layer explicit and avoids introducing additional global dependencies into the output path.
