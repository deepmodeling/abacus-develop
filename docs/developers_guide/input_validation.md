# Shared INPUT validation

INPUT contracts are declared on `Input_Item::domain`. The common engine checks
scalar integer, real, Boolean and string values, finite real numbers, choices,
numerical bounds, and conditional requirements. It never branches on parameter
names. For example:

```cpp
item.domain.set_choices({"1", "2"});
item.domain.set_rule("value>0 and value<=1");
item.domain.add_requirement("noncolin==true or lspinorb==true", "nspin==4");
```

These are independent examples of declarations. Ordinary scalar readers register
their numeric/Boolean/string type automatically. A supplied ordinary string
parameter must have exactly one token, so a bare `stru_file` is rejected with
the expected type and manual link before its reader accesses the value. This
also rejects extra tokens that would otherwise be silently discarded.
The existing `sync_*` registrations
supply typed, read-only value access, including nested MD fields and vectors.
The checker does not infer token cardinality from a prose documentation type,
or format final values and parse them back. Custom input grammars remain in
their readers until their contracts are declared; legal empty strings and
vectors must remain accepted.

Conditions reuse the parsed availability grammar (`and`, `or`, comparisons,
`in`, `contains`). Runtime requirements are explicit declarations: descriptive
availability alone does not imply that a disabled switch is an error.
References and operand types are checked when registrations are assembled.

## Preparation order

1. Check all explicit domains before readers or reset hooks can hide invalid
   input. Complete tokens and scalar cardinality are checked; real numbers
   must be finite. Independent errors are collected.
2. Read typed values, retain snapshots of supplied values, and run existing
   default/reset hooks and canonicalization.
3. Evaluate conditional requirements with supplied values taking precedence
   over prepared defaults, so reset hooks cannot hide a conflicting request.
4. Check resolved values and requirements using the same engine.

Errors identify the file, parameter, offending value, legal choices/rules, and
its manual URL. All ranks receive the error before parameter broadcasts.
File-open failures, a missing `INPUT_PARAMETERS` header, unknown names and
duplicate parameters use the same broadcast-and-exit path. Reader code must
return these errors to `read_parameters` rather than exit only on rank zero.
Normal execution and `--check-input` share this preparation boundary. Existing
checks needing structure, build capabilities or calculation context continue
at their existing boundary.

A private unset initializer does not belong to the public domain. For example,
a `{1, 2}` domain rejects explicit `-1` before defaults; omission skips the
first check and the derived value must satisfy the same domain. A documented
negative disabled mode or negative physical value has different semantics.
An effective domain can have an activation condition for unused defaults;
explicit supplied values are still checked.

## Documentation and defaults

The normal numeric/Boolean/string reader registrations enforce scalar type contracts
automatically. Explicit choice/range/requirement declarations cover system,
electronic, relaxation, MD, vdW, output, DeePKS, DFPT and TDDFT inputs, including
`scf_thr_type`, `symmetry`, `mem_saver` and the spin dependency.

When writing reusable INPUT, an inactive value outside its public domain is
omitted: internal unset markers must not be presented as legal explicit input.
The value in memory is not changed.

Choices, value predicates and conditional requirements are exported to built-in
help, generated YAML and the INPUT manual from these declarations. Generated
YAML and Markdown are build artifacts and must not be committed. Recommended
values and bounds in a default formula are not automatically legal limits.

Validation does not itself resolve defaults or make every legacy reset safe.
Default resolution should distinguish presence from numerical value, and fill
only omitted inputs. In particular `0`, `false` and a legal negative value are
explicit values. Defaults requiring structure/orbitals must be resolved when
those dependencies exist. Runtime numeric fields need not become nullable:
a separate input-preparation state can represent missing and pending values.
