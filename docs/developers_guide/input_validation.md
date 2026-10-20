# Shared INPUT domain validation

`Input_Item::domain` declares legal values once for runtime checks, built-in
help and the generated INPUT reference. The initial operator supports scalar
integer choices. For example:

```cpp
item.type = "Integer";
item.domain = InputDomain(&Input_para::nspin, {1, 2, 4});
```

The field binding is typed and read-only. The validator does not branch on
parameter names and does not format final values and parse them back. Adding
another integer-choice parameter requires its binding and choices, not a new
validation callback. A domain must have a non-null binding, at least one choice,
no duplicate choices and matching `Integer` metadata.

## Preparation order

After collecting INPUT tokens, the reader calls `prepare_input` once:

1. Reject missing values for supplied keywords. Validate all explicitly supplied
   domains, including complete integer tokens and exactly one value. Collect
   independent errors before running any value reader or default/reset hook.
2. Run the existing readers, contextual default/reset hooks and alias
   normalization.
3. Validate the effective values, including omitted parameters' defaults,
   through the same domains.

Errors identify the file, parameter, supplied or resolved value, legal choices
and manual section. An invalid explicit value cannot be made acceptable by a
later reset. The domain validator never changes values. Domain errors are
broadcast before all ranks stop, ahead of parameter broadcasts. Normal startup
and `--check-input` use this same path.

`nspin` and `gga_grad` are the first migrated parameters. Their legal choices
and defaults are unchanged. Inputs such as `nspin 2abc`, `nspin 2 4` and a bare
`nspin` are now rejected instead of accepting a numeric prefix, ignoring extra
tokens or overlooking a final keyword without a value. Legal explicit zero
remains a value, not a request for a default.

## Documentation

The C++ registration is the source of truth. `--generate-parameters-yaml`
exports a typed `allowed_values` sequence, and `generate_input_main.py` renders
it as **Allowed values**. `abacus -h <parameter>` uses the same domain. Generated
YAML and Markdown must not be committed.

Only hard legal values belong in the domain. Recommended values and bounds
appearing in default formulas must not be converted into restrictions without
checking their documented meaning.

## Incremental scope

This first slice establishes shared domain checks, not a replacement for all
existing input processing. Parameters without a declared domain retain their
existing behavior. Legacy cross-parameter checks, including the requirement
for `nspin=4` with `noncolin` or `lspinorb`, still run after the required context
is available. [Availability](input_availability.md) describes applicability;
it is not automatically a hard runtime constraint.

Further migrations should extend common domain/condition operators and add
rule data. They should also separate explicit inputs from computed defaults so
that defaults can only fill omitted fields; the existing reset callbacks are
not yet restricted by this PR. In particular, this initial domain check does
not detect a legacy reset changing one legal value into a different legal value.
Ranges, conditional requirements, vector lengths and structured source-line
diagnostics are subsequent extensions of the shared engine, not reasons to
introduce parameter-specific validators.
