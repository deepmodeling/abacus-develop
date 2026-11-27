# Shared INPUT validation

Declare scalar parameter constraints on `Input_Item::domain` when registering
an INPUT parameter. Standard scalar readers already check the type. Parameters
with custom input syntax keep their own parsing rules.

| Constraint | Example |
| --- | --- |
| A fixed set of choices | `item.domain.set_choices({"cpu", "gpu"});` |
| A numerical range | `item.domain.set_rule("value>0");` |
| A dependency on another parameter | `item.domain.add_requirement("device_xc==gpu", "device==gpu");` |

Use `item.set_availability(...)` to describe when a parameter applies. Use
`add_requirement(...)` when an incompatible combination should produce an
error. Both use the [same condition syntax](https://github.com/deepmodeling/abacus-develop/blob/develop/docs/developers_guide/input_availability.md).

## Preparation order

Normal runs and `--check-input` use the same validation. The reader checks
supplied values before applying defaults. It then checks dependencies against
the user's requested values and validates the final values and their dependencies.
Keep checks that need the structure, build capabilities or calculation state
where that information is available.

Return input errors to `read_parameters` so all MPI ranks receive the same
diagnostic. Shared validation reports the offending value, the expected values
or rule, and a link to the parameter manual.

## Documentation and defaults

Declare values that users may supply. Keep internal markers for missing values
out of the list of choices. Apply defaults only to omitted parameters and
preserve explicit values such as `0` or `false`.

Built-in help and the generated INPUT manual use the same declarations.
Update the C++ registration when a constraint changes.
