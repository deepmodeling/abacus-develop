# RFC: structured INPUT availability

Status: Phase 1 contract proposed by PR #7783.

## Purpose

`Input_Item::availability` describes when an INPUT parameter is applicable. It
is metadata for documentation and tooling; Phase 1 does not reject or alter a
user's INPUT from this condition. A later runtime workflow may emit warnings
for explicitly supplied parameters whose conditions are false.

## Invariants

- The C++ `Input_Item` registration is the source of truth. YAML and Markdown
  are generated artifacts.
- Every non-empty registration is canonical. `set_availability()` rejects
  syntax errors and non-canonical spelling.
- The AST is the stored representation; the exported string is serialized from
  it, so the two forms cannot diverge.
- Every expression is a complete, independently evaluable predicate. It must
  include enclosing requirements rather than inheriting them implicitly from a
  referenced parameter.
- After all INPUT items are registered, every referenced label, operator and
  literal is checked against machine-readable parameter type information.

## Grammar and meaning

```text
expression := or-expression
or-expression := and-expression ("or" and-expression)*
and-expression := primary (("and" | ",") primary)*
primary := condition | "(" expression ")"
condition := parameter comparison value
           | parameter "in" "[" value ("," value)* "]"
           | parameter "contains" value
comparison := "==" | "!=" | ">" | ">=" | "<" | "<="
```

`and` binds more tightly than `or`. `==` has exactly one right-hand value;
alternatives use `in [...]`. `/` is an ordinary value character, not another
spelling of membership. Ordered comparisons require a numeric scalar, and
`contains` requires a vector.

Some existing INPUT parameters accept multiple whitespace-separated tokens as
one value. In particular, PR #6517 changed `relax_method` from `string` to
`vector<string>`, so `relax_method in [cg 2]` preserves that existing two-token
INPUT form. It does not introduce a new runtime interpretation in Phase 1.

Examples:

```cpp
item.set_availability("basis_type==pw");
item.set_availability("vdw_method in [d2, d3_0]");
item.set_availability("td_ttype contains 2");
item.set_availability("esolver_type==sdft and method_sto==2");
```

## Registration and validation workflow

1. Parse and require canonical spelling in `set_availability()`.
2. Finish registering all `Input_Item` objects.
3. Validate referenced parameter names, operator compatibility and literals.
4. Serialize the AST into `docs/parameters.yaml`.
5. Generate `docs/advanced/input_files/input-main.md` from that YAML.

Runtime evaluation is intentionally deferred. Before it is enabled, the
project must define evaluation timing, treatment of defaults and reset values,
and warning behavior for explicitly supplied parameters.
