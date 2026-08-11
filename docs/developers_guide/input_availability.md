# Structured INPUT availability

## Purpose

`Input_Item::availability` describes when an INPUT parameter is applicable. It
is metadata for documentation and tooling; it does not reject or alter a
user's INPUT based on this condition. Runtime diagnostics, if added, must define
how explicitly supplied parameters whose conditions are false are handled.

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
one value. Such existing INPUT forms are preserved and are not assigned a new
runtime interpretation by this metadata contract.

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

Runtime evaluation is outside this metadata contract. Any implementation must
define evaluation timing, treatment of defaults and reset values, and warning
behavior for explicitly supplied parameters.
