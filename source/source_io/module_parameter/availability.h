#ifndef AVAILABILITY_H
#define AVAILABILITY_H

#include <string>
#include <vector>

namespace ModuleIO
{

/// A single condition `param op values` (e.g. `basis_type==lcao` or
/// `vdw_method in [d2, d3_0]`).
struct AvailabilityCondition
{
    std::string param;                  ///< parameter identifier
    std::string op;                     ///< ==, !=, >, >=, <, <=, "in" or "contains"
    std::vector<std::string> values;    ///< one value for comparisons, several for "in"

    std::string to_string() const;
};

/// Boolean expression tree. A leaf holds a single condition; a non-leaf node
/// holds `op` ("and"|"or") and `children`.
struct AvailabilityExpr
{
    std::string op;                       ///< "" (leaf) | "and" | "or"
    AvailabilityCondition condition;      ///< valid when leaf
    std::vector<AvailabilityExpr> children; ///< valid when non-leaf

    bool is_leaf() const
    {
        return op.empty();
    }

    std::string to_string() const;
};

/// Parse an availability string into its boolean-expression tree.
///
/// Accepts the canonical grammar (`param==value`, `param in [a, b]`, the
/// comparison operators ==, !=, >, >=, <, <=, `and`/`or`/`,` combinators and
/// `(...)` grouping). An empty string yields an empty (always-available)
/// expression.
AvailabilityExpr parse_availability(const std::string& raw);

} // namespace ModuleIO

#endif // AVAILABILITY_H
