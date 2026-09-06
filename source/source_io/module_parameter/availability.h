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
    std::string op;                     ///< ==, !=, >, >=, <, <=, "in" or vector containment
    std::vector<std::string> values; ///< one value for comparisons, several for "in"

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

/// Parse a canonical availability string into its boolean-expression tree.
///
/// The input must exactly match the AST serialization (`param==value`,
/// `param in [a, b]`, comparison operators, `and`/`or` combinators and `(...)`
/// grouping). An empty string yields an empty (always-available) expression.
/// Invalid or non-canonical input throws std::invalid_argument.
AvailabilityExpr parse_availability(const std::string& raw);

} // namespace ModuleIO

#endif // AVAILABILITY_H
