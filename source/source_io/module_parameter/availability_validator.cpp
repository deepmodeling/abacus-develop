#include "source_io/module_parameter/availability_validator.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <set>
#include <stdexcept>

namespace ModuleIO
{
namespace
{

bool starts_with(const std::string& value, const std::string& prefix)
{
    return value.compare(0, prefix.size(), prefix) == 0;
}

bool is_vector(const AvailabilityValueKind kind)
{
    return kind == AvailabilityValueKind::BooleanVector
           || kind == AvailabilityValueKind::IntegerVector
           || kind == AvailabilityValueKind::RealVector
           || kind == AvailabilityValueKind::StringVector;
}

AvailabilityValueKind element_kind(const AvailabilityValueKind kind)
{
    switch (kind)
    {
        case AvailabilityValueKind::BooleanVector:
            return AvailabilityValueKind::Boolean;
        case AvailabilityValueKind::IntegerVector:
            return AvailabilityValueKind::Integer;
        case AvailabilityValueKind::RealVector:
            return AvailabilityValueKind::Real;
        case AvailabilityValueKind::StringVector:
            return AvailabilityValueKind::String;
        default:
            return kind;
    }
}

bool is_integer(const std::string& value)
{
    if (value.empty())
    {
        return false;
    }
    char* end = nullptr;
    errno = 0;
    std::strtol(value.c_str(), &end, 10);
    return errno != ERANGE && end != value.c_str() && *end == '\0';
}

bool is_real(const std::string& value)
{
    if (value.empty())
    {
        return false;
    }
    char* end = nullptr;
    errno = 0;
    std::strtod(value.c_str(), &end);
    return errno != ERANGE && end != value.c_str() && *end == '\0';
}

bool literal_matches(const std::string& value, const AvailabilityValueKind kind)
{
    switch (kind)
    {
        case AvailabilityValueKind::Boolean:
            return value == "true" || value == "false" || value == "0" || value == "1";
        case AvailabilityValueKind::Integer:
            return is_integer(value);
        case AvailabilityValueKind::Real:
            return is_real(value);
        case AvailabilityValueKind::String:
            return !value.empty();
        default:
            return false;
    }
}

void fail(const std::string& owner, const std::string& message)
{
    throw std::invalid_argument("Invalid availability for '" + owner + "': " + message);
}

void validate_condition(
    const std::string& owner,
    const AvailabilityCondition& condition,
    const std::map<std::string, AvailabilityValueKind>& parameter_types)
{
    const auto parameter = parameter_types.find(condition.param);
    if (parameter == parameter_types.end())
    {
        fail(owner, "unknown parameter '" + condition.param + "'");
    }
    if (parameter->second == AvailabilityValueKind::Unknown)
    {
        fail(owner, "parameter '" + condition.param + "' has no machine-readable type");
    }

    const AvailabilityValueKind kind = parameter->second;
    if (condition.op == "contains" && !is_vector(kind))
    {
        fail(owner, "operator 'contains' requires a vector parameter, but '" + condition.param
                        + "' is scalar");
    }
    if ((condition.op == ">" || condition.op == ">=" || condition.op == "<"
         || condition.op == "<=")
        && kind != AvailabilityValueKind::Integer && kind != AvailabilityValueKind::Real)
    {
        fail(owner, "ordered comparison requires a numeric scalar parameter, but '"
                        + condition.param + "' is not one");
    }

    const AvailabilityValueKind literal_kind = element_kind(kind);
    for (const std::string& value : condition.values)
    {
        const bool legacy_string_vector
            = kind == AvailabilityValueKind::StringVector && condition.op == "in";
        if (!legacy_string_vector && !literal_matches(value, literal_kind))
        {
            fail(owner,
                 "value '" + value + "' is incompatible with parameter '" + condition.param
                     + "'");
        }
    }
}

void validate_node(
    const std::string& owner,
    const AvailabilityExpr& expression,
    const std::map<std::string, AvailabilityValueKind>& parameter_types)
{
    if (expression.is_leaf())
    {
        if (!expression.condition.param.empty())
        {
            validate_condition(owner, expression.condition, parameter_types);
        }
        return;
    }
    for (const AvailabilityExpr& child : expression.children)
    {
        validate_node(owner, child, parameter_types);
    }
}

} // namespace

AvailabilityValueKind availability_value_kind(const std::string& type)
{
    if (starts_with(type, "Vector of Boolean"))
    {
        return AvailabilityValueKind::BooleanVector;
    }
    if (starts_with(type, "Vector of Integer") || starts_with(type, "A number(ntype) of Integers")
        || starts_with(type, "Integer \\[Integer\\]"))
    {
        return AvailabilityValueKind::IntegerVector;
    }
    if (starts_with(type, "Vector of Real"))
    {
        return AvailabilityValueKind::RealVector;
    }
    if (starts_with(type, "Vector of String") || starts_with(type, "Vector of string"))
    {
        return AvailabilityValueKind::StringVector;
    }
    if (type == "Boolean")
    {
        return AvailabilityValueKind::Boolean;
    }
    if (type == "Integer")
    {
        return AvailabilityValueKind::Integer;
    }
    if (type == "Real" || type == "Float")
    {
        return AvailabilityValueKind::Real;
    }
    if (type == "String")
    {
        return AvailabilityValueKind::String;
    }
    return AvailabilityValueKind::Unknown;
}

void validate_availability_expr(
    const std::string& owner,
    const AvailabilityExpr& expression,
    const std::map<std::string, AvailabilityValueKind>& parameter_types)
{
    validate_node(owner, expression, parameter_types);
}

namespace
{

using Conjunction = std::set<std::string>;

std::string structural_key(const AvailabilityExpr& expression);

void collect_operands(const AvailabilityExpr& expression,
                      const std::string& op,
                      std::vector<std::string>& operands)
{
    if (!expression.is_leaf() && expression.op == op)
    {
        for (const AvailabilityExpr& child : expression.children)
        {
            collect_operands(child, op, operands);
        }
        return;
    }
    operands.push_back(structural_key(expression));
}

/// Build an order-independent key while flattening only associative uses of the
/// same boolean operator. No distributive or condition-level inference is done.
std::string structural_key(const AvailabilityExpr& expression)
{
    if (expression.is_leaf())
    {
        const std::string condition = expression.condition.to_string();
        return "leaf:" + std::to_string(condition.size()) + ":" + condition;
    }

    std::vector<std::string> operands;
    collect_operands(expression, expression.op, operands);
    std::sort(operands.begin(), operands.end());

    std::string result = expression.op + ":";
    for (const std::string& operand : operands)
    {
        result += std::to_string(operand.size()) + ":" + operand;
    }
    return result;
}

void add_conjuncts(const AvailabilityExpr& expression, Conjunction& conjunction)
{
    if (!expression.is_leaf() && expression.op == "and")
    {
        for (const AvailabilityExpr& child : expression.children)
        {
            add_conjuncts(child, conjunction);
        }
        return;
    }
    conjunction.insert(structural_key(expression));
}

bool conjunction_implies(const AvailabilityExpr& prerequisite,
                         const Conjunction& conjunction)
{
    if (prerequisite.is_leaf() && prerequisite.condition.param.empty())
    {
        return true;
    }
    if (conjunction.count(structural_key(prerequisite)) != 0)
    {
        return true;
    }
    if (prerequisite.is_leaf())
    {
        return false;
    }
    if (prerequisite.op == "and")
    {
        for (const AvailabilityExpr& child : prerequisite.children)
        {
            if (!conjunction_implies(child, conjunction))
            {
                return false;
            }
        }
        return true;
    }
    for (const AvailabilityExpr& child : prerequisite.children)
    {
        if (conjunction_implies(child, conjunction))
        {
            return true;
        }
    }
    return false;
}

/// Check every reference against the explicit conditions in its enclosing
/// conjunction. An AND prerequisite needs every operand, while any satisfied OR
/// branch is sufficient. Different leaf conditions are not related.
void validate_node_self_contained(
    const std::string& owner,
    const AvailabilityExpr& expression,
    const std::map<std::string, AvailabilityExpr>& expressions,
    const Conjunction& enclosing_conjunction)
{
    if (expression.is_leaf())
    {
        if (expression.condition.param.empty())
        {
            return;
        }
        const std::string& referenced = expression.condition.param;
        if (referenced == owner)
        {
            fail(owner, "references itself");
        }
        const auto it = expressions.find(referenced);
        if (it != expressions.end())
        {
            if (!conjunction_implies(it->second, enclosing_conjunction))
            {
                fail(owner,
                     "references '" + referenced + "', whose availability requires '"
                         + it->second.to_string()
                         + "'; ensure that prerequisite is implied on the same path");
            }
        }
        return;
    }

    if (expression.op == "and")
    {
        for (std::size_t i = 0; i < expression.children.size(); ++i)
        {
            Conjunction child_conjunction = enclosing_conjunction;
            for (std::size_t j = 0; j < expression.children.size(); ++j)
            {
                if (i != j)
                {
                    add_conjuncts(expression.children[j], child_conjunction);
                }
            }
            validate_node_self_contained(owner,
                                         expression.children[i],
                                         expressions,
                                         child_conjunction);
        }
        return;
    }

    for (const AvailabilityExpr& child : expression.children)
    {
        validate_node_self_contained(owner, child, expressions, enclosing_conjunction);
    }
}

} // namespace

void validate_availability_self_contained(
    const std::string& owner,
    const AvailabilityExpr& expression,
    const std::map<std::string, AvailabilityExpr>& expressions)
{
    validate_node_self_contained(owner, expression, expressions, Conjunction{});
}

} // namespace ModuleIO
