#include "input_domain.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace ModuleIO
{
namespace
{
using Kind = AvailabilityValueKind;

std::string join(const std::vector<std::string>& values, const std::string& separator)
{
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0) { out << separator; }
        out << values[i];
    }
    return out.str();
}

InputValue parse_scalar(const Kind kind, const std::string& token)
{
    size_t consumed = 0;
    if (kind == Kind::String) { return InputValue(token); }
    if (kind == Kind::Boolean)
    {
        std::string lower = token;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
        const std::vector<std::string> yes = {"true", "1", "t", "yes", "y", "on", ".true."};
        const std::vector<std::string> no = {"false", "0", "f", "no", "n", "off", ".false."};
        if (std::find(yes.begin(), yes.end(), lower) != yes.end()) { return InputValue(true); }
        if (std::find(no.begin(), no.end(), lower) != no.end()) { return InputValue(false); }
    }
    else if (kind == Kind::Integer)
    {
        const int parsed = std::stoi(token, &consumed);
        if (consumed == token.size()) { return InputValue(parsed); }
    }
    else if (kind == Kind::Real)
    {
        const double parsed = std::stod(token, &consumed);
        if (consumed == token.size() && std::isfinite(parsed)) { return InputValue(parsed); }
    }
    throw std::invalid_argument("Invalid scalar INPUT value: " + token);
}

Kind element_kind(const Kind kind)
{
    if (kind == Kind::IntegerVector) { return Kind::Integer; }
    if (kind == Kind::RealVector) { return Kind::Real; }
    if (kind == Kind::StringVector) { return Kind::String; }
    return kind;
}

bool compare(const InputValue& actual, const std::string& op, const std::string& literal)
{
    const Kind scalar_kind = element_kind(actual.kind);
    const InputValue rhs = parse_scalar(scalar_kind, literal);
    if (scalar_kind == Kind::String)
    {
        if (op == "contains")
        {
            return std::find(actual.strings.begin(), actual.strings.end(), literal) != actual.strings.end();
        }
        const bool equal = actual.strings.size() == 1 && actual.strings.front() == literal;
        return op == "!=" ? !equal : equal;
    }
    const double wanted = rhs.numbers.front();
    if (op == "contains")
    {
        return std::find(actual.numbers.begin(), actual.numbers.end(), wanted) != actual.numbers.end();
    }
    if (actual.numbers.size() != 1) { throw std::invalid_argument("Scalar INPUT comparison requires one value."); }
    const double number = actual.numbers.front();
    if (op == ">") { return number > wanted; }
    if (op == ">=") { return number >= wanted; }
    if (op == "<") { return number < wanted; }
    if (op == "<=") { return number <= wanted; }
    if (op == "!=") { return number != wanted; }
    return number == wanted;
}

void validate_runtime_condition(const std::string& owner, const AvailabilityExpr& expression,
                                const std::map<std::string, Kind>& types)
{
    validate_availability_expr(owner, expression, types);
    if (expression.is_leaf())
    {
        const auto& condition = expression.condition;
        if (condition.param.empty()) { return; }
        const Kind kind = types.at(condition.param);
        const Kind scalar_kind = element_kind(kind);
        if (kind != scalar_kind && condition.op != "contains")
        {
            throw std::invalid_argument("Runtime INPUT vector conditions require contains: " + owner);
        }
        for (const auto& literal : condition.values) { parse_scalar(scalar_kind, literal); }
    }
    for (const auto& child : expression.children) { validate_runtime_condition(owner, child, types); }
}

void collect_references(const AvailabilityExpr& expression, std::set<std::string>& references)
{
    if (expression.is_leaf() && !expression.condition.param.empty())
    {
        references.insert(expression.condition.param);
    }
    for (const auto& child : expression.children) { collect_references(child, references); }
}

std::string describe_values(const InputValues& values, const InputRequirement& rule)
{
    std::set<std::string> references;
    collect_references(rule.when, references);
    collect_references(rule.require, references);
    std::vector<std::string> rendered;
    for (const auto& name : references) { rendered.push_back(name + "=" + values.at(name).format()); }
    return join(rendered, ", ");
}
} // namespace

InputValue::InputValue(int value) : kind(Kind::Integer)
{
    numbers.push_back(static_cast<double>(value));
}
InputValue::InputValue(double value) : kind(Kind::Real), numbers{value} {}
InputValue::InputValue(bool value) : kind(Kind::Boolean)
{
    numbers.push_back(value ? 1 : 0);
}
InputValue::InputValue(const std::string& value) : kind(Kind::String), strings{value} {}
InputValue::InputValue(const std::vector<int>& value) : kind(Kind::IntegerVector)
{
    numbers.assign(value.begin(), value.end());
}
InputValue::InputValue(const std::vector<double>& value) : kind(Kind::RealVector), numbers(value) {}
InputValue::InputValue(const std::vector<std::string>& value) : kind(Kind::StringVector), strings(value) {}

std::string InputValue::format() const
{
    if (kind == Kind::String || kind == Kind::StringVector) { return join(strings, " "); }
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (size_t i = 0; i < numbers.size(); ++i)
    {
        if (i != 0) { out << ' '; }
        if (kind == Kind::Boolean) { out << (numbers[i] != 0 ? "true" : "false"); }
        else { out << numbers[i]; }
    }
    return out.str();
}

bool evaluate_input_condition(const AvailabilityExpr& expression, const InputValues& values)
{
    if (expression.is_leaf())
    {
        const auto& condition = expression.condition;
        if (condition.param.empty()) { return true; }
        const auto& actual = values.at(condition.param);
        for (const auto& literal : condition.values)
        {
            if (compare(actual, condition.op, literal)) { return true; }
        }
        return false;
    }
    for (const auto& child : expression.children)
    {
        const bool result = evaluate_input_condition(child, values);
        if (expression.op == "and" && !result) { return false; }
        if (expression.op == "or" && result) { return true; }
    }
    return expression.op == "and";
}

std::string input_manual_link(const std::string& label)
{
    std::string anchor = label;
    std::replace(anchor.begin(), anchor.end(), '_', '-');
    return "https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#" + anchor;
}

std::string input_diagnostic(const std::string& label, const std::string& value,
                             const std::string& origin, const std::string& expected)
{
    return "Invalid INPUT parameter '" + label + "': " + origin + " '" + value + "'. " + expected
           + "\nSee " + input_manual_link(label) + "\n";
}

void InputDomain::bind(const Kind kind, const std::function<InputValue(const Parameter&)>& getter)
{
    kind_ = kind;
    getter_ = getter;
}
void InputDomain::scalar() { scalar_ = true; }
void InputDomain::set_choices(const std::vector<std::string>& choices)
{
    const std::set<std::string> unique(choices.begin(), choices.end());
    if (choices.empty() || unique.size() != choices.size())
    {
        throw std::invalid_argument("INPUT choices must be nonempty and unique.");
    }
    choices_ = choices;
    scalar_ = true;
}
void InputDomain::set_rule(const std::string& rule)
{
    rule_ = parse_availability(rule);
    scalar_ = true;
}
void InputDomain::set_active_when(const std::string& condition) { active_when_ = parse_availability(condition); }
void InputDomain::add_requirement(const std::string& when, const std::string& require)
{
    InputRequirement rule;
    rule.when = parse_availability(when);
    rule.require = parse_availability(require);
    requirements_.push_back(rule);
}
bool InputDomain::empty() const { return !scalar_; }
bool InputDomain::bound() const { return static_cast<bool>(getter_); }
bool InputDomain::can_write(const Parameter& input) const
{
    if (empty() || !bound()) { return true; }
    const InputValue resolved = value(input);
    return accepts(resolved);
}
Kind InputDomain::kind() const { return kind_; }
InputValue InputDomain::value(const Parameter& input) const { return getter_(input); }
std::string InputDomain::format_choices() const { return join(choices_, ", "); }
std::string InputDomain::format_rule() const { return rule_.to_string(); }
const std::vector<std::string>& InputDomain::choices() const { return choices_; }

std::string InputDomain::format_requirements() const
{
    std::vector<std::string> rendered;
    for (const auto& rule : requirements_)
    {
        const std::string when = rule.when.to_string();
        const std::string condition = rule.require.to_string();
        rendered.push_back(when.empty() ? condition : "When " + when + ": require " + condition);
    }
    return join(rendered, "; ");
}

void InputDomain::validate_schema(const std::string& label, const std::map<std::string, Kind>& types) const
{
    if (scalar_)
    {
        if (!bound() || (kind_ != Kind::Integer && kind_ != Kind::Real && kind_ != Kind::String && kind_ != Kind::Boolean))
        {
            throw std::invalid_argument("Missing scalar INPUT binding for " + label);
        }
        for (const auto& choice : choices_) { parse_scalar(kind_, choice); }
        const std::map<std::string, Kind> domain_types = {{"value", kind_}};
        validate_runtime_condition(label, rule_, domain_types);
    }
    validate_runtime_condition(label, active_when_, types);
    for (const auto& rule : requirements_)
    {
        validate_runtime_condition(label, rule.when, types);
        validate_runtime_condition(label, rule.require, types);
    }
}

bool InputDomain::accepts(const InputValue& value) const
{
    if (kind_ == Kind::Real && !std::isfinite(value.numbers.front())) { return false; }
    if (!choices_.empty())
    {
        bool found = false;
        for (const auto& choice : choices_) { found = found || compare(value, "==", choice); }
        if (!found) { return false; }
    }
    const InputValues values = {{"value", value}};
    return evaluate_input_condition(rule_, values);
}

std::string InputDomain::expected() const
{
    std::string type;
    if (kind_ == Kind::Integer) { type = "Integer"; }
    else if (kind_ == Kind::Real) { type = "finite Real"; }
    else if (kind_ == Kind::Boolean) { type = "Boolean (true/false, 1/0, t/f, yes/no, y/n, on/off, .true./.false.)"; }
    else { type = "String"; }
    std::string result = "Expected one " + type;
    if (!choices_.empty()) { result += " from {" + format_choices() + "}"; }
    if (!format_rule().empty()) { result += " satisfying " + format_rule(); }
    return result + ".";
}

std::string InputDomain::validate_explicit(const std::string& label, const std::vector<std::string>& tokens) const
{
    if (empty()) { return ""; }
    if (tokens.size() == 1)
    {
        try
        {
            const InputValue parsed = parse_scalar(kind_, tokens.front());
            if (accepts(parsed)) { return ""; }
        }
        catch (const std::invalid_argument&) {}
        catch (const std::out_of_range&) {}
    }
    const std::string supplied = join(tokens, " ");
    const std::string requirement = expected();
    return input_diagnostic(label, supplied, "supplied", requirement);
}

std::string InputDomain::validate_resolved(const std::string& label, const InputValues& values) const
{
    if (empty() || !evaluate_input_condition(active_when_, values)) { return ""; }
    const auto& actual = values.at(label);
    if (accepts(actual)) { return ""; }
    const std::string rendered = actual.format();
    const std::string requirement = expected();
    return input_diagnostic(label, rendered, "resolved value", requirement);
}

std::string InputDomain::validate_requirements(const std::string& label, const InputValues& values) const
{
    std::string errors;
    for (const auto& rule : requirements_)
    {
        if (evaluate_input_condition(rule.when, values) && !evaluate_input_condition(rule.require, values))
        {
            const std::string actual = describe_values(values, rule);
            const std::string expected = "Requires " + rule.require.to_string() + ".";
            errors += input_diagnostic(label, actual, "configuration", expected);
        }
    }
    return errors;
}
} // namespace ModuleIO
