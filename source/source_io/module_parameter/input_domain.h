#ifndef MODULE_IO_INPUT_DOMAIN_H
#define MODULE_IO_INPUT_DOMAIN_H

#include "availability_validator.h"
#include <functional>
#include <map>
#include <string>
#include <vector>

class Parameter;

namespace ModuleIO
{
// Read-only typed snapshots. Integer INPUT fields fit exactly in double storage.
struct InputValue
{
    AvailabilityValueKind kind = AvailabilityValueKind::Unknown;
    std::vector<double> numbers;
    std::vector<std::string> strings;
    InputValue() = default;
    explicit InputValue(int value);
    explicit InputValue(double value);
    explicit InputValue(bool value);
    explicit InputValue(const std::string& value);
    explicit InputValue(const std::vector<int>& value);
    explicit InputValue(const std::vector<double>& value);
    explicit InputValue(const std::vector<std::string>& value);
    std::string format() const;
};
using InputValues = std::map<std::string, InputValue>;

bool evaluate_input_condition(const AvailabilityExpr& expression, const InputValues& values);
std::string input_manual_link(const std::string& label);
std::string input_diagnostic(const std::string& label, const std::string& value,
                             const std::string& origin, const std::string& expected);

struct InputRequirement
{
    AvailabilityExpr when;
    AvailabilityExpr require;
};

class InputDomain
{
  public:
    InputDomain() = default;
    void bind(AvailabilityValueKind kind, const std::function<InputValue(const Parameter&)>& getter);
    void scalar();
    void set_choices(const std::vector<std::string>& choices);
    // Predicates use "value", e.g. "value>0 and value<=1".
    void set_rule(const std::string& rule);
    // Only effective defaults may be inactive; explicit input is always checked.
    void set_active_when(const std::string& condition);
    void add_requirement(const std::string& when, const std::string& require);

    bool empty() const;
    bool bound() const;
    bool can_write(const Parameter& input) const;
    AvailabilityValueKind kind() const;
    InputValue value(const Parameter& input) const;
    std::string format_choices() const;
    std::string format_rule() const;
    std::string format_requirements() const;
    const std::vector<std::string>& choices() const;
    void validate_schema(const std::string& label,
                         const std::map<std::string, AvailabilityValueKind>& types) const;
    std::string validate_explicit(const std::string& label, const std::vector<std::string>& tokens) const;
    std::string validate_resolved(const std::string& label, const InputValues& values) const;
    std::string validate_requirements(const std::string& label, const InputValues& values) const;

  private:
    AvailabilityValueKind kind_ = AvailabilityValueKind::Unknown;
    std::function<InputValue(const Parameter&)> getter_;
    bool scalar_ = false;
    std::vector<std::string> choices_;
    AvailabilityExpr rule_;
    AvailabilityExpr active_when_;
    std::vector<InputRequirement> requirements_;
    bool accepts(const InputValue& value) const;
    std::string expected() const;
};
} // namespace ModuleIO
#endif // MODULE_IO_INPUT_DOMAIN_H
