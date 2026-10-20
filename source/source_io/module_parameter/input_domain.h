#ifndef MODULE_IO_INPUT_DOMAIN_H
#define MODULE_IO_INPUT_DOMAIN_H

#include <string>
#include <vector>

struct Input_para;

namespace ModuleIO
{

/**
 * Declarative legal values and a typed, read-only binding to an INPUT field.
 * The initial domain operator is membership in a set of scalar integers.
 * An empty domain leaves a parameter on its existing validation path.
 */
class InputDomain
{
  public:
    InputDomain() = default;
    InputDomain(int Input_para::* member, const std::vector<int>& choices);

    bool empty() const;
    const std::vector<int>& integer_choices() const;
    std::string format_choices() const;

    // Return an empty string on success, otherwise a diagnostic with legal
    // values and a manual link. Neither operation modifies input values.
    std::string validate_explicit(const std::string& label, const std::vector<std::string>& tokens) const;
    std::string validate_resolved(const std::string& label, const Input_para& input) const;

  private:
    int Input_para::* member_ = nullptr;
    std::vector<int> choices_;

    bool contains(int value) const;
    std::string diagnostic(const std::string& label, const std::string& value, const std::string& origin) const;
};

} // namespace ModuleIO

#endif // MODULE_IO_INPUT_DOMAIN_H
