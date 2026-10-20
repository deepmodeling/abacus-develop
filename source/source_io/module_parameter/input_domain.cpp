#include "input_domain.h"
#include "input_parameter.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace ModuleIO
{

InputDomain::InputDomain(int Input_para::* member, const std::vector<int>& choices)
    : member_(member), choices_(choices)
{
    if (member_ == nullptr || choices_.empty())
    {
        throw std::invalid_argument("An integer INPUT domain requires a field binding and nonempty legal values.");
    }
    std::vector<int> unique_choices = choices_;
    std::sort(unique_choices.begin(), unique_choices.end());
    if (std::adjacent_find(unique_choices.begin(), unique_choices.end()) != unique_choices.end())
    {
        throw std::invalid_argument("An integer INPUT domain must not contain duplicate legal values.");
    }
}

bool InputDomain::empty() const
{
    return member_ == nullptr;
}

const std::vector<int>& InputDomain::integer_choices() const
{
    return choices_;
}

std::string InputDomain::format_choices() const
{
    std::ostringstream out;
    for (size_t i = 0; i < choices_.size(); ++i)
    {
        if (i != 0)
        {
            out << ", ";
        }
        out << choices_[i];
    }
    return out.str();
}

bool InputDomain::contains(int value) const
{
    return std::find(choices_.begin(), choices_.end(), value) != choices_.end();
}

std::string InputDomain::diagnostic(const std::string& label,
                                  const std::string& value,
                                  const std::string& origin) const
{
    std::string anchor = label;
    std::replace(anchor.begin(), anchor.end(), '_', '-');
    return "Invalid INPUT parameter '" + label + "': " + origin + " '" + value
           + "'. Expected one Integer from {" + format_choices() + "}.\n"
           + "See https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#" + anchor + "\n";
}

std::string InputDomain::validate_explicit(const std::string& label, const std::vector<std::string>& tokens) const
{
    if (empty())
    {
        return "";
    }
    if (tokens.size() == 1)
    {
        try
        {
            size_t consumed = 0;
            const int value = std::stoi(tokens.front(), &consumed);
            if (consumed == tokens.front().size() && contains(value))
            {
                return "";
            }
        }
        catch (const std::invalid_argument&)
        {
            // Report malformed input using the same domain diagnostic.
        }
        catch (const std::out_of_range&)
        {
            // An integer overflow is an input error, not a default request.
        }
    }
    std::ostringstream supplied;
    for (size_t i = 0; i < tokens.size(); ++i)
    {
        if (i != 0)
        {
            supplied << ' ';
        }
        supplied << tokens[i];
    }
    const std::string value = supplied.str();
    return diagnostic(label, value, "supplied");
}

std::string InputDomain::validate_resolved(const std::string& label, const Input_para& input) const
{
    if (empty())
    {
        return "";
    }
    const int value = input.*member_;
    if (contains(value))
    {
        return "";
    }
    const std::string rendered = std::to_string(value);
    return diagnostic(label, rendered, "resolved value");
}

} // namespace ModuleIO
