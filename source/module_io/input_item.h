#include "module_parameter/parameter.h"

#include <functional>
#include <map>
#include <string>
#include <vector>
namespace ModuleIO
{

enum ParameterType
{
    BOOL,
    INT,
    DOUBLE,
    STRING,
    VECTOR_I,
    VECTOR_D
};

class Input_Item
{
  public:
     Input_Item(){};
    // call these functions
    Input_Item(const std::string& label_in)
    {
        label = label_in;
    }

    // template <typename T>
    // void default_1(const T& value);
    // template <typename T_head, typename... T_tail>
    // void default_1(const T_head& value_head, const T_tail... values_tail);

    // set these variables and functions
    std::string annotation;
    std::function<void(const Input_Item&, Parameter&)> readvalue = [](const Input_Item& item, Parameter& param) {};
    std::function<void(const Input_Item&, const Parameter&)> checkvalue = nullptr;
    std::function<void(const Input_Item&, const Parameter&)> resetvalue = nullptr;

    std::string label;
    std::vector<std::string> str_values;
    ParameterType values_type;
};

} // namespace ModuleIO