#include "module_parameter/parameter.h"

#include <functional>
#include <map>
#include <string>
#include <vector>
namespace ModuleIO
{
class Input_Item
{
  public:
    Input_Item(){};

    Input_Item(const std::string& label_in)
    {
        label = label_in;
    }

    std::string label;                   ///< label of the input item
    std::vector<std::string> str_values; ///< string values of the input item

    std::string annotation; ///< annotation of the input item

    /// read value function from input file
    std::function<void(const Input_Item&, Parameter&)> readvalue = [](const Input_Item& item, Parameter& param) {};
    /// check the value read from input file
    std::function<void(const Input_Item&, const Parameter&)> checkvalue = nullptr;
    /// reset other values if this value is read in and set.
    std::function<void(const Input_Item&, Parameter&)> resetvalue = nullptr;
};

} // namespace ModuleIO