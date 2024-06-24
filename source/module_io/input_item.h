#include "module_parameter/parameter.h"

#include <functional>
#include <map>
#include <string>
#include <sstream>
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

    Input_Item(const Input_Item& item)
    {
        label = item.label;
        str_values = item.str_values;
        final_value.str(item.final_value.str());
        annotation = item.annotation;
        readvalue = item.readvalue;
        checkvalue = item.checkvalue;
        resetvalue = item.resetvalue;
        getfinalvalue = item.getfinalvalue;
    }

    std::string label;                   ///< label of the input item
    std::vector<std::string> str_values; ///< string values of the input item
    std::stringstream final_value;       ///< final value for writing to output INPUT file

    std::string annotation; ///< annotation of the input item

    /// read value function from input file
    std::function<void(const Input_Item&, Parameter&)> readvalue = [](const Input_Item& item, Parameter& param) {};
    /// check the value read from input file
    std::function<void(const Input_Item&, const Parameter&)> checkvalue = nullptr;
    /// reset other values if this value is read in and set.
    std::function<void(const Input_Item&, Parameter&)> resetvalue = nullptr;
    /// get final_value function
    std::function<void(Input_Item&, const Parameter&)> getfinalvalue = nullptr;
};

} // namespace ModuleIO