#ifndef INPUT_ITEM_H
#define INPUT_ITEM_H
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "module_parameter/parameter.h"
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

    bool is_read() const ///< check if the input item is read
    {
        return !str_values.empty();
    }

    size_t get_size() const ///< get size of the input item
    {
        if (str_values.empty())
            return 0;
        else if (str_values.size() == 1 && str_values[0].empty())
            return 0;
        else
            return str_values.size();
        return str_values.size();
    }

    std::string annotation; ///< annotation of the input item

    // ====== !!! These functions are complete.        ======
    // ====== !!! Do not add any more functions here.  ======
    /// read value function
    std::function<void(const Input_Item&, Parameter&)> readvalue = [](const Input_Item& item, Parameter& param) {};
    /// check value function
    std::function<void(const Input_Item&, const Parameter&)> checkvalue = nullptr;
    /// reset this value when some conditions are met
    std::function<void(const Input_Item&, Parameter&)> resetvalue = nullptr;
    /// get final_value function for output INPUT file
    std::function<void(Input_Item&, const Parameter&)> getfinalvalue = nullptr;
    // ====== !!! Do not add any more functions here.  ======
};

} // namespace ModuleIO
#endif // INPUT_ITEM_H