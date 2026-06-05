#ifndef INPUT_ITEM_H
#define INPUT_ITEM_H
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "source_io/module_parameter/parameter.h"
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
        // Documentation fields
        category = item.category;
        type = item.type;
        description = item.description;
        default_value = item.default_value;
        unit = item.unit;
        availability = item.availability;
        annotation = item.annotation;
        read_value = item.read_value;
        check_value = item.check_value;
        reset_value = item.reset_value;
        get_final_value = item.get_final_value;
    }

    std::string label;                   ///< label of the input item
    std::vector<std::string> str_values; ///< string values of the input item
    std::stringstream final_value;       ///< final value for writing to output INPUT file

    // Documentation fields for built-in help system
    std::string category;      ///< category for grouping (e.g., "System variables")
    std::string type;          ///< data type ("Integer", "Real", "String", "Boolean")
    std::string description;   ///< full description (supports multi-line, lists, notes)
    std::string default_value; ///< default value as string
    std::string unit;          ///< unit of measurement (empty if none)
    std::string availability;  ///< availability conditions (empty if always)

    bool is_read() const ///< check if the input item is read
    {
        return !str_values.empty();
    }

    size_t get_size() const ///< get size of the input item
    {
        return str_values.size();
    }

    std::string annotation; ///< brief annotation (kept for backward compatibility)

    // ====== !!! These functions are complete.        ======
    // ====== !!! Do not add any more functions here.  ======
    /// read value function
    std::function<void(const Input_Item&, Parameter&)> read_value = [](const Input_Item& item, Parameter& param) {};
    /// check value function
    std::function<void(const Input_Item&, const Parameter&)> check_value = nullptr;
    /// reset this value when some conditions are met
    /// e.g. should only reset the value of this item
    std::function<void(const Input_Item&, Parameter&)> reset_value = nullptr;
    /// get final_value function for output INPUT file
    std::function<void(Input_Item&, const Parameter&)> get_final_value = nullptr;
    // ====== !!! Do not add any more functions here.  ======
};

} // namespace ModuleIO
#endif // INPUT_ITEM_H