#ifndef READ_INPUT_H
#define READ_INPUT_H
#include "input_item.h"
#include "module_parameter/parameter.h"

#include <sstream>
#include <string>

namespace ModuleIO
{
class ReadInput
{
  public:
    ReadInput(const int& rank);
    ~ReadInput(){};
    /**
     * @brief clear all input items
     */
    void clear()
    {
        for (auto& item: input_lists)
        {
            item.second.final_value.str("");
            item.second.str_values.clear();
        }
        readvalue_items.clear();
        checkvalue_items.clear();
        resetvalue_items.clear();
    }
    /**
     * @brief read in parameters from input file
     *
     * @param param parameters of ABACUS
     * @param filename_in read INPUT file name
     */
    void read_parameters(Parameter& param, const std::string& filename_in);

    /**
     * @brief Create a directory for output files
     * 
     * @param param parameters of ABACUS
     */
    void create_directory(const Parameter& param);

    /**
     * @brief write out parameters to output file
     *
     * @param param parameters of ABACUS
     * @param filename_out write output file name
     */
    void write_parameters(const Parameter& param, const std::string& filename_out);
    static bool check_mode;
    bool check_ntype_flag = true; ///< check ntype from STRU file

  private:
    /**
     * @brief read INPUT file of txt format
     *
     * @param param parameters of ABACUS
     * @param filename INPUT
     */
    void read_txt_input(Parameter& param, const std::string& filename);
    /**
     * @brief write INPUT file of txt format
     *
     * @param param parameters of ABACUS
     * @param filename output file name
     */
    void write_txt_input(const Parameter& param, const std::string& filename);
    /**
     * @brief count_nype from STRU file
     *
     */
    void check_ntype(const std::string& fn, int& param_ntype);
    /**
     * @brief add item to input list
     *
     * @param item input_item
     */
    void add_item(const Input_Item& item);
    // general items
    void item_general();
    // items for pw
    void item_pw();
    // items for sdft
    void item_sdft();
    // items for relax
    void item_relax();
    // items for lcao
    void item_lcao();
    // items for postprocess
    void item_postprocess();
    // items for md
    void item_md();
    // items for others
    void item_others();

  private:
    int rank = 0;
    // All input items
    // std::map<std::string, Input_Item> input_lists;
    // use vector instead of map to keep the order of input items
    std::vector<std::pair<std::string, Input_Item>> input_lists;
    //----These functions are done only when INPUT file has them.------
    // read value if INPUT file has this item
    std::vector<Input_Item*> readvalue_items;
    // check value if INPUT file has this item
    std::vector<Input_Item*> checkvalue_items;
    // reset some values if INPUT file has this item
    std::vector<Input_Item*> resetvalue_items;
    //-----------------------------------------------------------------

    //----These functions must be done----------------------
    /**
     * @brief autoset some values
     *        For "default" inputs, e.g. ks_esolver = "default", force_thr = -1, etc.
     * @note "autosetfuncs" can also serve as a fallback function for "resetvalue_items" or "checkvalue_items" because
     * it will definitely execute, but it is recommended to use "autosetfuncs" as much as possible. This will help
     * you understand the relationships between input parameters.
     */
    std::vector<std::function<void(Parameter&)>> autosetfuncs;
    /// bcast all values function
    /// if no MPI, this function will resize the vector
    std::vector<std::function<void(Parameter&)>> bcastfuncs;
    //------------------------------------------------------
};

void strtolower(char* sa, char* sb);
bool convert_bool(std::string str);
bool find_str(const std::vector<std::string>& strings, const std::string& strToFind);

} // namespace ModuleIO

#endif