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
     * @brief read in parameters from input file
     *
     * @param param parameters of ABACUS
     * @param filename_in read INPUT file name
     * @param filename_out write output INPUT file name
     */
    void readin_parameters(Parameter& param, const std::string& filename_in, const std::string& filename_out = "");

    static bool check_mode;

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
    std::map<std::string, Input_Item> input_lists;

    // read value items for readin parameters
    std::vector<Input_Item*> readvalue_items;
    // check value items for readin parameters
    std::vector<Input_Item*> checkvalue_items;
    // reset value items for readin parameters
    std::vector<Input_Item*> resetvalue_items;
#ifdef __MPI
    /// bcast value function
    std::vector<std::function<void(Parameter&)>> bcastfuncs;
#endif
};

} // namespace ModuleIO

#endif