#ifndef READ_INPUT_H
#define READ_INPUT_H

#include "input_item.h"
#include "module_parameter/parameter.h"

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
     * @brief determine the md step in restart case
     *
     * @param file_dir directory of Restart_md.dat
     * @return md step
     */
    int current_md_step(const std::string& file_dir);
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
    //set globalv parameters
    void set_globalv(Parameter& para);
    // add bcast functions for global values
    void set_globalv_bcast();
    // system items
    void item_system();
    // items for electronic structure
    void item_elec_stru();
    // items for lcao
    void item_lcao();
    // items for relax
    void item_relax();
    // items for md
    void item_md();
    // items for ofdft
    void item_ofdft();
    // items for sdft
    void item_sdft();
    // items for deepks
    void item_deepks();
    // items for real time tddft
    void item_rt_tddft();
    // items for linear response tddft
    void item_lr_tddft();
    // items for output
    void item_output();
    // items for postprocess
    void item_postprocess();
    // items for some models
    void item_model();
    // items for exx
    void item_exx();
    // items for dft+u
    void item_dftu();
    // items for other
    void item_others();

  private:
    int rank = 0; ///< rank of MPI

    // All input items
    // std::map<std::string, Input_Item> input_lists;
    // use vector instead of map to keep the order of input items
    std::vector<std::pair<std::string, Input_Item>> input_lists;

    // read value if INPUT file has this item
    // This function will be done only when INPUT file has them.
    std::vector<Input_Item*> readvalue_items;
    //-----------------------------------------------------------------

    // bcast all values function
    // if no MPI, this function will resize the vector
    // This function must be done, no matter INPUT file has them or not.
    std::vector<std::function<void(Parameter&)>> bcastfuncs;
};

// convert string vector to a long string
std::string longstring(const std::vector<std::string>& str_values);
// convert string to bool
bool assume_as_boolean(const std::string& val);
// convert to directory format
std::string to_dir(const std::string& str);
// return a warning string if the string is not found in the vector
std::string nofound_str(std::vector<std::string> init_chgs, const std::string& str);

} // namespace ModuleIO

#endif