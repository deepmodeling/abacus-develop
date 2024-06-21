#ifndef READ_INPUT_H
#define READ_INPUT_H
#include "input_item.h"
#include "module_parameter/parameter.h"

#include <string>
#include <sstream>
#ifdef __MPI
#include "mpi.h"
#endif

namespace ModuleIO
{

class ReadInput
{
  public:
    ReadInput();
    ~ReadInput(){};
    /**
     * @brief read INPUT file of txt format
     *
     * @param param parameters of ABACUS
     * @param filename INPUT file name
     */
    void ReadTxtInput(Parameter& param, const std::string& filename);
#ifdef __MPI
    /**
     * @brief initialize MPI
     * @param rank_in rank
     * @param comm_world_in MPI communicator
     *
     */
    void InitMPI(const int rank_in, MPI_Comm comm_world_in)
    {
        rank = rank_in;
        comm_world = comm_world_in;
    }
#endif
  private:
    void add_item(const Input_Item& item);
    void item_general();
    void item_md();

  private:
#ifdef __MPI
    MPI_Comm comm_world;
#endif
    int rank = 0;
    // All input items
    std::map<std::string, Input_Item> input_lists;

    // read value items for readin parameters
    std::vector<Input_Item*> readvalue_items;
    // check value items for readin parameters
    std::vector<Input_Item*> checkvalue_items;
    // reset value items for readin parameters
    std::vector<Input_Item*> resetvalue_items;
};

template <class T>
T convertstr(const std::string& in)
{
  T out;
  std::stringstream ss;
  ss << in;
  ss >> out;
  return out;
}
#define strvalue item.str_values[0]
#define intvalue convertstr<int>(item.str_values[0])
#define doublevalue convertstr<double>(item.str_values[0])
#define boolvalue convertstr<bool>(item.str_values[0])

} // namespace ModuleIO

#endif