#ifndef DRIVER_H
#define DRIVER_H

#include <string>

namespace ModuleIO
{
struct CommandLineArgs;
}

class Driver
{
  public:
    Driver();
    ~Driver();

    /**
     * @brief An interface function.
     * This function calls "this->reading()", "this->atomic_world()" in order.
     *
     * The selected INPUT path is passed explicitly without global state.
     * Parsed variables are reserved for PART II and are not applied here.
     */
    void init(const ModuleIO::CommandLineArgs& cli);

  private:
    /**
     * @brief Print the start information.
     *
     */
    void print_start_info(const std::string& input_card);
    /**
     * @brief Reading the parameters and split the MPI world.
     * This function read the parameter in "INPUT", "STRU" etc,
     * and split the MPI world into different groups.
     *
     * The selected INPUT path is forwarded to ReadInput here.
     */
    void reading(const ModuleIO::CommandLineArgs& cli);

    /**
     * @brief An interface function.
     * This function calls "this->driver_run()" to do calculation,
     * and log the time and  memory consumed during calculation.
     */
    void atomic_world(const ModuleIO::CommandLineArgs& cli);

    // the actual calculations
    void driver_run();

    // Init harewares according to Input parameters
    void init_hardware();
    void finalize_hardware();
};

#endif
