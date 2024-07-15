#ifndef PARAMETER_H
#define PARAMETER_H
#include "input_parameter.h"
#include "system_parameter.h"
namespace ModuleIO
{
class ReadInput;
}
class Parameter
{
  public:
    // Construct a new Parameter object
    Parameter(){};
    // Destruct the Parameter object
    ~Parameter(){};
    
  public:
    // ---------------------------------------------------------------
    // --------------          Getters                ----------------
    // ---------------------------------------------------------------
    
    // We can only read the value of input, but cannot modify it.
    const Input_para& inp = input;
    // We can only read the value of mdp, but cannot modify it.
    const MD_para& mdp = input.mdp;
    // We can only read the value of globalv parameters, but cannot modify it.
    const System_para& globalv = sys;

    // Set the rank & nproc
    void set_rank_nproc(const int& myrank, const int& nproc);
    // Set the start time
    void set_start_time(const std::time_t& start_time);

  private:
    // Only ReadInput can modify the value of Parameter.
    friend class ModuleIO::ReadInput;
    // INPUT parameters
    Input_para input;
    // System parameters
    System_para sys;
};

extern Parameter PARAM;

// temperarily put here
namespace GlobalV
{
extern int NPROC;
extern int MY_RANK;
extern std::ofstream ofs_running;
extern std::ofstream ofs_warning;
} // namespace GlobalV
#endif