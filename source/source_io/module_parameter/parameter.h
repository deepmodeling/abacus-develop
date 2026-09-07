#ifndef PARAMETER_H
#define PARAMETER_H
#include "input_parameter.h"
#include "system_parameter.h"
namespace ModuleIO
{
   class ReadInput;
}

namespace elecstate
{
   class ParamUpdater;
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

    // Set the rank & nproc & nthreads_per_proc
    // changed from set_rank_nproc in 2024-1018
    void set_pal_param(const int& myrank, const int& nproc, const int& nthread_per_proc);
    // Set the start time
    void set_start_time(const std::time_t& start_time);

    // ---------------------------------------------------------------
    // --------------      Test-only accessors        ----------------
    // ---------------------------------------------------------------
    // Unit tests frequently need to drive the code under test through a
    // specific INPUT configuration. These two accessors are the only
    // sanctioned way to do that: they replace the historical
    // `#define private public` hack, which reinterprets access control for
    // every declaration in the translation unit (including standard library
    // headers) and makes the test TU disagree with the rest of the build.
    //
    // Production code must use the read-only views above (inp / mdp /
    // globalv). The governance checker rejects these two names outside
    // test directories.
    Input_para& input_for_test()
    {
        return input;
    }
    System_para& sys_for_test()
    {
        return sys;
    }

  private:
    friend class ModuleIO::ReadInput; // ReadInput read INPUT file and give the value to Parameter
    friend class elecstate::ParamUpdater; // ParamUpdater updates Parameter values from atoms_info

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
