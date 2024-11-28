#ifndef PARAMETER_H
#define PARAMETER_H
#include "input_parameter.h"
#include "system_parameter.h"
namespace ModuleIO
{
class ReadInput;
}
class CalAtomInfo;
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

    // set input.nbands 
    void set_input_nbands(const int& nbands);
    // set sys.nlocal
    void set_sys_nlocal(const int& nlocal);

  private:
    // Only ReadInput and CalAtomInfo can modify the value of Parameter.
    // Do not add extra friend class here!!!
    friend class ModuleIO::ReadInput; // ReadInput read INPUT file and give the value to Parameter
    friend class CalAtomsInfo; // CalAtomInfo calculate the atom information from pseudopotential and give the value to
                               // Parameter

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