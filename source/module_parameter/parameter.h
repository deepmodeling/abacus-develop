#ifndef PARAMETER_H
#define PARAMETER_H
#include "input_parameter.h"
#include "other_parameter.h"
namespace ModuleIO
{
class ReadInput;
}
class Parameter
{
  public:
    Parameter(){};
    ~Parameter(){};
    // ---------------------------------------------------------------
    // --------------          Getters                ----------------
    // ---------------------------------------------------------------
  public:
    // We can only read the value of input, but cannot modify it.
    const Input_para& get() const;
    // We can only read the value of mdp, but cannot modify it.
    const MD_para& get_mdp() const;
    // We can only read the value of other parameters, but cannot modify it.
    const Other_para& globalV() const;

  private:
    // Only ReadInput can modify the value of Parameter.
    friend class ModuleIO::ReadInput;
    // INPUT parameters
    Input_para input;
    // Other parameters
    Other_para gv;
};

extern Parameter PARAM;
#endif