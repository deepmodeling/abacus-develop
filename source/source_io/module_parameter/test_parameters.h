#ifndef TEST_PARAMETERS_H
#define TEST_PARAMETERS_H

#include "parameter.h"

/// @brief Test-only mutable access to a Parameter object.
///
/// `Parameter` exposes its INPUT and derived-system data to the rest of the code
/// through const references (`PARAM.inp`, `PARAM.mdp`, `PARAM.globalv`), so
/// production code cannot mutate them; the underlying `input` and `sys` members
/// are private and only `ModuleIO::ReadInput`, `elecstate::ParamUpdater` and this
/// class are friends of `Parameter`.
///
/// Unit tests do legitimately need to set parameter values before exercising the
/// code under test. They used to do that with
///
///     #define private public
///     #include "source_io/module_parameter/parameter.h"
///     #undef private
///
/// which compiles `Parameter` with a different definition in the test translation
/// unit than in the rest of the library that the test links against -- one class
/// with two definitions in one program, i.e. an ODR violation. It also disabled
/// access control for every other header that happened to be included inside the
/// same region.
///
/// Include this header instead. `Parameter` declares this class a friend for
/// exactly this purpose.
///
/// @warning Test code only. Do not include this header from production sources.
class TestParameters
{
  public:
    /// Mutable view of the global PARAM's INPUT parameters.
    static Input_para& input()
    {
        return PARAM.input;
    }

    /// Mutable view of the global PARAM's derived system parameters.
    static System_para& sys()
    {
        return PARAM.sys;
    }

    /// Mutable view of a local Parameter's INPUT parameters.
    static Input_para& input(Parameter& param)
    {
        return param.input;
    }

    /// Mutable view of a local Parameter's derived system parameters.
    static System_para& sys(Parameter& param)
    {
        return param.sys;
    }
};

#endif // TEST_PARAMETERS_H
