#include "module_base/parallel_common.h"

#include <sstream>
#include <string>

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

#ifdef __MPI
#define bcast_string_param(PARAMETER)                                                                                  \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_string(para.PARAMETER); });                  \
    }
#define bcast_int_param(PARAMETER)                                                                                     \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_int(para.PARAMETER); });                     \
    }
#define bcast_double_param(PARAMETER)                                                                                  \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_double(para.PARAMETER); });                  \
    }
#define bcast_bool_param(PARAMETER)                                                                                    \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_bool(para.PARAMETER); });                    \
    }
#define bcast_doublevec_param(PARAMETER, N)                                                                            \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_double(para.PARAMETER, N); });               \
    }
#define bcast_intvec_param(PARAMETER, N)                                                                               \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_int(para.PARAMETER, N); });                  \
    }

#else
#define bcast_string_param(PARAMETER)
#define bcast_int_param(PARAMETER)
#define bcast_double_param(PARAMETER)
#define bcast_bool_param(PARAMETER)
#define bcast_doublevec_param(PARAMETER, N)
#define bcast_intvec_param(PARAMETER, N)

#endif

#ifdef __MPI
#define read_bcast_string(PARAMETER)                                                                                   \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.PARAMETER = strvalue; };                   \
        bcast_string_param(PARAMETER);                                                                                 \
    }
#define read_bcast_int(PARAMETER)                                                                                      \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.PARAMETER = intvalue; };                   \
        bcast_int_param(PARAMETER);                                                                                    \
    }
#define read_bcast_double(PARAMETER)                                                                                   \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.PARAMETER = doublevalue; };                \
        bcast_double_param(PARAMETER);                                                                                 \
    }
#define read_bcast_bool(PARAMETER)                                                                                     \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.PARAMETER = boolvalue; };                  \
        bcast_bool_param(PARAMETER);                                                                                   \
    }

#else
#define read_bcast_string(PARAMETER)                                                                                   \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.PARAMETER = strvalue; };                   \
    }
#define read_bcast_int(PARAMETER)                                                                                      \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.PARAMETER = intvalue; };                   \
    }
#define read_bcast_double(PARAMETER)                                                                                   \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.PARAMETER = doublevalue; };                \
    }
#define read_bcast_bool(PARAMETER)                                                                                     \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.PARAMETER = boolvalue; };                  \
    }

#endif