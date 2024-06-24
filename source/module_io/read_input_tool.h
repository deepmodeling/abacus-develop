
#include <sstream>
#include <string>

#ifdef __MPI
#include "module_base/parallel_common.h"
#endif

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
#define add_double_bcast(PARAMETER)                                                                                    \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_double(para.input.PARAMETER); });            \
    }
#define add_int_bcast(PARAMETER)                                                                                       \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_int(para.input.PARAMETER); });               \
    }
#define add_bool_bcast(PARAMETER)                                                                                      \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_bool(para.input.PARAMETER); });              \
    }
#define add_string_bcast(PARAMETER)                                                                                    \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_string(para.input.PARAMETER); });            \
    }
#define add_doublevec_bcast(PARAMETER, N)                                                                              \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_double(para.input.PARAMETER.data(), N); });  \
    }
#define add_intvec_bcast(PARAMETER, N)                                                                                 \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_int(para.input.PARAMETER.data(), N); });     \
    }
#define add_stringvec_bcast(PARAMETER, N)                                                                              \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_string(para.input.PARAMETER.data(), N); });  \
    }

#else
#define add_double_bcast(PARAMETER)
#define add_int_bcast(PARAMETER)
#define add_bool_bcast(PARAMETER)
#define add_string_bcast(PARAMETER)
#define add_doublevec_bcast(PARAMETER, N)
#define add_intvec_bcast(PARAMETER, N)
#define add_stringvec_bcast(PARAMETER, N)

#endif

#define sync_string(PARAMETER)                                                                                         \
    {                                                                                                                  \
        item.getfinalvalue                                                                                             \
            = [](Input_Item& item, const Parameter& para) { item.final_value << para.input.PARAMETER; };               \
        add_string_bcast(PARAMETER);                                                                                   \
    }
#define sync_int(PARAMETER)                                                                                            \
    {                                                                                                                  \
        item.getfinalvalue                                                                                             \
            = [](Input_Item& item, const Parameter& para) { item.final_value << para.input.PARAMETER; };               \
        add_int_bcast(PARAMETER);                                                                                      \
    }
#define sync_double(PARAMETER)                                                                                         \
    {                                                                                                                  \
        item.getfinalvalue                                                                                             \
            = [](Input_Item& item, const Parameter& para) { item.final_value << para.input.PARAMETER; };               \
        add_double_bcast(PARAMETER);                                                                                   \
    }
#define sync_bool(PARAMETER)                                                                                           \
    {                                                                                                                  \
        item.getfinalvalue                                                                                             \
            = [](Input_Item& item, const Parameter& para) { item.final_value << para.input.PARAMETER; };               \
        add_bool_bcast(PARAMETER);                                                                                     \
    }
#define sync_doublevec(PARAMETER, N)                                                                                   \
    {                                                                                                                  \
        item.getfinalvalue = [](Input_Item& item, const Parameter& para) {                                             \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.input.PARAMETER[i] << " ";                                                    \
            }                                                                                                          \
        };                                                                                                             \
        add_doublevec_bcast(PARAMETER, N);                                                                             \
    }
#define sync_intvec(PARAMETER, N)                                                                                      \
    {                                                                                                                  \
        item.getfinalvalue = [](Input_Item& item, const Parameter& para) {                                             \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.input.PARAMETER[i] << " ";                                                    \
            }                                                                                                          \
        };                                                                                                             \
        add_intvec_bcast(PARAMETER, N);                                                                                \
    }
#define sync_stringvec(PARAMETER, N)                                                                                   \
    {                                                                                                                  \
        item.getfinalvalue = [](Input_Item& item, const Parameter& para) {                                             \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.input.PARAMETER[i] << " ";                                                    \
            }                                                                                                          \
        };                                                                                                             \
        add_stringvec_bcast(PARAMETER, N);                                                                             \
    }

#define read_sync_string(PARAMETER)                                                                                    \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.input.PARAMETER = strvalue; };             \
        sync_string(PARAMETER);                                                                                        \
    }
#define read_sync_int(PARAMETER)                                                                                       \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.input.PARAMETER = intvalue; };             \
        sync_int(PARAMETER);                                                                                           \
    }
#define read_sync_double(PARAMETER)                                                                                    \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.input.PARAMETER = doublevalue; };          \
        sync_double(PARAMETER);                                                                                        \
    }
#define read_sync_bool(PARAMETER)                                                                                      \
    {                                                                                                                  \
        item.readvalue = [](const Input_Item& item, Parameter& para) { para.input.PARAMETER = boolvalue; };            \
        sync_bool(PARAMETER);                                                                                          \
    }
