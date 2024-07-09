
#include <sstream>
#include <string>

#ifdef __MPI
#include "module_base/parallel_common.h"
#endif

#define strvalue item.str_values[0]
#define intvalue std::stoi(item.str_values[0])
#define doublevalue std::stod(item.str_values[0])
#define boolvalue convert_bool(item.str_values[0])

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
#define add_doublevec_bcast(PARAMETER, N, FILL)                                                                        \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) {                                                                     \
            if (para.input.PARAMETER.size() != N)                                                                      \
                para.input.PARAMETER.resize(N, FILL);                                                                  \
            Parallel_Common::bcast_double(para.input.PARAMETER.data(), N);                                             \
        });                                                                                                            \
    }
#define add_intvec_bcast(PARAMETER, N, FILL)                                                                           \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) {                                                                     \
            if (para.input.PARAMETER.size() != N)                                                                      \
                para.input.PARAMETER.resize(N, FILL);                                                                  \
            Parallel_Common::bcast_int(para.input.PARAMETER.data(), N);                                                \
        });                                                                                                            \
    }
#define add_stringvec_bcast(PARAMETER, N, FILL)                                                                        \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) {                                                                     \
            if (para.input.PARAMETER.size() != N)                                                                      \
                para.input.PARAMETER.resize(N, FILL);                                                                  \
            Parallel_Common::bcast_string(para.input.PARAMETER.data(), N);                                             \
        });                                                                                                            \
    }

#else
#define add_double_bcast(PARAMETER)
#define add_int_bcast(PARAMETER)
#define add_bool_bcast(PARAMETER)
#define add_string_bcast(PARAMETER)
#define add_doublevec_bcast(PARAMETER, N, FILL)                                                                        \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) {                                                                     \
            if (para.input.PARAMETER.size() != N)                                                                      \
                para.input.PARAMETER.resize(N, FILL);                                                                  \
        });                                                                                                            \
    }

#define add_intvec_bcast(PARAMETER, N, FILL) add_doublevec_bcast(PARAMETER, N, FILL)
#define add_stringvec_bcast(PARAMETER, N, FILL) add_doublevec_bcast(PARAMETER, N, FILL)

#endif

#define sync_string(PARAMETER)                                                                                         \
    {                                                                                                                  \
        item.get_final_value                                                                                             \
            = [](Input_Item& item, const Parameter& para) { item.final_value << para.input.PARAMETER; };               \
        add_string_bcast(PARAMETER);                                                                                   \
    }
#define sync_int(PARAMETER)                                                                                            \
    {                                                                                                                  \
        item.get_final_value                                                                                             \
            = [](Input_Item& item, const Parameter& para) { item.final_value << para.input.PARAMETER; };               \
        add_int_bcast(PARAMETER);                                                                                      \
    }
#define sync_double(PARAMETER)                                                                                         \
    {                                                                                                                  \
        item.get_final_value                                                                                             \
            = [](Input_Item& item, const Parameter& para) { item.final_value << para.input.PARAMETER; };               \
        add_double_bcast(PARAMETER);                                                                                   \
    }
#define sync_bool(PARAMETER)                                                                                           \
    {                                                                                                                  \
        item.get_final_value                                                                                             \
            = [](Input_Item& item, const Parameter& para) { item.final_value << para.input.PARAMETER; };               \
        add_bool_bcast(PARAMETER);                                                                                     \
    }
#define sync_doublevec(PARAMETER, N, FILL)                                                                             \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) {                                             \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.input.PARAMETER[i] << " ";                                                    \
            }                                                                                                          \
        };                                                                                                             \
        add_doublevec_bcast(PARAMETER, N, FILL);                                                                       \
    }
#define sync_intvec(PARAMETER, N, FILL)                                                                                \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) {                                             \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.input.PARAMETER[i] << " ";                                                    \
            }                                                                                                          \
        };                                                                                                             \
        add_intvec_bcast(PARAMETER, N, FILL);                                                                          \
    }
#define sync_stringvec(PARAMETER, N, FILL)                                                                             \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) {                                             \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.input.PARAMETER[i] << " ";                                                    \
            }                                                                                                          \
        };                                                                                                             \
        add_stringvec_bcast(PARAMETER, N, FILL);                                                                       \
    }

#define read_sync_string(PARAMETER)                                                                                    \
    {                                                                                                                  \
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.PARAMETER = strvalue; };             \
        sync_string(PARAMETER);                                                                                        \
    }
#define read_sync_int(PARAMETER)                                                                                       \
    {                                                                                                                  \
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.PARAMETER = intvalue; };             \
        sync_int(PARAMETER);                                                                                           \
    }
#define read_sync_double(PARAMETER)                                                                                    \
    {                                                                                                                  \
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.PARAMETER = doublevalue; };          \
        sync_double(PARAMETER);                                                                                        \
    }
#define read_sync_bool(PARAMETER)                                                                                      \
    {                                                                                                                  \
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.PARAMETER = boolvalue; };            \
        sync_bool(PARAMETER);                                                                                          \
    }
