#include <string>
#include <vector>
#include <stdexcept>
#ifdef __MPI
#include "module_base/parallel_common.h"
#endif

#define strvalue item.str_values[0]
#define intvalue std::stoi(item.str_values[0])
#define doublevalue std::stod(item.str_values[0])
#define boolvalue assume_as_boolean(item.str_values[0])

#ifdef __MPI
#define add_double_bcast(PARAMETER)                                                                                    \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_double(para.PARAMETER); });                  \
    }
#define add_int_bcast(PARAMETER)                                                                                       \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_int(para.PARAMETER); });                     \
    }
#define add_bool_bcast(PARAMETER)                                                                                      \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_bool(para.PARAMETER); });                    \
    }
#define add_string_bcast(PARAMETER)                                                                                    \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) { Parallel_Common::bcast_string(para.PARAMETER); });                  \
    }
#define add_doublevec_bcast(PARAMETER, N, FILL)                                                                        \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) {                                                                     \
            int _vec_size = N;                                                                                         \
            Parallel_Common::bcast_int(_vec_size);                                                                     \
            if (para.PARAMETER.size() != _vec_size)                                                                    \
                para.PARAMETER.resize(_vec_size, FILL);                                                                \
            Parallel_Common::bcast_double(para.PARAMETER.data(), _vec_size);                                           \
        });                                                                                                            \
    }
#define add_intvec_bcast(PARAMETER, N, FILL)                                                                           \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) {                                                                     \
            int _vec_size = N;                                                                                         \
            Parallel_Common::bcast_int(_vec_size);                                                                     \
            if (para.PARAMETER.size() != _vec_size)                                                                    \
                para.PARAMETER.resize(_vec_size, FILL);                                                                \
            Parallel_Common::bcast_int(para.PARAMETER.data(), _vec_size);                                              \
        });                                                                                                            \
    }
#define add_stringvec_bcast(PARAMETER, N, FILL)                                                                        \
    {                                                                                                                  \
        bcastfuncs.push_back([](Parameter& para) {                                                                     \
            int _vec_size = N;                                                                                         \
            Parallel_Common::bcast_int(_vec_size);                                                                     \
            if (para.PARAMETER.size() != _vec_size)                                                                    \
                para.PARAMETER.resize(_vec_size, FILL);                                                                \
            Parallel_Common::bcast_string(para.PARAMETER.data(), _vec_size);                                           \
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
            if (para.PARAMETER.size() != N)                                                                            \
                para.PARAMETER.resize(N, FILL);                                                                        \
        });                                                                                                            \
    }

#define add_intvec_bcast(PARAMETER, N, FILL) add_doublevec_bcast(PARAMETER, N, FILL)
#define add_stringvec_bcast(PARAMETER, N, FILL) add_doublevec_bcast(PARAMETER, N, FILL)

#endif

#define sync_string(PARAMETER)                                                                                         \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) { item.final_value << para.PARAMETER; };    \
        add_string_bcast(PARAMETER);                                                                                   \
    }
#define sync_int(PARAMETER)                                                                                            \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) { item.final_value << para.PARAMETER; };    \
        add_int_bcast(PARAMETER);                                                                                      \
    }
#define sync_double(PARAMETER)                                                                                         \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) { item.final_value << para.PARAMETER; };    \
        add_double_bcast(PARAMETER);                                                                                   \
    }
#define sync_bool(PARAMETER)                                                                                           \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) { item.final_value << para.PARAMETER; };    \
        add_bool_bcast(PARAMETER);                                                                                     \
    }
#define sync_doublevec(PARAMETER, N, FILL)                                                                             \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) {                                           \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.PARAMETER[i] << " ";                                                          \
            }                                                                                                          \
        };                                                                                                             \
        add_doublevec_bcast(PARAMETER, N, FILL);                                                                       \
    }
#define sync_intvec(PARAMETER, N, FILL)                                                                                \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) {                                           \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.PARAMETER[i] << " ";                                                          \
            }                                                                                                          \
        };                                                                                                             \
        add_intvec_bcast(PARAMETER, N, FILL);                                                                          \
    }
#define sync_stringvec(PARAMETER, N, FILL)                                                                             \
    {                                                                                                                  \
        item.get_final_value = [](Input_Item& item, const Parameter& para) {                                           \
            for (int i = 0; i < N; i++)                                                                                \
            {                                                                                                          \
                item.final_value << para.PARAMETER[i] << " ";                                                          \
            }                                                                                                          \
        };                                                                                                             \
        add_stringvec_bcast(PARAMETER, N, FILL);                                                                       \
    }

#define read_sync_string(PARAMETER)                                                                                    \
    {                                                                                                                  \
        item.read_value = [](const Input_Item& item, Parameter& para) { para.PARAMETER = strvalue; };                  \
        sync_string(PARAMETER);                                                                                        \
    }
#define read_sync_int(PARAMETER)                                                                                       \
    {                                                                                                                  \
        item.read_value = [](const Input_Item& item, Parameter& para) { para.PARAMETER = intvalue; };                  \
        sync_int(PARAMETER);                                                                                           \
    }
#define read_sync_double(PARAMETER)                                                                                    \
    {                                                                                                                  \
        item.read_value = [](const Input_Item& item, Parameter& para) { para.PARAMETER = doublevalue; };               \
        sync_double(PARAMETER);                                                                                        \
    }
#define read_sync_bool(PARAMETER)                                                                                      \
    {                                                                                                                  \
        item.read_value = [](const Input_Item& item, Parameter& para) { para.PARAMETER = boolvalue; };                 \
        sync_bool(PARAMETER);                                                                                          \
    }

/**
 * @brief To parse input parameters as expressions into vectors
 *
 * @tparam T
 * @param expressions  (vector<string>): expressions such as "3*1 0 2*0.5 3*0"
 * @param result (vector): stores parsing results,
 *            for example, "3*1 0 2*0.5 1*1.5" can be parsed as
 *            [1, 1, 1, 0, 0.5, 0.5, 1.5]
 */
template <typename T>
void parse_expression(const std::vector<std::string>& expressions, std::vector<T>& result)
{
    result.clear(); // Clear the output vector to prepare for new entries
    if (expressions.empty())
    {
        return;
    }
    else if (expressions[0].empty())
    {
        return;
    }

    for (const auto& expr: expressions)
    {
        size_t first_star_pos = expr.find('*');
        size_t second_star_pos = expr.rfind('*'); // rfind finds the rightmost '*'

        // e.g. "3", "3.5"
        // If no '*' found, convert the whole expression to double/int and add to result
        if (first_star_pos == std::string::npos)
        {
            T T_value = static_cast<T>(std::stof(expr));
            result.push_back(T_value);
        }
        // e.g. "2*3", "2*3.5"
        // If only one '*' found, split the expression into two parts and convert them to double/int
        else if (first_star_pos == second_star_pos)
        {
            std::string int_part = expr.substr(0, first_star_pos);
            std::string T_part = expr.substr(first_star_pos + 1);

            int num = std::stoi(int_part);
            T T_value = static_cast<T>(std::stof(T_part));
            for(int i = 0 ; i < num; ++i)
            {
                result.push_back(T_value);
            }
        }
        // e.g. "2*3*3" 
        // If more than one '*' found, output an error message
        else
        {
            throw std::runtime_error("Invalid expression: " + expr + " - More than one '*' found.");
        }
    }
}