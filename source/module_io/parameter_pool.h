#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "module_base/constants.h"
#include "module_base/global_file.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_common.h"
#include "module_base/timer.h"
#include "module_io/input.h"
#include "module_io/parameter_string.h"
#include "module_io/parameter_vector.h"
#include "version.h"
enum ParameterType
{
    BOOL,
    INT,
    DOUBLE,
    STRING,
    VECTOR_I,
    VECTOR_D
};
namespace ModuleIO
{

class InputParameter
{
  public:
    ParameterType type; // Parameter Type Enumeration value
    union param_value { // Parameter value association
        bool b;
        int i;
        double d;
        SimpleString s;          // Simple string type value
        SimpleVector<int> vi;    // Simple integer vector type value
        SimpleVector<double> vd; // Simple double precision floating-point vector type value

        param_value(){};
        ~param_value(){};

    } value;

    InputParameter()
    {
    }
    /**
     * @brief Constructor to set parameter types and default parameter values
     *
     * @param t Parameter Type Enumeration value
     */
    InputParameter(ParameterType t)
    {
        type = t;
        switch (type)
        {
        case BOOL:
            value.b = false;
            break;
        case INT:
            value.i = 0;
            break;
        case DOUBLE:
            value.d = 0.0;
            break;
        case STRING:
            value.s = SimpleString();
            break;
        case VECTOR_I:
            value.vi = SimpleVector<int>();
            break;
        case VECTOR_D:
            value.vd = SimpleVector<double>();
            break;
        default:
            break;
        }
    }
    /**
     * @brief Overloads the assignment operator to copy the values of another InputParameter object to this object.
     *
     * @param other The InputParameter object to be copied.
     * @return A reference to this InputParameter object.
     */
    InputParameter& operator=(const InputParameter& other)
    {
        if (this != &other)
        {
            type = other.type;
            switch (type)
            {
            case BOOL:
                value.b = other.value.b;
                break;
            case INT:
                value.i = other.value.i;
                break;
            case DOUBLE:
                value.d = other.value.d;
                break;
            case STRING:
                value.s = other.value.s;
                break;
            case VECTOR_I:
                value.vi = other.value.vi;
                break;
            case VECTOR_D:
                value.vd = other.value.vd;
                break;
            default:
                break;
            }
        }
        return *this;
    }
    /**
     * @brief Copy constructor to create a new InputParameter object with the same values as another InputParameter
     * object.
     *
     * @param other The InputParameter object to be copied.
     */
    InputParameter(const InputParameter& other)
    {
        type = other.type;
        switch (type)
        {
        case BOOL:
            value.b = other.value.b;
            break;
        case INT:
            value.i = other.value.i;
            break;
        case DOUBLE:
            value.d = other.value.d;
            break;
        case STRING:
            value.s = other.value.s;
            break;
        case VECTOR_I:
            value.vi = other.value.vi;
            break;
        case VECTOR_D:
            value.vd = other.value.vd;
            break;
        default:
            break;
        }
    }
    /**
     * @brief   to free any dynamically allocated memory and destroy the InputParameter object.
     *          If the parameter type is STRING, manually call the destructor of SimpleString.
     *          If the parameter type is VECTOR_I or VECTOR_D, manually call the destructor of SimpleVector.
     *          No need to do anything for BOOL, INT, and DOUBLE.
     */
    ~InputParameter()
    {
        // If the parameter type is STRING, manually call the destructor of SimpleString
        if (type == STRING)
        {
            value.s.~SimpleString();
        }
        // If the parameter type is VECTOR_I or VECTOR_D, manually call the destructor of SimpleVector
        else if (type == VECTOR_I)
        {
            value.vi.~SimpleVector();
        }
        else if (type == VECTOR_D)
        {
            value.vd.~SimpleVector();
        }
        // No need to do anything for BOOL, INT, and DOUBLE
    }

    /**
     * @brief Set parameter values
     *
     * @param v Parameter value pointer
     */
    void set(void* v)
    {
        switch (type)
        {
        case BOOL: {
            value.b = *(bool*)(v);
            break;
        }
        case INT: {
            value.i = *(int*)(v);
            break;
        }
        case DOUBLE: {
            value.d = *(double*)(v);
            break;
        }
        case STRING: {
            // if(value.s == NULL) value.s = new std::string;
            value.s = *static_cast<SimpleString*>(v);
            break;
        }
        case VECTOR_I: {
            value.vi = *static_cast<SimpleVector<int>*>(v);
            break;
        }
        case VECTOR_D: {
            value.vd = *static_cast<SimpleVector<double>*>(v);
            break;
        }
        default:;
        }
    }
    /**
     * @brief Gets a pointer to the parameter value
     *
     * @return void* pointer
     */
    void* get()
    {
        switch (type)
        {
        case BOOL:
            return (void*)(&value.b);
        case INT:
            return (void*)(&value.i);
        case DOUBLE:
            return (void*)(&value.d);
        case STRING:
            return (void*)(&value.s);
        case VECTOR_I:
            return (void*)(&value.vi);
        case VECTOR_D:
            return (void*)(&value.vd);
        default:
            return NULL;
        }
    }
};
bool Init(const std::string& default_type_path,
          const std::string& default_value_path,
          const std::string& input_value_path);
bool default_parametes_reader(const std::string& fn, std::map<std::string, std::string>& default_parametes_type);
bool input_parameters_get(const std::string& fn, std::map<std::string, InputParameter>& input);
bool input_parameters_set(std::map<std::string, InputParameter> input_parameters);

extern std::map<std::string, InputParameter> input_parameters;
extern std::map<std::string, std::string> default_parametes_type;
extern std::map<std::string, InputParameter> default_parametes_value;
} // namespace ModuleIO

// std::string para_key = "lcao_ecut";
// intpu_parameters["lcao_ecut"].get()