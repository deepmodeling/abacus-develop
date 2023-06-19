#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "module_io/parameter_string.h"
#include "module_io/parameter_vector.h"
enum ParameterType
{
    BOOL,
    INT,
    DOUBLE,
    STRING,
    VECTOR_I,
    VECTOR_D
};
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
        default:
            return NULL;
        }
    }
};
map<std::string, InputParameter> input_parameters;
map<std::string, string> default_parametes_type;
// std::string para_key = "lcao_ecut";
// intpu_parameters["lcao_ecut"].get()