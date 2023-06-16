#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <vector>
#include<map>
#include "module_io/parameter_vector.h"
#include "module_io/parameter_string.h"
enum ParameterType{
    BOOL,
    INT,
    DOUBLE,
    STRING,
    VECTOR_I,
    VECTOR_D
};
class InputParameter{
    public:
        ParameterType type; // 参数类型枚举值
        union param_value{ // 参数值联合体
            bool b; // 布尔类型值
            int i; // 整型值
            double d; // 双精度浮点型值
            SimpleString s; // 简单字符串类型值
            SimpleVector<int> vi; // 简单整型向量类型值
            SimpleVector<double> vd; // 简单双精度浮点型向量类型值

            param_value(){};
            ~param_value(){};

        } value; // 参数值

        InputParameter(){
        }
        /**
         * @brief 构造函数，用于设置参数类型和默认参数值
         * 
         * @param t 参数类型枚举值
         */
        InputParameter(ParameterType t){
            type = t;
            switch(type){
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
         * @brief 设置参数值
         * 
         * @param v 参数值指针
         */
        void set(void* v){
            switch(type){
                case BOOL:{
                    value.b = *(bool*)(v);
                    break;
                } case INT:{
                    value.i = *(int*)(v);
                    break;
                } case DOUBLE:{
                    value.d = *(double*)(v);
                    break;
                } case STRING:{
                    //if(value.s == NULL) value.s = new std::string;
                    value.s = *static_cast<SimpleString*>(v);
                    break;
                } case VECTOR_I:{
                    value.vi = *static_cast<SimpleVector<int>*>(v);
                    break;
                } case VECTOR_D:{
                    value.vd = *static_cast<SimpleVector<double>*>(v);
                    break;
                }
                default: ;
            }
        }
        /**
         * @brief 获取参数值指针
         * 
         * @return void* 参数值指针
         */
        void* get(){
            switch(type){
                case BOOL: return (void*)(&value.b);
                case INT: return (void*)(&value.i);
                case DOUBLE: return (void*)(&value.d);
                case STRING: return (void*)(&value.s);
                default: return NULL;
            }
        }
};
map<std::string, InputParameter> input_parameters;
map<std::string, string> default_parametes_type;
// std::string para_key = "lcao_ecut";
// intpu_parameters["lcao_ecut"].get()