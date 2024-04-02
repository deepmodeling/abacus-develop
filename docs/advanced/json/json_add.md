# Table of Contents

1. [Abacus-Json Usage Instructions](#1-abacus-json-usage-instructions)
   - [Normal Usage](#normal-usage)
     - [Add/Modify a value to object json node](#addmodify-a-value-to-object-json-node-key2-is-a-object-node)
     - [Pushback a value to array json node](#pushback-a-value-to-array-json-node-key2-is-a-array-node)
   - [Initialization and Assignment Functions for Different Value Types in Arrays](#initialization-and-assignment-functions-for-different-value-types-in-arrays)
     - [Object Type](#object-type)
     - [Array Type](#array-type)
   - [Array Modification Instructions](#array-modification-instructions)
2. [Json Codes Addition Guidelines](#2-json-codes-addition-guidelines)
   - [Abacus JSON Functionality Code Structure](#abacus-json-functionality-code-structure)
   - [Add JSON code principles](#add-json-code-principles)



# 1. Abacus-Json Usage Instructions

In Abacus, the main utility functions for manipulating JSON trees are outlined below. These functions are used to add objects to Abacus JSON trees.

Function signature:
void AbacusJson::add_json (std::vector<std::string> keys, const T& value,bool IsArray)
Where:
- `keys` is a vector of string dictionaries, representing the node paths where values are to be added in the JSON tree.
- `value` is a generic value type, including int, bool, double, string, or rapidjson value type, indicating the value to be added to the Abacus JSON.
- `IsArray` is a boolean object, indicating whether the current node being added is an array. `true` represents an array node, while `false` represents a non-array node.

Example usage:
const std::string version = "v3.5.2";
AbacusJson::add_json({"general_info", "version"}, version, false);



## Normal usage

### Add/Modify a value to object json node (key2 is a object node): 
```cpp
Json::AbacusJson::add_json({"key1","key2"}, 3.1415,false);
```

### Pushback a value to array json node (key2 is a array node):
```cpp
Json::AbacusJson::add_json({"key1","key2"}, 3.1415,true);
```

Through this function alone, the addition of the majority of JSON parameters can be achieved. However, for complex array types, additional operations are required.


## Initialization and Assignment Functions for Different Value Types in Arrays

### Object Type:
Since the object type consists of key-value pairs, four member functions are divided based on whether key and val are of type std::string.

- JaddStringV(str,val): key is not string, val is string
- JaddStringK(str,val): key is string, val is not string
- JaddStringKV(str,val): both key and val are string
- JaddNormal(str,val): both key and val are not string


### Array Type:
For array types, the following member functions are used directly.

- JPushBack(val): val is not string
- JPushBackString(val): val is string

For example, to add nodes to a JSON tree with multiple arrays in Abacus, the following code is needed:

```cpp
// add key-val to an object array
for(int i=0;i<1;i++){
    Json::jsonValue object(JobjectType);
    std::string str = std::to_string(i*100);  
    
    object.JaddNormal("int",i);
    object.JaddStringV("string", str);
    Json::AbacusJson::add_json({"array"}, object,true);
}
```

```cpp
// add array in array
Json::jsonValue object0(JarrayType);

object0.JPushBack(1);
object0.JPushBack(2);
object0.JPushBack(3);

Json::AbacusJson::add_json({"Darray"}, object0,true);
```





## Array Modification Instructions

For values that need to be modified in arrays, the following method can be used:
- The index number of the array starts at 0, if it's negative, it's going from back to front. eg. If the index is -1, it means that the last element of the array is modified:
- If the path contains an array, use the array index directly.

```cpp
AbacusJson::add_json({"path",index }, value, is_array);
```

Here, index is a number. index >= 0 indicates the index from the beginning of the array, while index < 0 indicates traversal from the end of the array.

For example, to modify the value of "vasp" to "cp2k" in the following JSON tree:

```json
"Json":{
    "key6": {
        "key7": [
            {
                "a":1,
                "new":2
            }
            "vasp",
            "abacus"
        ]
    }
}
```
The relative path of "vasp" in layman's terms is Json - key6 - key7[0]. To use the JSON modification method in abacus, simply change the index [0] to "0".

```cpp
AbacusJson::add_json({"Json","key6","key7",1}, "cp2k" , false);
```

If traversal is done from the end:
```cpp
AbacusJson::add_json({"Json","key6","key7",-2}, "cp2k", false);
```

An error is reported if index exceeds the array length!
```cpp
AbacusJson::add_json({"Json","key6","key7",3}, "cp2k", false);
```




# 2. Abacus Json Codes Addition Guidelines

## Abacus JSON Functionality Code Structure

The current code structure of JSON functionality in Abacus is roughly as follows:

- source/module_io
  - para_json.cpp: Contains JSON generation and output interfaces directly called by the device in Abacus.
  - json_output/: Contains the functionality encapsulation class `abacusjson.cpp` of RapidJSON in Abacus and code classes for parameter generation in various JSON modules.
    - test: Code testing files in `json_output`.


## Add JSON code principles:
In Abacus JSON addition, the following principles need to be followed:

1. Whenever possible, code to be added in the module should be written in the `json_output` module (there may also be cases where parameters cannot be directly obtained through parameter passing in `json_output`), and then called in the path `para_json.cpp` -> `device.cpp` or other main execution paths. (Ensure minimal impact on other modules as much as possible)
   
2. For parameters that can be obtained without depending on other modules, do not reference parameter values saved in other modules. (Such as `mpi_num`, `start_time`)

3. Use classes as function parameters as much as possible instead of using global classes for obtained parameters. (For example, in `gen_general_info`, `Input`)

4. After adding parameters, supplement test code in `module_io/json_output/test`.

For the current JSON file, there are two JSON modules: `init` and `general_info`,  `output_info`.
Taking `general_info` as an example, the code to be added is as follows:

```cpp
namespace Json
{

#ifdef __RAPIDJSON
void gen_general_info(Input *input)
{

#ifdef VERSION
    const std::string version = VERSION;
#else
    const std::string version = "unknown";
#endif
#ifdef COMMIT
    const std::string commit = COMMIT;
#else
    const std::string commit = "unknown";
#endif

    // start_time
    std::time_t start_time = input->get_start_time();
    std::string start_time_str;
    convert_time(start_time, start_time_str);

    // end_time
    std::time_t time_now = std::time(NULL);
    std::string end_time_str;
    convert_time(time_now, end_time_str);

#ifdef __MPI
    int mpi_num = Parallel_Global::mpi_number;
    int omp_num = Parallel_Global::omp_number;
#elif
    int mpi_num = 1;
    int omp_num = 1;
#endif

    AbacusJson::add_json({"general_info", "version"}, version,false);
    AbacusJson::add_json({"general_info", "commit"}, commit,false);
    AbacusJson::add_json({"general_info", "device"}, input->device,false);
    AbacusJson::add_json({"general_info", "mpi_num"}, mpi_num,false);
    AbacusJson::add_json({"general_info", "omp_num"}, omp_num,false);
    AbacusJson::add_json({"general_info", "pseudo_dir"}, input->pseudo_dir,false);
    AbacusJson::add_json({"general_info", "orbital_dir"}, input->orbital_dir,false);
    AbacusJson::add_json({"general_info", "stru_file"}, input->stru_file,false);
    AbacusJson::add_json({"general_info", "kpt_file"}, input->kpoint_file,false);
    AbacusJson::add_json({"general_info", "start_time"}, start_time_str,false);
    AbacusJson::add_json({"general_info", "end_time"}, end_time_str,false);
}
#endif
} // namespace Json
```