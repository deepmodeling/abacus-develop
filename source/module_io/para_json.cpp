
#include "para_json.h"

#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

#include "module_base/global_variable.h"
#ifdef __RAPIDJSON
#include "json_output/abacusjson.h"
#include "json_output/general_info.h"
#include "module_io/input.h"
#endif // __RAPIDJSON

namespace Json
{
void json_output()
{
#ifdef __RAPIDJSON
    gen_general_info();
#ifdef __MPI
    if (GlobalV::MY_RANK == 0)
        AbacusJson::write_to_json("abacus.json");
#elif
    AbacusJson::write_to_json("abacus.json");
#endif

#endif // __RAPIDJSON
}

void convert_time(std::time_t time_now, std::string& time_str)
{
    std::tm* tm = std::localtime(&time_now);
    std::ostringstream oss;
    oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
    time_str = oss.str();
}

} // namespace Json
