#ifndef PARA_JSON_H
#define PARA_JSON_H

#include <ctime>
#include <string>

class Parameter;
class UnitCell;
struct Input_para;

namespace Json
{

void create_Json(UnitCell* ucell, const Parameter& param);
void json_output();
void convert_time(std::time_t time_now, std::string& time_str);
void gen_stru_wrapper(UnitCell* ucell, const Input_para& inp);

} // namespace Json

#endif
