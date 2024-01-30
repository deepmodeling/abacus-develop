#include <ctime>
#include <string>

namespace Json
{

// Output the json to abacus.json file
void json_output();

// Convert time_t to string
void convert_time(std::time_t time_now, std::string& time_str);

} // namespace Json
