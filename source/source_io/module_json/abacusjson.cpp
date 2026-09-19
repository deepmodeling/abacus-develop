#include "abacusjson.h"

#ifdef __JSON
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>

namespace Json
{

jsonValue AbacusJson::doc = jsonValue::object();

jsonValue& AbacusJson::document()
{
    return doc;
}

void AbacusJson::write_to_json(const std::string& filename)
{
    const auto content = doc.dump(4);
    std::ofstream file(filename);
    if (!file)
    {
        throw std::runtime_error("Cannot open JSON output file: " + filename);
    }
    file << content;
    file.close();
    if (!file)
    {
        throw std::runtime_error("Cannot write JSON output file: " + filename);
    }
}

} // namespace Json
#endif // __JSON
