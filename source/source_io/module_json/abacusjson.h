#ifndef ABACUS_JSON_H
#define ABACUS_JSON_H

#include <string>

#ifdef __JSON
// Keep the implementation-heavy json.hpp out of this header.
#include <nlohmann/json_fwd.hpp>

namespace Json
{

using jsonValue = nlohmann::ordered_json;

class AbacusJson
{
  public:
    // Shared document for the schema generators in module_json; keep its root an object.
    static jsonValue& document();
    static void write_to_json(const std::string& filename);

  private:
    static jsonValue doc;
};

} // namespace Json
#endif // __JSON

#endif
