#ifndef ABACUS_JSON_H
#define ABACUS_JSON_H

#include <string>
#include <vector>
#include "json_node.h"

#ifdef __JSON
// Keep the implementation-heavy json.hpp out of this header.
#include <nlohmann/json_fwd.hpp>

namespace Json
{

using jsonValue = nlohmann::ordered_json;

class AbacusJsonTestAccess;

class AbacusJson
{
  public:
    static void write_to_json(const std::string& filename);

    /**
     * Replace a value at a named or indexed path, including whole containers.
     * Missing named parents are created as objects. Integer indices must refer
     * to existing array elements; negative indices count from the end.
     * An empty path leaves the document unchanged.
     */
    static void set_json(const std::vector<jsonKeyNode>& keys, jsonValue value);

    /**
     * Append one value to an array at the path, without flattening that value.
     * A missing named destination is created as an array. An existing
     * destination must be an array, including when selected by an integer
     * index; nulls, objects and scalars are rejected. Path rules match set_json.
     */
    static void append_json(const std::vector<jsonKeyNode>& keys, jsonValue value);

  private:
    friend class AbacusJsonTestAccess;
    static jsonValue doc;
};

} // namespace Json
#endif // __JSON

#endif
