#include "abacusjson.h"

#ifdef __JSON
#include <nlohmann/json.hpp>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <utility>

namespace Json
{
namespace
{
// Only missing named nodes are created. Indexed access never grows an array.
jsonValue* resolve_path(jsonValue& root,
                        const std::vector<jsonKeyNode>& keys,
                        jsonValue initial_value)
{
    if (keys.empty())
    {
        return nullptr;
    }

    jsonValue* parent = &root;
    for (std::size_t i = 0; i < keys.size(); ++i)
    {
        const jsonKeyNode& key = keys[i];
        if (key.is_index)
        {
            if (!parent->is_array())
            {
                throw std::invalid_argument("JSON output: an integer path component requires an array");
            }
            const std::ptrdiff_t size = static_cast<std::ptrdiff_t>(parent->size());
            std::ptrdiff_t index = static_cast<std::ptrdiff_t>(key.i);
            if (index < 0)
            {
                index += size;
            }
            if (index < 0 || index >= size)
            {
                throw std::out_of_range("JSON output: array index out of range");
            }
            parent = &parent->at(static_cast<jsonValue::size_type>(index));
        }
        else
        {
            if (!parent->is_object())
            {
                throw std::invalid_argument("JSON output: a named path component requires an object");
            }
            jsonValue::iterator child = parent->find(key.key);
            if (child == parent->end())
            {
                jsonValue initial = i + 1 == keys.size() ? std::move(initial_value) : jsonValue::object();
                child = parent->emplace(key.key, std::move(initial)).first;
            }
            parent = &child.value();
        }
    }
    return parent;
}
} // namespace

jsonValue AbacusJson::doc = jsonValue::object();

void AbacusJson::set_json(const std::vector<jsonKeyNode>& keys, jsonValue value)
{
    jsonValue* target = resolve_path(doc, keys, nullptr);
    if (target != nullptr)
    {
        *target = std::move(value);
    }
}

void AbacusJson::append_json(const std::vector<jsonKeyNode>& keys, jsonValue value)
{
    jsonValue* target = resolve_path(doc, keys, jsonValue::array());
    if (target == nullptr)
    {
        return;
    }
    if (!target->is_array())
    {
        throw std::invalid_argument("JSON output: append requires an array");
    }
    target->push_back(std::move(value));
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
