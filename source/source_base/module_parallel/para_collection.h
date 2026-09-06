#ifndef PARA_COLLECTION_H
#define PARA_COLLECTION_H

#include <memory>
#include <string>
#include <vector>

#include "para_world.h"

namespace Parallel
{

/**
 * @brief Container for all parallel communication domains.
 *
 * A ParaCollection owns a set of ParaWorld objects (base class pointers),
 * each describing one communication domain (see ParaTag). Callers look up
 * domains by tag via find(); a missing tag yields a static empty (invalid)
 * domain as a safe degradation, never an exception.
 *
 * The collection is passed explicitly to functions that need communicator
 * access, replacing reads of loose globals such as GlobalV::POOL_WORLD.
 */
class ParaCollection
{
public:
    ParaCollection() = default;

    /**
     * @brief Append a domain to the collection.
     *
     * Duplicate tags are rejected (the existing entry is kept).
     *
     * @param[in] world  domain to add (ownership transferred)
     */
    void add(std::unique_ptr<ParaWorld> world);

    /**
     * @brief Look up a domain by tag.
     *
     * @param[in] tag  domain tag string
     * @return the matching ParaWorld, or a static empty domain if not found
     */
    const ParaWorld& find(const std::string& tag) const;

    /**
     * @brief Look up a domain by tag and cast to the requested subclass.
     *
     * @tparam T   expected subclass (e.g. ParaKmeshWorld)
     * @param[in] tag  domain tag string
     * @return pointer to the domain if found and type matches, nullptr otherwise
     */
    template <typename T>
    const T* find_as(const std::string& tag) const;

    /**
     * @brief Number of domains in the collection.
     */
    size_t size() const
    {
        return worlds_.size();
    }

private:
    std::vector<std::unique_ptr<ParaWorld>> worlds_; ///< owned domains
};

template <typename T>
const T* ParaCollection::find_as(const std::string& tag) const
{
    for (const auto& world : worlds_)
    {
        if (world->tag() == tag)
        {
            return dynamic_cast<const T*>(world.get());
        }
    }
    return nullptr;
}

} // namespace Parallel

#endif // PARA_COLLECTION_H
