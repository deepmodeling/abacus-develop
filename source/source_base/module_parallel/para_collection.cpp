#include "para_collection.h"

namespace Parallel
{

void ParaCollection::add(std::unique_ptr<ParaWorld> world)
{
    for (const auto& existing : worlds_)
    {
        if (existing->tag() == world->tag())
        {
            return;
        }
    }
    worlds_.push_back(std::move(world));
}

const ParaWorld& ParaCollection::find(const std::string& tag) const
{
    for (const auto& world : worlds_)
    {
        if (world->tag() == tag)
        {
            return *world;
        }
    }
    static const ParaWorld empty = ParaWorld::serial("");
    return empty;
}

} // namespace Parallel
