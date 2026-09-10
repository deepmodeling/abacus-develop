#include "para_worlds_global.h"

#include <memory>

#include "source_base/module_parallel/para_collection.h"
#include "source_base/module_parallel/para_setup.h"
#include "source_base/module_parallel/para_tag.h"
#include "source_base/tool_quit.h"

namespace Parallel
{

namespace
{
/// Owning storage for the process-wide collection, plus an initialization
/// latch. Initialized once at startup, then read-only, so this does not
/// introduce mutable cross-module workflow state.
std::unique_ptr<ParaCollection> g_para_worlds;
bool g_initialized = false;
} // namespace

ParaCollection& init_global_para_worlds(int nproc, int my_rank, int nimage)
{
    if (g_initialized)
    {
        ModuleBase::WARNING_QUIT("init_global_para_worlds",
                                 "global ParaCollection is already initialized");
    }

    auto collection = std::unique_ptr<ParaCollection>(new ParaCollection());

    // Argument validation is independent of the MPI build.
    if (nimage < 1 || nproc < nimage)
    {
        ModuleBase::WARNING_QUIT("init_global_para_worlds",
                                 "require 1 <= nimage <= nproc");
    }

#ifdef __MPI
    int image_id = 0;
    int rank_in_esolver = 0;
    int esolver_size = 0;
    ParaWorld esolver_world = ParaWorld::serial(ParaTag::esolver);
    ParaWorld images_world = ParaWorld::serial(ParaTag::images);
    split_images(nproc, my_rank, nimage, image_id, rank_in_esolver,
                 esolver_size, esolver_world, images_world);

    collection->add(std::unique_ptr<ParaWorld>(new ParaWorld(esolver_world)));
    collection->add(std::unique_ptr<ParaWorld>(new ParaWorld(images_world)));
#else
    (void)nproc;
    (void)my_rank;
    (void)nimage;
    collection->add(ParaWorld::make_serial(ParaTag::esolver));
    collection->add(ParaWorld::make_serial(ParaTag::images));
#endif

    g_para_worlds = std::move(collection);
    g_initialized = true;
    return *g_para_worlds;
}

const ParaCollection& global_para_worlds()
{
    if (!g_initialized)
    {
        ModuleBase::WARNING_QUIT("global_para_worlds",
                                 "init_global_para_worlds must be called first");
    }
    return *g_para_worlds;
}

void reset_global_para_worlds_for_test()
{
    g_para_worlds.reset();
    g_initialized = false;
}

} // namespace Parallel
