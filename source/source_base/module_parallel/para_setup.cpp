#include "para_setup.h"

#include <cassert>
#include <iostream>

#ifdef __MPI
#include <mpi.h>
#endif

namespace Parallel
{

void divide_mpi_groups(int nproc, int num_groups, int rank, bool even,
                       int& procs_in_group, int& my_group, int& rank_in_group)
{
    assert(num_groups > 0);
    assert(nproc >= num_groups);

    procs_in_group = nproc / num_groups;
    int extra_procs = nproc % num_groups;

    if (even && extra_procs != 0)
    {
        std::cerr << "Error: " << nproc << " processes not evenly divisible by "
                  << num_groups << " groups." << std::endl;
        assert(false);
    }

    if (rank < extra_procs * (procs_in_group + 1))
    {
        procs_in_group++;
        my_group = rank / procs_in_group;
        rank_in_group = rank % procs_in_group;
    }
    else
    {
        my_group = (rank - extra_procs) / procs_in_group;
        rank_in_group = (rank - extra_procs) % procs_in_group;
    }
}

#ifdef __MPI

namespace {

// Helper: split a parent communicator into ngroup sub-communicators.
// Mirrors MPICommGroup::divide_group_comm in parallel_comm.cpp:
//   - group_comm: intra-group communicator (color = my_group)
//   - inter_comm: communicator of same-rank processes across groups
//                 (color = rank_in_group); MPI_COMM_NULL for a single
//                 group or an uneven split, exactly like KP_WORLD.
struct GroupSplitResult
{
    MPI_Comm group_comm = MPI_COMM_NULL;
    MPI_Comm inter_comm = MPI_COMM_NULL;
    int ngroups = 0;
    int nprocs_in_group = 0;
    int my_group = 0;
    int rank_in_group = 0;
};

GroupSplitResult split_comm_group(MPI_Comm parent, int ngroup, bool even)
{
    GroupSplitResult res;
    res.ngroups = ngroup;

    int gsize = 0;
    int grank = 0;
    MPI_Comm_size(parent, &gsize);
    MPI_Comm_rank(parent, &grank);

    divide_mpi_groups(gsize, ngroup, grank, even,
                      res.nprocs_in_group, res.my_group, res.rank_in_group);

    // Intra-group communicator: one sub-communicator per group.
    MPI_Comm_split(parent, res.my_group, res.rank_in_group, &res.group_comm);

    // Inter-group communicator: processes with the same rank inside
    // their group talk to each other. Only valid for an even split;
    // an uneven split leaves some groups without a corresponding rank.
    const bool is_even = (gsize % ngroup == 0);
    if (ngroup > 1 && is_even)
    {
        MPI_Comm_split(parent, res.rank_in_group, res.my_group, &res.inter_comm);
    }

    return res;
}

} // anonymous namespace

void split_images(int nproc, int my_rank, int nimage,
                  int& image_id, int& rank_in_esolver, int& esolver_size,
                  ParaWorld& esolver_world, ParaWorld& images_world)
{
    assert(nimage > 0);
    assert(nproc >= nimage);

    int procs_in_image = 0;
    divide_mpi_groups(nproc, nimage, my_rank, false,
                      procs_in_image, image_id, rank_in_esolver);
    esolver_size = procs_in_image;

    // Intra-image domain: all processes of one esolver.
    MPI_Comm esolver_comm;
    MPI_Comm_split(MPI_COMM_WORLD, image_id, rank_in_esolver, &esolver_comm);
    esolver_world = ParaWorld::make_mpi(ParaTag::esolver, esolver_comm);

    // Inter-image domain: same rank_in_esolver across images. Follows the
    // KP_WORLD convention: absent for a single image or an uneven split.
    const bool is_even = (nproc % nimage == 0);
    if (nimage > 1 && is_even)
    {
        MPI_Comm images_comm;
        MPI_Comm_split(MPI_COMM_WORLD, rank_in_esolver, image_id, &images_comm);
        images_world = ParaWorld::make_mpi(ParaTag::images, images_comm);
    }
    else
    {
        images_world = ParaWorld::make_mpi(ParaTag::images, MPI_COMM_NULL);
    }
}

void split_pools(int parent_size, int parent_rank, int bndpar, int kpar,
                 const MPI_Comm& parent_comm,
                 int& nproc_in_pool, int& rank_in_pool, int& my_pool,
                 int& nproc_in_bndgroup, int& rank_in_bpgroup, int& my_bndgroup,
                 ParaWorld& pw_world, ParaWorld& kmesh_world,
                 ParaWorld& bgroup_int, ParaWorld& bgroup_bp)
{
    if (bndpar > 1 && parent_size % (bndpar * kpar) != 0)
    {
        std::cerr << "Error: " << parent_size
                  << " processes in the parent domain must be divisible by "
                  << "BNDPAR*KPAR (" << bndpar * kpar << ")." << std::endl;
        assert(false);
    }

    // k-point parallelization: split the parent domain into kpar pools.
    GroupSplitResult kpar_res = split_comm_group(parent_comm, kpar, false);

    // band parallelization: split each pool into bndpar groups.
    GroupSplitResult bndpar_res = split_comm_group(kpar_res.group_comm, bndpar, true);

    // Set output indices.
    nproc_in_pool = bndpar_res.nprocs_in_group;
    rank_in_pool = bndpar_res.rank_in_group;
    my_pool = kpar_res.my_group;

    // POOL_WORLD: processes with the same k point and the same bands
    // (plane-wave distribution lives inside it).
    MPI_Comm pool_comm;
    MPI_Comm_dup(bndpar_res.group_comm, &pool_comm);
    pw_world = ParaWorld::make_mpi(ParaTag::pw, pool_comm);

    // KP_WORLD: inter-pool communicator (same rank across pools).
    if (kpar_res.inter_comm != MPI_COMM_NULL)
    {
        MPI_Comm kp_comm;
        MPI_Comm_dup(kpar_res.inter_comm, &kp_comm);
        kmesh_world = ParaWorld::make_mpi(ParaTag::kmesh, kp_comm);
    }
    else
    {
        kmesh_world = ParaWorld::make_mpi(ParaTag::kmesh, MPI_COMM_NULL);
    }

    // Band group communicators.
    if (bndpar > 1)
    {
        nproc_in_bndgroup = kpar_res.ngroups * bndpar_res.nprocs_in_group;
        rank_in_bpgroup = kpar_res.my_group * bndpar_res.nprocs_in_group + bndpar_res.rank_in_group;
        my_bndgroup = bndpar_res.my_group;

        // INT_BGROUP: same bands across pools (bsame_kdiff).
        MPI_Comm int_bgroup;
        MPI_Comm_split(parent_comm, my_bndgroup, rank_in_bpgroup, &int_bgroup);
        bgroup_int = ParaWorld::make_mpi(ParaTag::bsame_kdiff, int_bgroup);

        // BP_WORLD: same k point across band groups (bdiff_ksame).
        MPI_Comm bp_comm;
        MPI_Comm_dup(bndpar_res.inter_comm, &bp_comm);
        bgroup_bp = ParaWorld::make_mpi(ParaTag::bdiff_ksame, bp_comm);
    }
    else
    {
        nproc_in_bndgroup = parent_size;
        rank_in_bpgroup = parent_rank;
        my_bndgroup = 0;

        // No band parallelism: INT_BGROUP spans the whole parent domain,
        // BP_WORLD degenerates to one process per rank.
        MPI_Comm int_bgroup;
        MPI_Comm_dup(parent_comm, &int_bgroup);
        bgroup_int = ParaWorld::make_mpi(ParaTag::bsame_kdiff, int_bgroup);

        MPI_Comm bp_comm;
        MPI_Comm_split(parent_comm, parent_rank, 0, &bp_comm);
        bgroup_bp = ParaWorld::make_mpi(ParaTag::bdiff_ksame, bp_comm);
    }
}

ParaWorld split_diag_world(int diag_np, int parent_size, int parent_rank,
                           const MPI_Comm& parent_comm,
                           int& drank, int& dsize, int& dcolor)
{
    assert(diag_np > 0);

    int procs_in_group = 0;
    int my_group = 0;
    int rank_in_group = 0;
    divide_mpi_groups(parent_size, diag_np, parent_rank, false,
                      procs_in_group, my_group, rank_in_group);

    MPI_Comm diag_comm;
    MPI_Comm_split(parent_comm, my_group, rank_in_group, &diag_comm);

    MPI_Comm_rank(diag_comm, &drank);
    MPI_Comm_size(diag_comm, &dsize);
    dcolor = my_group;

    return ParaWorld::make_mpi(ParaTag::diag, diag_comm);
}

ParaWorld split_grid_world(int diag_np, int parent_size, int parent_rank,
                           const MPI_Comm& parent_comm,
                           int& grank, int& gsize)
{
    assert(diag_np > 0);

    int procs_in_group = 0;
    int my_group = 0;
    int rank_in_group = 0;
    divide_mpi_groups(parent_size, diag_np, parent_rank, false,
                      procs_in_group, my_group, rank_in_group);

    MPI_Comm grid_comm;
    MPI_Comm_split(parent_comm, my_group, rank_in_group, &grid_comm);

    MPI_Comm_rank(grid_comm, &grank);
    MPI_Comm_size(grid_comm, &gsize);

    return ParaWorld::make_mpi(ParaTag::rgrid, grid_comm);
}

ParaCollection setup_para_worlds(int nproc, int my_rank, int nimage,
                                 int bndpar, int kpar, int diag_np)
{
    ParaCollection worlds;

    // 0. Top-level split: independent images.
    //    esolver_world contains all processes of one esolver instance;
    //    images_world connects corresponding ranks across images.
    int image_id = 0;
    int rank_in_esolver = 0;
    int esolver_size = 0;
    ParaWorld esolver_world = ParaWorld::make_mpi(ParaTag::esolver, MPI_COMM_NULL);
    ParaWorld images_world = ParaWorld::make_mpi(ParaTag::images, MPI_COMM_NULL);
    split_images(nproc, my_rank, nimage,
                 image_id, rank_in_esolver, esolver_size,
                 esolver_world, images_world);
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::esolver, esolver_world.comm()));
    // images_world may be an invalid domain (nimage == 1 or uneven split);
    // it is still registered so that find(ParaTag::images) returns it and
    // callers can test valid().
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::images, images_world.comm()));

    // All solver domains are derived from the esolver domain, never from
    // MPI_COMM_WORLD directly (see the hierarchy diagram in para_setup.h).
    const MPI_Comm esolver_comm = esolver_world.comm();

    // 1. k-pools and band groups.
    int nproc_in_pool = 0;
    int rank_in_pool = 0;
    int my_pool = 0;
    int nproc_in_bndgroup = 0;
    int rank_in_bpgroup = 0;
    int my_bndgroup = 0;

    ParaWorld pw_world = ParaWorld::make_mpi(ParaTag::pw, MPI_COMM_NULL);
    ParaWorld kmesh_world = ParaWorld::make_mpi(ParaTag::kmesh, MPI_COMM_NULL);
    ParaWorld bgroup_int = ParaWorld::make_mpi(ParaTag::bsame_kdiff, MPI_COMM_NULL);
    ParaWorld bgroup_bp = ParaWorld::make_mpi(ParaTag::bdiff_ksame, MPI_COMM_NULL);

    split_pools(esolver_size, rank_in_esolver, bndpar, kpar, esolver_comm,
                nproc_in_pool, rank_in_pool, my_pool,
                nproc_in_bndgroup, rank_in_bpgroup, my_bndgroup,
                pw_world, kmesh_world, bgroup_int, bgroup_bp);

    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::pw, pw_world.comm()));
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::kmesh, kmesh_world.comm()));
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::bsame_kdiff, bgroup_int.comm()));
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::bdiff_ksame, bgroup_bp.comm()));

    // 2. Diagonalization domain.
    int drank = 0;
    int dsize = 0;
    int dcolor = 0;
    ParaWorld diag_world = split_diag_world(diag_np, esolver_size, rank_in_esolver,
                                            esolver_comm, drank, dsize, dcolor);
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::diag, diag_world.comm()));

    // 3. Real-space grid domain.
    int grank = 0;
    int gsize = 0;
    ParaWorld grid_world = split_grid_world(diag_np, esolver_size, rank_in_esolver,
                                            esolver_comm, grank, gsize);
    worlds.add(ParaWorld::make_mpi_ptr(ParaTag::rgrid, grid_world.comm()));

    // 4. Matrix domain: serial for now until its own 2D-grid split lands.
    worlds.add(ParaWorld::make_serial(ParaTag::matrix));

    return worlds;
}

#endif // __MPI

} // namespace Parallel
