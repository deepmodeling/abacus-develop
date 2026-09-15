#ifndef PARA_SETUP_H
#define PARA_SETUP_H

#include "para_collection.h"
#include "para_tag.h"
#include "para_world.h"

#ifdef __MPI
#include "mpi.h"
#endif

namespace Parallel
{

/**
 * @file para_setup.h
 *
 * @brief Construction of the parallel communication domain hierarchy.
 *
 * The domains form a tree rooted at MPI_COMM_WORLD. The first split
 * separates independent calculation images (e.g. NEB replicas, each with
 * its own unit cell); every existing solver domain is then derived from
 * the intra-image domain instead of being hard-wired to MPI_COMM_WORLD:
 *
 * @code
 * MPI_COMM_WORLD
 * |
 * +-- split by nimage (same MPICommGroup pattern as k-parallelism)
 * |
 * |   para_esolver_world   [color = image_id]
 * |   |                    all processes belonging to one esolver
 * |   |
 * |   +-- all existing domains are derived from it
 * |   |       (instead of being hard-wired to MPI_COMM_WORLD):
 * |   |     +-- kmesh / pw               (split by kpar)
 * |   |     +-- bsame_kdiff / bdiff_ksame (split by bndpar)
 * |   |     +-- rgrid / diag             (split by diago_proc)
 * |   |     +-- matrix
 * |   |
 * |   para_images_world   [color = rank_in_esolver]
 * |                        cross-image: processes with the same
 * |                        rank inside their esolver
 * @endcode
 *
 * Degenerate cases follow the existing KP_WORLD convention:
 *   - nimage == 1: esolver_world is a dup of MPI_COMM_WORLD and
 *     images_world is MPI_COMM_NULL (invalid domain).
 *   - Uneven image split: images_world is MPI_COMM_NULL as well,
 *     because corresponding ranks do not exist across all images.
 */

/**
 * @brief Divide nproc processes into num_groups groups.
 *
 * Replaces Parallel_Global::divide_mpi_groups. Works in both serial and
 * MPI builds.
 *
 * @param[in] nproc        total number of processes
 * @param[in] num_groups   desired number of groups
 * @param[in] rank         global rank
 * @param[in] even         if true, require an even split (assert on failure)
 * @param[out] procs_in_group  processes per group
 * @param[out] my_group        which group this rank belongs to
 * @param[out] rank_in_group   rank within the group
 */
void divide_mpi_groups(int nproc, int num_groups, int rank, bool even,
                       int& procs_in_group, int& my_group, int& rank_in_group);

#ifdef __MPI

/**
 * @brief Split MPI_COMM_WORLD into nimage independent images.
 *
 * Produces the two top-level domains:
 *   - esolver_world: intra-image communicator (all processes of one esolver),
 *     split with color = image_id;
 *   - images_world: inter-image communicator (same rank_in_esolver across
 *     images), split with color = rank_in_esolver; MPI_COMM_NULL when
 *     nimage == 1 or the split is uneven.
 *
 * @param[in]  nproc               total number of MPI processes
 * @param[in]  my_rank             global rank
 * @param[in]  nimage              number of images (>= 1)
 * @param[out] image_id            image this rank belongs to
 * @param[out] rank_in_esolver     rank inside the esolver domain
 * @param[out] esolver_size        number of processes in the esolver domain
 * @param[out] esolver_world       intra-image domain (tag: ParaTag::esolver)
 * @param[out] images_world        inter-image domain (tag: ParaTag::images,
 *                                 invalid when nimage == 1 or uneven split)
 */
void split_images(int nproc, int my_rank, int nimage,
                  int& image_id, int& rank_in_esolver, int& esolver_size,
                  ParaWorld& esolver_world, ParaWorld& images_world);

/**
 * @brief Split a parent domain into k-pools and band groups.
 *
 * Replaces Parallel_Global::divide_pools. All splits are performed inside
 * @p parent_comm (normally the esolver domain), so the same code works
 * for both single-image and multi-image runs.
 *
 * @param[in]  parent_size     number of processes in the parent domain
 * @param[in]  parent_rank     rank of this process in the parent domain
 * @param[in]  bndpar          number of band groups
 * @param[in]  kpar            number of k-pools
 * @param[in]  parent_comm     communicator the pools are split from
 * @param[out] nproc_in_pool      processes per pool
 * @param[out] rank_in_pool       rank within pool
 * @param[out] my_pool            pool index
 * @param[out] nproc_in_bndgroup  processes per band group
 * @param[out] rank_in_bpgroup    rank within the band group
 * @param[out] my_bndgroup        band group index
 * @param[out] pw_world           domain for POOL_WORLD
 * @param[out] kmesh_world        domain for KP_WORLD (invalid when kpar == 1)
 * @param[out] bgroup_int         domain for INT_BGROUP
 * @param[out] bgroup_bp          domain for BP_WORLD
 */
void split_pools(int parent_size, int parent_rank, int bndpar, int kpar,
                 const MPI_Comm& parent_comm,
                 int& nproc_in_pool, int& rank_in_pool, int& my_pool,
                 int& nproc_in_bndgroup, int& rank_in_bpgroup, int& my_bndgroup,
                 ParaWorld& pw_world, ParaWorld& kmesh_world,
                 ParaWorld& bgroup_int, ParaWorld& bgroup_bp);

/**
 * @brief Split a parent domain for diagonalization.
 *
 * Replaces Parallel_Global::split_diag_world.
 *
 * @param[in]  diag_np     number of diag groups
 * @param[in]  parent_size number of processes in the parent domain
 * @param[in]  parent_rank rank of this process in the parent domain
 * @param[in]  parent_comm communicator the diag domain is split from
 * @param[out] drank       rank in the diag domain
 * @param[out] dsize       size of the diag domain
 * @param[out] dcolor      color (diag group index)
 * @return  ParaWorld for the diag domain
 */
ParaWorld split_diag_world(int diag_np, int parent_size, int parent_rank,
                           const MPI_Comm& parent_comm,
                           int& drank, int& dsize, int& dcolor);

/**
 * @brief Split a parent domain for the real-space grid.
 *
 * Replaces Parallel_Global::split_grid_world.
 *
 * @param[in]  diag_np     number of grid groups (same parameter as diag)
 * @param[in]  parent_size number of processes in the parent domain
 * @param[in]  parent_rank rank of this process in the parent domain
 * @param[in]  parent_comm communicator the grid domain is split from
 * @param[out] grank       rank in the grid domain
 * @param[out] gsize       size of the grid domain
 * @return  ParaWorld for the rgrid domain
 */
ParaWorld split_grid_world(int diag_np, int parent_size, int parent_rank,
                           const MPI_Comm& parent_comm,
                           int& grank, int& gsize);

/**
 * @brief Assemble the full parallel domain hierarchy.
 *
 * Top-level initialization: splits images first, then derives every
 * solver domain from the esolver domain (see the tree diagram in the
 * file header). Replaces the old divide_pools + split_diag_world +
 * split_grid_world sequence in driver.cpp.
 *
 * @param[in] nproc     total MPI processes (GlobalV::NPROC)
 * @param[in] my_rank   global rank (GlobalV::MY_RANK)
 * @param[in] nimage    number of independent images (1 for a normal run)
 * @param[in] bndpar    number of band groups
 * @param[in] kpar      number of k-pools
 * @param[in] diag_np   number of diag/grid groups
 * @return  ParaCollection containing all domains
 */
ParaCollection setup_para_worlds(int nproc, int my_rank, int nimage,
                                 int bndpar, int kpar, int diag_np);

#endif // __MPI

} // namespace Parallel

#endif // PARA_SETUP_H
