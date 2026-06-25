#ifndef GRID_H
#define GRID_H

#include "atom_pack.h"
#include "source_cell/unitcell.h"
#include "sltk_atom.h"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

typedef std::vector<FAtom> AtomMap;

class Grid
{
  public:
    enum class NeighborBuildMode
    {
        // Production route: build only AtomPack/GridStorage neighbor data.
        // The legacy FAtom* containers stay empty in this mode.
        AtomPackOnly,
        // Regression route: also build atoms_in_box/all_adj_info so tests can
        // compare the new integer-indexed path with the original Grid output.
        AtomPackAndLegacy
    };

    enum class NeighborSearchMode
    {
        // Visit the complete 3x3x3 box stencil for every center atom.
        Full27,
        // Visit one half-domain box stencil and restore directed symmetry.
        Half14
    };

    enum class NeighborReferenceMode
    {
        // Production route: do not build the full 27-direction reference list.
        None,
        // Regression route: also build neighbor_pairs_27 for correctness checks.
        Full27
    };

    // Constructors and destructor
    // Grid is Global class,so init it with constant number
    Grid() : test_grid(0){};
    Grid(const int& test_grid_in);
    virtual ~Grid();

    Grid& operator=(Grid&&) = default;

    // Build a replicated Grid. radius_in and all stored coordinates use lat0
    // units. The default path builds only the Half14 AtomPack representation.
    void init(std::ofstream& ofs,
              const UnitCell& ucell,
              const double radius_in,
              const bool boundary = true,
              const NeighborBuildMode build_mode = NeighborBuildMode::AtomPackOnly,
              const NeighborSearchMode search_mode = NeighborSearchMode::Half14,
              const NeighborReferenceMode reference_mode = NeighborReferenceMode::None);
    // Build a spatially distributed Grid containing only locally owned centers
    // plus received ghost candidates. Every rank must call this collectively
    // with the same cell, radius, boundary mode and communicator.
    void init_mpi(std::ofstream& ofs,
                  const UnitCell& ucell,
                  const double radius_in,
                  const bool boundary,
                  ModuleNeighbor::NeighborMpiComm communicator,
                  ModuleNeighbor::MpiGhostExchangeStats* stats = nullptr);

    bool mpi_mode() const
    {
        return mpi_mode_;
    }
    // In MPI mode, return whether this rank owns the requested physical atom.
    // Replicated grids treat every valid center as local.
    bool is_local_center(const int type, const int natom) const;
    const ModuleNeighbor::MpiDomain& mpi_domain() const
    {
        return mpi_domain_;
    }

    // Data
    bool pbc=false; // When pbc is set to false, periodic boundary conditions are explicitly ignored.
    double sradius2=0.0; // searching radius squared (unit:lat0)
    double sradius=0.0;  // searching radius (unit:lat0)
    // The Grid may be built with an additional Verlet skin. query_radius is
    // the physical cutoff applied when cached candidates are returned.
    double query_radius2=0.0;
    double query_radius=0.0;

    // Select the physical cutoff used to filter a candidate list built at
    // sradius. The value must be finite and lie in [0, sradius].
    void set_query_radius(const double radius)
    {
        if (!std::isfinite(radius) || radius < 0.0 || radius > sradius)
        {
            throw std::invalid_argument("Grid query radius must be within the build radius.");
        }
        query_radius = radius;
        query_radius2 = radius * radius;
    }
    
    // coordinate range of the input atom (unit:lat0)
    double x_min=0.0;
    double y_min=0.0;
    double z_min=0.0;
    double x_max=0.0;
    double y_max=0.0;
    double z_max=0.0;

    // The algorithm for searching neighboring atoms uses a "box" partitioning method. 
    // Each box has an edge length of sradius + 0.1, and the number of boxes in each direction is recorded here.
    double box_edge_length=0.0;
    int box_nx=0;
    int box_ny=0;
    int box_nz=0;

    void getBox(int& bx, int& by, int& bz, const double& x, const double& y, const double& z)
    {
        bx = std::floor((x - x_min) / box_edge_length);
        by = std::floor((y - y_min) / box_edge_length);
        bz = std::floor((z - z_min) / box_edge_length);
    }
    // Stores the atoms after box partitioning.
    std::vector<std::vector<std::vector<AtomMap>>> atoms_in_box;

    // Stores the adjacent information of atoms. [ntype][natom][adj list]
    std::vector<std::vector< std::vector<FAtom *> >> all_adj_info;

    // Flat integer-indexed search data. neighbor_pairs is the production result
    // queried by Grid_Driver; neighbor_pairs_27 is an optional test reference.
    ModuleNeighbor::AtomPack atom_pack;
    ModuleNeighbor::GridStorage grid_storage;
    std::vector<ModuleNeighbor::NeighborPair> neighbor_pairs;
    std::vector<ModuleNeighbor::NeighborPair> neighbor_pairs_27;
    ModuleNeighbor::PagedNeighborList paged_neighbor_list;
    // MPI ownership mask indexed by [type][natom]. It remains empty for a
    // replicated Grid, where every center atom is queryable.
    std::vector<std::vector<bool>> local_center_mask;

    void clear_atoms()
    {
        // we have to clear the all_adj_info
        // because the pointers point to the memory in vector atoms_in_box
        all_adj_info.clear();

        atoms_in_box.clear();
        atom_pack.clear();
        grid_storage.clear();
        neighbor_pairs.clear();
        neighbor_pairs_27.clear();
        paged_neighbor_list.clear();
        local_center_mask.clear();
        mpi_mode_ = false;
        mpi_domain_ = ModuleNeighbor::MpiDomain();
    }
    void clear_adj_info()
    {
        // here dont need to free the memory, 
        // because the pointers point to the memory in vector atoms_in_box
        all_adj_info.clear();
    }
    int getGlayerX() const
    {
        return glayerX;
    }
    int getGlayerY() const
    {
        return glayerY;
    }
    int getGlayerZ() const
    {
        return glayerZ;
    }
    int getGlayerX_minus() const
    {
        return glayerX_minus;
    }
    int getGlayerY_minus() const
    {
        return glayerY_minus;
    }
    int getGlayerZ_minus() const
    {
        return glayerZ_minus;
    }
  private:
    int test_grid;

    void setMemberVariables(std::ofstream& ofs_in, const UnitCell& ucell);

    void Construct_Adjacent(const UnitCell& ucell);
    void Construct_Adjacent_near_box(const FAtom& fatom);
    void Construct_Adjacent_final(const FAtom& fatom1, FAtom* fatom2);

    // Build the selected integer-indexed pair list and an optional Full27
    // correctness reference from the current AtomPack/GridStorage.
    void Build_AtomPack_Search_Path(const UnitCell& ucell,
                                    const NeighborSearchMode search_mode,
                                    const NeighborReferenceMode reference_mode);

    // Convert neighbor_pairs into per-center fixed-size page chains.
    void Build_Paged_Neighbor_List(const UnitCell& ucell);

    // Build the original pointer-based cell list only when explicitly requested.
    void Build_Legacy_Search_Path(const UnitCell& ucell);

    // Expand received MPI records by the periodic images required for distance
    // checks while preserving their physical atom identity and ghost status.
    void Append_Periodic_Images(const UnitCell& ucell,
                                const std::vector<ModuleNeighbor::MpiAtomRecord>& records);

    void Check_Expand_Condition(const UnitCell& ucell);
    bool mpi_mode_ = false;
    ModuleNeighbor::MpiDomain mpi_domain_;
    int glayerX=0;
    int glayerX_minus=0;
    int glayerY=0;
    int glayerY_minus=0;
    int glayerZ=0;
    int glayerZ_minus=0;
};

#endif
