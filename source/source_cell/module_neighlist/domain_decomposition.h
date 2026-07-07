#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#ifdef __MPI

#include "source_cell/module_neighlist/atom_provider.h"
#include "source_cell/module_neighlist/local_atom.h"

#include <array>
#include <vector>

#include <mpi.h>

/**
 * @brief MPI domain decomposition for distributed neighbor-search input.
 *
 * The decomposition is performed in fractional coordinates. Owned atoms are
 * selected by wrapped fractional position, and ghost atoms are exchanged as
 * shifted periodic images.
 */
class DomainDecomposition
{
public:
    DomainDecomposition();
    ~DomainDecomposition();

    void init(MPI_Comm comm,
              const ModuleBase::Matrix3& latvec,
              double lat0,
              double cutoff,
              double skin);

    int owner_rank_from_frac(const ModuleBase::Vector3<double>& frac) const;

    void split_owned_atoms_from_ucell(const AtomProvider& ucell,
                                      std::vector<LocalAtom>& owned_atoms) const;

    void exchange_ghost_atoms(const std::vector<LocalAtom>& owned_atoms,
                              std::vector<LocalAtom>& ghost_atoms) const;

    const std::array<int, 3>& dims() const;
    const std::array<int, 3>& coords() const;
    int rank() const;
    int size() const;

private:
    struct PackedAtom
    {
        double frac[3];
        int image_shift[3];
        int type;
        int type_index;
        ModuleNeighList::GlobalAtomId global_id;
        int owner_rank;
    };

    struct GhostExchangeSlot
    {
        std::array<int, 3> offset;
        std::array<int, 3> target_coords;
        std::array<int, 3> image_shift;
        std::array<int, 3> recv_image_shift;
        int send_rank;
        int recv_rank;
    };

    MPI_Comm comm_;
    MPI_Comm cart_comm_;
    bool owns_cart_comm_;
    int rank_;
    int size_;
    std::array<int, 3> dims_;
    std::array<int, 3> coords_;
    std::array<double, 3> margin_;
    ModuleBase::Matrix3 latvec_;
    ModuleBase::Matrix3 inv_latvec_;
    double lat0_;
    double cutoff_;
    double skin_;

    static double wrap_fractional(double value);
    static int floor_div(int value, int divisor);
    static int positive_mod(int value, int divisor);
    static double dot_product(const ModuleBase::Vector3<double>& a,
                              const ModuleBase::Vector3<double>& b);
    static ModuleBase::Vector3<double> cross_product(const ModuleBase::Vector3<double>& a,
                                                     const ModuleBase::Vector3<double>& b);
    static double norm(const ModuleBase::Vector3<double>& value);

    ModuleBase::Vector3<double> wrapped_frac_from_cart(const ModuleBase::Vector3<double>& cart) const;
    int rank_from_coords(const std::array<int, 3>& coords) const;
    void target_for_offset(const std::array<int, 3>& offset,
                           std::array<int, 3>& target_coords,
                           std::array<int, 3>& image_shift) const;
    bool atom_overlaps_target_halo(const LocalAtom& atom,
                                   const std::array<int, 3>& target_coords,
                                   const std::array<int, 3>& image_shift) const;
    int neighbor_layer(int dim) const;
    void build_ghost_exchange_slots(std::vector<GhostExchangeSlot>& slots) const;
    PackedAtom pack_atom(const LocalAtom& atom, const std::array<int, 3>& image_shift) const;
    LocalAtom unpack_ghost_atom(const PackedAtom& packed) const;
};

#endif // __MPI

#endif // DOMAIN_DECOMPOSITION_H
