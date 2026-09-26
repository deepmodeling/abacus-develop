#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_cell/module_neighlist/local_atom.h"

#include <array>
#include <cstdint>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

class UnitCell;
class MDCell;
namespace ModuleBase
{
class CommunicationDomain;
}

/**
 * @brief Domain decomposition for distributed and single-domain neighbor input.
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
    DomainDecomposition(const DomainDecomposition&) = delete;
    DomainDecomposition& operator=(const DomainDecomposition&) = delete;

    void init(const ModuleBase::CommunicationDomain& comm_domain,
              const ModuleBase::Matrix3& latvec,
              double lat0,
              double cutoff,
              double skin);

    std::vector<LocalAtom> split_owned_atoms_from_ucell(const UnitCell& ucell) const;
    void migrate_owned_atoms(MDCell& cell);
    void exchange_ghost_atoms(MDCell& cell);
    void update_ghost_atom_positions(MDCell& cell);
    void accumulate_ghost_forces(MDCell& cell);
    void prepare_neighbors(MDCell& cell);

    int owner_rank_from_frac(const ModuleBase::Vector3<double>& frac) const;

    const std::array<int, 3>& dims() const;
    const std::array<int, 3>& coords() const;

private:
    struct PackedAtom
    {
        double frac[3];
        double vel[3];
        double force[3];
        int mbl[3];
        double mass;
        int image_shift[3];
        int type;
        std::int64_t type_index;
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
        std::vector<int> send_atom_indices;
        std::size_t ghost_begin;
        int ghost_count;
    };

    struct ForceRecord
    {
        int type;
        std::int64_t type_index;
        double force[3];
    };

    // Geometry is cached here; the physical state and neighbor data belong to MDCell.
    void update_geometry_(const ModuleBase::Matrix3& latvec, double lat0, double cutoff, double skin);
    void synchronize_geometry_(MDCell& cell);

    // Transport wrappers: only these helpers and communicator lifetime handling
    // depend on MPI. Atom ownership and halo algorithms are shared by all builds.
    void exchange_bytes_(const void* send, std::size_t send_bytes,
                         void* recv, std::size_t recv_bytes,
                         int send_rank, int recv_rank, int tag) const;
    void exchange_migration_(const std::array<int, 6>& neighbors,
                             const std::array<std::vector<PackedAtom>, 6>& send_atoms,
                             const std::array<int, 6>& send_counts,
                             std::array<std::vector<PackedAtom>, 6>& recv_atoms) const;

    static double wrap_fractional(double value);
    static int floor_div(int value, int divisor);
    static int positive_mod(int value, int divisor);
    ModuleBase::Vector3<double> wrapped_frac_from_cart(const ModuleBase::Vector3<double>& cart) const;
    int rank_from_coords(const std::array<int, 3>& coords) const;
    int neighbor_layer(int dim) const;
    void target_for_offset(const std::array<int, 3>& offset,
                           std::array<int, 3>& target_coords,
                           std::array<int, 3>& image_shift) const;
    void build_ghost_exchange_slots(std::vector<GhostExchangeSlot>& slots) const;
    PackedAtom pack_atom(const LocalAtom& atom, const std::array<int, 3>& image_shift) const;
    LocalAtom unpack_ghost_atom(const PackedAtom& packed) const;
    LocalAtom unpack_owned_atom(const PackedAtom& packed) const;

    bool initialized_ = false;
#ifdef __MPI
    MPI_Comm comm_ = MPI_COMM_NULL;
    MPI_Comm cart_comm_ = MPI_COMM_NULL;
    bool owns_cart_comm_ = false;
#endif
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
    std::vector<GhostExchangeSlot> ghost_slots_;
    bool ghost_layout_valid_ = false;
};

#endif // DOMAIN_DECOMPOSITION_H
