#ifndef MD_CELL_H
#define MD_CELL_H

#include "source_cell/base_cell.h"
#include "source_cell/md_stru_file_metadata.h"
#include "source_cell/module_neighlist/local_atom.h"

#ifdef __MPI
#include "source_cell/module_neighlist/domain_decomposition.h"
#endif

#include <string>
#include <cstdint>
#include <memory>
#include <vector>

class UnitCell;
class NeighborSearch;
namespace ModuleBase
{
class CommunicationDomain;
}

class MDCell : public BaseCell
{
public:
    MDCell();
    ~MDCell();
    MDCell(const MDCell&) = delete;
    MDCell& operator=(const MDCell&) = delete;
    MDCell(MDCell&&);
    MDCell& operator=(MDCell&&);

    void initialize_from_unitcell(UnitCell& ucell,
                                  double cutoff,
                                  double skin,
                                  const ModuleBase::CommunicationDomain& comm_domain);
    void initialize_from_owned_atoms(const ModuleBase::Matrix3& latvec,
                                     const ModuleBase::Matrix3& gt,
                                     double lat0,
                                     double omega,
                                     std::int64_t nat,
                                     const std::vector<LocalAtom>& owned_atoms,
                                     const std::vector<std::string>& type_labels,
                                     const std::vector<double>& type_masses,
                                     const std::vector<std::int64_t>& type_atom_counts,
                                     double cutoff,
                                     double skin,
                                     const ModuleBase::CommunicationDomain& comm_domain);

#ifdef __MPI
    int mpi_rank() const;
    int mpi_size() const;
    MPI_Comm communicator() const { return comm_; }

#endif

    void exchange_ghost_atoms();
    void accumulate_ghost_forces();
    void migrate_owned_atoms();
    void prepare_neighbors();
    bool has_neighbor_search() const;
    const NeighborSearch& neighbor_search() const;
    void set_lattice_vectors(const ModuleBase::Matrix3& latvec);
    void refresh_cart_from_frac();

    const std::vector<LocalAtom>& owned_atoms() const { return owned_atoms_; }
    const std::vector<LocalAtom>& ghost_atoms() const;
    const std::vector<std::string>& type_labels() const { return type_labels_; }
    const std::vector<double>& type_masses() const { return type_masses_; }
    const std::vector<std::int64_t>& type_atom_counts() const { return type_atom_counts_; }
    MdStruFileMetadata& mutable_stru_file_metadata() { return stru_file_metadata_; }
    const MdStruFileMetadata& stru_file_metadata() const { return stru_file_metadata_; }
    std::vector<LocalAtom>& mutable_owned_atoms();
    std::vector<LocalAtom>& mutable_ghost_atoms();

    int nlocal() const { return static_cast<int>(owned_atoms_.size()); }
    int nghost() const { return static_cast<int>(ghost_atoms_.size()); }
    double cutoff() const;
    bool has_backing_unitcell() const;
    UnitCell& backing_unitcell();
    const UnitCell& backing_unitcell() const;
    void sync_backing_unitcell();

private:
    Kind get_kind() const override;
    std::int64_t get_nat() const override;
    double get_lat0() const override;
    double get_omega() const override;
    const ModuleBase::Matrix3& get_latvec() const override;
    const ModuleBase::Matrix3& get_GT() const override;

#ifdef __MPI
    void initialize_from_ucell_(UnitCell& ucell, MPI_Comm comm, double cutoff, double skin);
    void initialize_from_owned_atoms_(MPI_Comm comm, double cutoff, double skin);
#else
    void initialize_from_ucell_(UnitCell& ucell, double cutoff, double skin);
    void initialize_from_owned_atoms_(double cutoff, double skin);
#endif

    void sync_backing_unitcell_geometry_();
    void sync_backing_unitcell_owned_atoms_();
    void clear_forces_(std::vector<LocalAtom>& atoms);
    static double wrap_fractional_(double value);

    std::int64_t nat_ = 0;
    double lat0_ = 0.0;
    double omega_ = 0.0;
    ModuleBase::Matrix3 latvec_;
    ModuleBase::Matrix3 gt_;
    std::vector<LocalAtom> owned_atoms_;
    std::vector<LocalAtom> ghost_atoms_;
    std::vector<std::string> type_labels_;
    std::vector<double> type_masses_;
    std::vector<std::int64_t> type_atom_counts_;
    MdStruFileMetadata stru_file_metadata_;
    double cutoff_ = 0.0;
    double skin_ = 0.0;
    std::unique_ptr<NeighborSearch> neighbor_search_;
    std::vector<ModuleBase::Vector3<double> > neighbor_reference_frac_;
    bool neighbor_layout_valid_ = false;
    UnitCell* backing_unitcell_ = nullptr;

#ifdef __MPI
    MPI_Comm comm_ = MPI_COMM_NULL;
    int rank_ = 0;
    int size_ = 1;
    DomainDecomposition decomp_;
#endif
};

#endif
