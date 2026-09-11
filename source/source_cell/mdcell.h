#ifndef MDCELL_H
#define MDCELL_H

#include "source_cell/basecell.h"
#include "source_cell/strumeta.h"
#include "source_cell/module_neighlist/local_atom.h"
#include "source_base/matrix3.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <string>
#include <cstdint>
#include <memory>
#include <vector>

class UnitCell;
class NeighborSearch;
class DomainDecomposition;
namespace ModuleBase
{
class CommunicationDomain;
}

class MDCell : public BaseCell
{
public:
    MDCell();
    ~MDCell();
    MDCell(MDCell&&);
    MDCell(const MDCell&) = delete;
    MDCell& operator=(const MDCell&) = delete;
    MDCell& operator=(MDCell&&);

    Kind kind() const override;

    void initialize_from_owned_atoms(const ModuleBase::Matrix3& latvec,
                                     const ModuleBase::Matrix3& gt,
                                     double lat0,
                                     double omega,
                                     std::int64_t nat,
                                     const std::vector<LocalAtom>& owned_atoms,
                                     const std::vector<std::string>& type_labels,
                                     const std::vector<double>& type_masses,
                                     const std::vector<std::int64_t>& type_atom_counts,
                                     double skin,
                                     const ModuleBase::CommunicationDomain& comm_domain);

    void set_neighbor_cutoff(double cutoff);

    bool has_neighbor_search() const;
    const NeighborSearch& neighbor_search() const;
    void set_lattice_vectors(const ModuleBase::Matrix3& latvec);
    void refresh_cart_from_frac();

    void sync_backing_unitcell();

#ifdef __MPI
    int mpi_rank() const;
    int mpi_size() const;
    MPI_Comm communicator() const { return comm_; }

#endif

public:
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
    StruMeta stru_meta_;

    UnitCell* backing_unitcell_ = nullptr;

private:
    friend class DomainDecomposition;
    void sync_backing_unitcell_geometry_();
    void clear_forces_(std::vector<LocalAtom>& atoms);
    static double wrap_fractional_(double value);

    double cutoff_ = 0.0;
    double skin_ = 0.0;
    std::unique_ptr<NeighborSearch> neighbor_search_;
    std::vector<ModuleBase::Vector3<double> > neighbor_reference_frac_;
    bool neighbor_layout_valid_ = false;

#ifdef __MPI
    MPI_Comm comm_ = MPI_COMM_NULL;
    int rank_ = 0;
    int size_ = 1;
#endif
};

#endif
