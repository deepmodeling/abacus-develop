#ifndef MD_CELL_H
#define MD_CELL_H

#include "source_cell/base_cell.h"
#include "source_cell/module_neighlist/local_atom.h"

#ifdef __MPI
#include "source_cell/module_neighlist/domain_decomposition.h"
#endif

#include <string>
#include <vector>

class UnitCell;
class MDCell : public BaseCell
{
public:
    MDCell(UnitCell& ucell, double cutoff, double skin);
    MDCell(const MDCell&) = delete;
    MDCell& operator=(const MDCell&) = delete;
    MDCell(MDCell&&) = default;
    MDCell& operator=(MDCell&&) = default;

    MDCell(const ModuleBase::Matrix3& latvec,
           const ModuleBase::Matrix3& gt,
           double lat0,
           double omega,
           int nat,
           const std::vector<LocalAtom>& owned_atoms,
           const std::vector<std::string>& type_labels,
           const std::vector<double>& type_masses,
           double cutoff,
           double skin);

#ifdef __MPI
    int mpi_rank() const;
    int mpi_size() const;
    MPI_Comm communicator() const;

    const DomainDecomposition& decomposition() const;
#endif

    void exchange_ghost_atoms();
    void accumulate_ghost_forces();
    void migrate_owned_atoms();
    void set_lattice_vectors(const ModuleBase::Matrix3& latvec);
    void refresh_cart_from_frac();

    const std::vector<LocalAtom>& owned_atoms() const;
    const std::vector<LocalAtom>& ghost_atoms() const;
    const std::vector<std::string>& type_labels() const;
    const std::vector<double>& type_masses() const;
    std::vector<LocalAtom>& mutable_owned_atoms();
    std::vector<LocalAtom>& mutable_ghost_atoms();

    int nlocal() const;
    int nghost() const;
    bool init_vel() const;
    void set_init_vel(bool init_vel);
    double cutoff() const;
    double skin() const;
    bool has_backing_unitcell() const;
    UnitCell& backing_unitcell();
    const UnitCell& backing_unitcell() const;
    void sync_backing_unitcell();

private:
    Kind get_kind() const override;
    int get_nat() const override;
    double get_lat0() const override;
    double get_omega() const override;
    const ModuleBase::Matrix3& get_latvec() const override;
    const ModuleBase::Matrix3& get_GT() const override;

#ifdef __MPI
    void initialize_from_ucell_(UnitCell& ucell, MPI_Comm comm, double cutoff, double skin);
#endif
    void initialize_from_ucell_serial_(UnitCell& ucell, double cutoff, double skin);
    void sync_backing_unitcell_geometry_();
    void sync_backing_unitcell_owned_atoms_();
    void clear_forces_(std::vector<LocalAtom>& atoms);
    static double wrap_fractional_(double value);
#ifdef __MPI
    void initialize_from_owned_atoms_(MPI_Comm comm, double cutoff, double skin);
#else
    void initialize_from_owned_atoms_(double cutoff, double skin);
#endif

    int nat_ = 0;
    double lat0_ = 0.0;
    double omega_ = 0.0;
    ModuleBase::Matrix3 latvec_;
    ModuleBase::Matrix3 gt_;
    std::vector<LocalAtom> owned_atoms_;
    std::vector<LocalAtom> ghost_atoms_;
    std::vector<std::string> type_labels_;
    std::vector<double> type_masses_;
    bool init_vel_ = false;
    double cutoff_ = 0.0;
    double skin_ = 0.0;
    UnitCell* backing_unitcell_ = nullptr;

#ifdef __MPI
    MPI_Comm comm_ = MPI_COMM_NULL;
    int rank_ = 0;
    int size_ = 1;
    DomainDecomposition decomp_;
#endif
};

#endif
