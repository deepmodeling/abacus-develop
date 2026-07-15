#include "source_cell/md_cell.h"

#include "source_cell/unitcell.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

double MdCell::wrap_fractional_(double value)
{
    value -= std::floor(value);
    if (value >= 1.0 - 1.0e-12 || value < 1.0e-12)
    {
        return 0.0;
    }
    return value;
}

double MdCell::infer_cutoff_from_parameter_(const Parameter& param)
{
    double cutoff = 0.0;
    const std::vector<double>& lj_rcut = param.inp.mdp.lj_rcut;
    for (std::size_t i = 0; i < lj_rcut.size(); ++i)
    {
        cutoff = std::max(cutoff, lj_rcut[i] * ModuleBase::ANGSTROM_AU);
    }
    return cutoff;
}

void MdCell::clear_forces_(std::vector<LocalAtom>& atoms)
{
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        atoms[i].force.set(0.0, 0.0, 0.0);
    }
}

void MdCell::sync_backing_unitcell_geometry_()
{
    if (backing_unitcell_ == nullptr)
    {
        return;
    }

    backing_unitcell_->latvec = latvec_;
    backing_unitcell_->omega = omega_;
    backing_unitcell_->GT = gt_;
    backing_unitcell_->G = gt_.Transpose();
    backing_unitcell_->GGT = backing_unitcell_->G * backing_unitcell_->GT;
    backing_unitcell_->invGGT = backing_unitcell_->GGT.Inverse();
    backing_unitcell_->lat0_angstrom = lat0_ * ModuleBase::BOHR_TO_A;
    backing_unitcell_->tpiba = ModuleBase::TWO_PI / lat0_;
    backing_unitcell_->tpiba2 = backing_unitcell_->tpiba * backing_unitcell_->tpiba;
    backing_unitcell_->a1.set(latvec_.e11, latvec_.e12, latvec_.e13);
    backing_unitcell_->a2.set(latvec_.e21, latvec_.e22, latvec_.e23);
    backing_unitcell_->a3.set(latvec_.e31, latvec_.e32, latvec_.e33);
    backing_unitcell_->cell_parameter_updated = true;
}

void MdCell::sync_backing_unitcell_owned_atoms_()
{
    if (backing_unitcell_ == nullptr)
    {
        return;
    }

    for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
    {
        const LocalAtom& atom = owned_atoms_[i];
        backing_unitcell_->atoms[atom.type].tau[atom.type_index] = atom.cart;
        backing_unitcell_->atoms[atom.type].taud[atom.type_index] = atom.frac;
        backing_unitcell_->atoms[atom.type].vel[atom.type_index] = atom.vel;
        backing_unitcell_->atoms[atom.type].mbl[atom.type_index] = atom.mbl;
    }
}

void MdCell::initialize_from_ucell_serial_(const UnitCell& ucell, double cutoff, double skin)
{
    backing_unitcell_ = const_cast<UnitCell*>(&ucell);
    nat_ = ucell.nat;
    lat0_ = ucell.lat0;
    omega_ = ucell.omega;
    latvec_ = ucell.latvec;
    gt_ = ucell.GT;
    type_labels_ = ucell.atom_label;
    type_masses_ = ucell.atom_mass;
    init_vel_ = ucell.init_vel;
    cutoff_ = cutoff;
    skin_ = skin;
    owned_atoms_.clear();
    ghost_atoms_.clear();

    ModuleNeighList::GlobalAtomId global_id = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            owned_atoms_.push_back(LocalAtom(ucell.atoms[it].tau[ia],
                                             ucell.atoms[it].taud[ia],
                                             ucell.atoms[it].vel[ia],
                                             ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                             ucell.atoms[it].mbl[ia],
                                             ucell.atoms[it].mass / ModuleBase::AU_to_MASS,
                                             it,
                                             ia,
                                             global_id,
                                             0,
                                             false));
            ++global_id;
        }
    }
}

void MdCell::initialize_from_positions_serial_(const std::vector<ModuleBase::Vector3<double> >& positions,
                                               double cutoff,
                                               double skin)
{
    cutoff_ = cutoff;
    skin_ = skin;
    nat_ = static_cast<int>(positions.size());
    init_vel_ = false;
    owned_atoms_.clear();
    ghost_atoms_.clear();

    const ModuleBase::Matrix3 inv_latvec = latvec_.Inverse();
    for (std::size_t iat = 0; iat < positions.size(); ++iat)
    {
        ModuleBase::Vector3<double> frac = positions[iat] * inv_latvec;
        frac.x = wrap_fractional_(frac.x);
        frac.y = wrap_fractional_(frac.y);
        frac.z = wrap_fractional_(frac.z);
        owned_atoms_.push_back(LocalAtom(frac * latvec_,
                                         frac,
                                         ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                         ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                         ModuleBase::Vector3<int>(1, 1, 1),
                                         1.0,
                                         0,
                                         static_cast<int>(iat),
                                         static_cast<ModuleNeighList::GlobalAtomId>(iat),
                                         0,
                                         false));
    }
}

void MdCell::initialize_from_owned_atoms_serial_(double cutoff, double skin)
{
    cutoff_ = cutoff;
    skin_ = skin;
    if (nat_ == 0)
    {
        for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
        {
            nat_ = std::max(nat_, static_cast<int>(owned_atoms_[i].global_id + 1));
        }
    }
    clear_forces_(owned_atoms_);
    ghost_atoms_.clear();
}

#ifdef __MPI
void MdCell::initialize_from_ucell_(const UnitCell& ucell, MPI_Comm comm, double cutoff, double skin)
{
    backing_unitcell_ = const_cast<UnitCell*>(&ucell);
    nat_ = ucell.nat;
    lat0_ = ucell.lat0;
    omega_ = ucell.omega;
    latvec_ = ucell.latvec;
    gt_ = ucell.GT;
    type_labels_ = ucell.atom_label;
    type_masses_ = ucell.atom_mass;
    init_vel_ = ucell.init_vel;
    comm_ = comm;
    cutoff_ = cutoff;
    skin_ = skin;

    owned_atoms_.clear();
    ghost_atoms_.clear();
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    decomp_.init(comm_, latvec_, lat0_, cutoff_, skin_);
    decomp_.split_owned_atoms_from_ucell(ucell, owned_atoms_);
    clear_forces_(owned_atoms_);
    exchange_ghost_atoms();
}

void MdCell::initialize_from_positions_(const std::vector<ModuleBase::Vector3<double> >& positions,
                                        MPI_Comm comm,
                                        double cutoff,
                                        double skin)
{
    comm_ = comm;
    cutoff_ = cutoff;
    skin_ = skin;
    nat_ = static_cast<int>(positions.size());
    init_vel_ = false;

    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    decomp_.init(comm_, latvec_, lat0_, cutoff_, skin_);

    const ModuleBase::Matrix3 inv_latvec = latvec_.Inverse();
    owned_atoms_.clear();
    ghost_atoms_.clear();
    for (std::size_t iat = 0; iat < positions.size(); ++iat)
    {
        const ModuleBase::Vector3<double>& cart = positions[iat];
        ModuleBase::Vector3<double> frac = cart * inv_latvec;
        frac.x = wrap_fractional_(frac.x);
        frac.y = wrap_fractional_(frac.y);
        frac.z = wrap_fractional_(frac.z);

        const int owner = decomp_.owner_rank_from_frac(frac);
        if (owner == rank_)
        {
            owned_atoms_.push_back(LocalAtom(frac * latvec_,
                                             frac,
                                             ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                             ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                             ModuleBase::Vector3<int>(1, 1, 1),
                                             1.0,
                                             0,
                                             static_cast<int>(iat),
                                             static_cast<ModuleNeighList::GlobalAtomId>(iat),
                                             owner,
                                             false));
        }
    }

    exchange_ghost_atoms();
}

void MdCell::initialize_from_owned_atoms_(MPI_Comm comm, double cutoff, double skin)
{
    comm_ = comm;
    cutoff_ = cutoff;
    skin_ = skin;

    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    decomp_.init(comm_, latvec_, lat0_, cutoff_, skin_);
    clear_forces_(owned_atoms_);
    exchange_ghost_atoms();
}
#endif

MdCell::MdCell(const UnitCell& ucell, const Parameter& param)
{
    const double cutoff = infer_cutoff_from_parameter_(param);
    if (cutoff <= 0.0)
    {
        throw std::runtime_error("MdCell requires a positive LJ cutoff from Parameter.");
    }
#ifdef __MPI
    initialize_from_ucell_(ucell, MPI_COMM_WORLD, cutoff, 0.0);
#else
    initialize_from_ucell_serial_(ucell, cutoff, 0.0);
#endif
}

MdCell::MdCell(const UnitCell& ucell, double cutoff, double skin)
{
#ifdef __MPI
    initialize_from_ucell_(ucell, MPI_COMM_WORLD, cutoff, skin);
#else
    initialize_from_ucell_serial_(ucell, cutoff, skin);
#endif
}

MdCell::MdCell(const ModuleBase::Matrix3& latvec,
               const ModuleBase::Matrix3& gt,
               double lat0,
               double omega,
               const std::vector<ModuleBase::Vector3<double> >& positions,
               double cutoff,
               double skin)
{
    latvec_ = latvec;
    gt_ = gt;
    lat0_ = lat0;
    omega_ = omega;
    initialize_from_positions_serial_(positions, cutoff, skin);
}

MdCell::MdCell(const ModuleBase::Matrix3& latvec,
               const ModuleBase::Matrix3& gt,
               double lat0,
               double omega,
               const std::vector<LocalAtom>& owned_atoms,
               const std::vector<std::string>& type_labels,
               const std::vector<double>& type_masses,
               double cutoff,
               double skin)
{
    latvec_ = latvec;
    gt_ = gt;
    lat0_ = lat0;
    omega_ = omega;
    owned_atoms_ = owned_atoms;
    type_labels_ = type_labels;
    type_masses_ = type_masses;
    nat_ = 0;
    initialize_from_owned_atoms_serial_(cutoff, skin);
}

#ifdef __MPI
MdCell::MdCell(const UnitCell& ucell, const Parameter& param, MPI_Comm comm)
{
    const double cutoff = infer_cutoff_from_parameter_(param);
    if (cutoff <= 0.0)
    {
        throw std::runtime_error("MdCell requires a positive LJ cutoff from Parameter.");
    }
    initialize_from_ucell_(ucell, comm, cutoff, 0.0);
}

MdCell::MdCell(const UnitCell& ucell, MPI_Comm comm, double cutoff, double skin)
{
    initialize_from_ucell_(ucell, comm, cutoff, skin);
}

MdCell::MdCell(const ModuleBase::Matrix3& latvec,
               const ModuleBase::Matrix3& gt,
               double lat0,
               double omega,
               const std::vector<ModuleBase::Vector3<double> >& positions,
               MPI_Comm comm,
               double cutoff,
               double skin)
{
    latvec_ = latvec;
    gt_ = gt;
    lat0_ = lat0;
    omega_ = omega;
    initialize_from_positions_(positions, comm, cutoff, skin);
}

MdCell::MdCell(const ModuleBase::Matrix3& latvec,
               const ModuleBase::Matrix3& gt,
               double lat0,
               double omega,
               const std::vector<LocalAtom>& owned_atoms,
               const std::vector<std::string>& type_labels,
               const std::vector<double>& type_masses,
               MPI_Comm comm,
               double cutoff,
               double skin)
{
    latvec_ = latvec;
    gt_ = gt;
    lat0_ = lat0;
    omega_ = omega;
    owned_atoms_ = owned_atoms;
    type_labels_ = type_labels;
    type_masses_ = type_masses;
    nat_ = 0;
    for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
    {
        nat_ = std::max(nat_, static_cast<int>(owned_atoms_[i].global_id + 1));
    }
    MPI_Allreduce(MPI_IN_PLACE, &nat_, 1, MPI_INT, MPI_MAX, comm);
    init_vel_ = true;
    initialize_from_owned_atoms_(comm, cutoff, skin);
}

int MdCell::mpi_rank() const
{
    return rank_;
}

int MdCell::mpi_size() const
{
    return size_;
}

const DomainDecomposition& MdCell::decomposition() const
{
    return decomp_;
}
#endif

void MdCell::exchange_ghost_atoms()
{
#ifdef __MPI
    decomp_.exchange_ghost_atoms(owned_atoms_, ghost_atoms_);
    clear_forces_(ghost_atoms_);
    return;
#endif
    ghost_atoms_.clear();
}

void MdCell::migrate_owned_atoms()
{
#ifdef __MPI
    decomp_.migrate_owned_atoms(owned_atoms_);
    exchange_ghost_atoms();
    return;
#endif
    for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
    {
        LocalAtom& atom = owned_atoms_[i];
        atom.frac = atom.cart * gt_;
        atom.frac.x = wrap_fractional_(atom.frac.x);
        atom.frac.y = wrap_fractional_(atom.frac.y);
        atom.frac.z = wrap_fractional_(atom.frac.z);
        atom.cart = atom.frac * latvec_;
    }
    sync_backing_unitcell_owned_atoms_();
}

void MdCell::set_lattice_vectors(const ModuleBase::Matrix3& latvec)
{
    latvec_ = latvec;
    gt_ = latvec_.Inverse();
    omega_ = std::abs(latvec_.Det()) * lat0_ * lat0_ * lat0_;
#ifdef __MPI
    if (comm_ != MPI_COMM_NULL)
    {
        decomp_.init(comm_, latvec_, lat0_, cutoff_, skin_);
    }
#endif
    sync_backing_unitcell_geometry_();
}

void MdCell::refresh_cart_from_frac()
{
    for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
    {
        owned_atoms_[i].frac.x = wrap_fractional_(owned_atoms_[i].frac.x);
        owned_atoms_[i].frac.y = wrap_fractional_(owned_atoms_[i].frac.y);
        owned_atoms_[i].frac.z = wrap_fractional_(owned_atoms_[i].frac.z);
        owned_atoms_[i].cart = owned_atoms_[i].frac * latvec_;
    }
    sync_backing_unitcell_owned_atoms_();
    exchange_ghost_atoms();
}

const std::vector<LocalAtom>& MdCell::owned_atoms() const
{
    return owned_atoms_;
}

const std::vector<LocalAtom>& MdCell::ghost_atoms() const
{
    return ghost_atoms_;
}

const std::vector<std::string>& MdCell::type_labels() const
{
    return type_labels_;
}

const std::vector<double>& MdCell::type_masses() const
{
    return type_masses_;
}

std::vector<LocalAtom>& MdCell::mutable_owned_atoms()
{
    return owned_atoms_;
}

std::vector<LocalAtom>& MdCell::mutable_ghost_atoms()
{
    return ghost_atoms_;
}

int MdCell::nlocal() const
{
    return static_cast<int>(owned_atoms_.size());
}

int MdCell::nghost() const
{
    return static_cast<int>(ghost_atoms_.size());
}

bool MdCell::init_vel() const
{
    return init_vel_;
}

void MdCell::set_init_vel(bool init_vel)
{
    init_vel_ = init_vel;
}

double MdCell::cutoff() const
{
    return cutoff_;
}

double MdCell::skin() const
{
    return skin_;
}

bool MdCell::has_backing_unitcell() const
{
    return backing_unitcell_ != nullptr;
}

UnitCell& MdCell::backing_unitcell()
{
    assert(backing_unitcell_ != nullptr);
    return *backing_unitcell_;
}

const UnitCell& MdCell::backing_unitcell() const
{
    assert(backing_unitcell_ != nullptr);
    return *backing_unitcell_;
}

BaseCell::Kind MdCell::do_kind() const
{
    return Kind::md_cell;
}

int MdCell::do_nat() const
{
    return nat_;
}

double MdCell::do_lat0() const
{
    return lat0_;
}

double MdCell::do_omega() const
{
    return omega_;
}

const ModuleBase::Matrix3& MdCell::do_latvec() const
{
    return latvec_;
}

const ModuleBase::Matrix3& MdCell::do_GT() const
{
    return gt_;
}
