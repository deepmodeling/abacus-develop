#include "source_cell/mdcell.h"

#include "source_base/parallel_cell.h"
#include "source_cell/unitcell.h"
#include "source_cell/module_neighlist/neighbor_search.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

MDCell::MDCell() = default;
MDCell::~MDCell() = default;
MDCell::MDCell(MDCell&&) = default;
MDCell& MDCell::operator=(MDCell&&) = default;

double MDCell::wrap_fractional_(double value)
{
    value -= std::floor(value);
    if (value >= 1.0 - 1.0e-12 || value < 1.0e-12)
    {
        return 0.0;
    }
    return value;
}

void MDCell::clear_forces_(std::vector<LocalAtom>& atoms)
{
    for (std::size_t i = 0; i < atoms.size(); ++i)
    {
        atoms[i].force.set(0.0, 0.0, 0.0);
    }
}

void MDCell::sync_backing_unitcell_geometry_()
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
}

void MDCell::sync_backing_unitcell_owned_atoms_()
{
    if (backing_unitcell_ == nullptr)
    {
        return;
    }

    for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
    {
        const LocalAtom& atom = owned_atoms_[i];
        ModuleBase::Vector3<double> displacement = atom.frac - backing_unitcell_->atoms[atom.type].taud[atom.type_index];
        for (int k = 0; k < 3; ++k)
        {
            if (displacement[k] > 0.5)
            {
                displacement[k] -= 1.0;
            }
            else if (displacement[k] < -0.5)
            {
                displacement[k] += 1.0;
            }
        }
        backing_unitcell_->atoms[atom.type].tau[atom.type_index] = atom.cart;
        backing_unitcell_->atoms[atom.type].taud[atom.type_index] = atom.frac;
        backing_unitcell_->atoms[atom.type].dis[atom.type_index] = displacement;
        backing_unitcell_->atoms[atom.type].vel[atom.type_index] = atom.vel;
        backing_unitcell_->atoms[atom.type].mbl[atom.type_index] = atom.mbl;
    }
}

#ifdef __MPI
void MDCell::initialize_from_ucell_(UnitCell& ucell, MPI_Comm comm, double cutoff, double skin)
{
    backing_unitcell_ = &ucell;
    nat_ = ucell.nat;
    lat0_ = ucell.lat0;
    omega_ = ucell.omega;
    latvec_ = ucell.latvec;
    gt_ = ucell.GT;
    type_labels_.resize(static_cast<std::size_t>(ucell.ntype));
    type_masses_.resize(static_cast<std::size_t>(ucell.ntype));
    type_atom_counts_.resize(static_cast<std::size_t>(ucell.ntype));
    for (int it = 0; it < ucell.ntype; ++it)
    {
        type_labels_[static_cast<std::size_t>(it)] = ucell.atoms[it].label;
        type_masses_[static_cast<std::size_t>(it)] = ucell.atoms[it].mass;
        type_atom_counts_[static_cast<std::size_t>(it)] = ucell.atoms[it].na;
    }
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

void MDCell::initialize_from_owned_atoms_(MPI_Comm comm, double cutoff, double skin)
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
#else
void MDCell::initialize_from_ucell_(UnitCell& ucell, double cutoff, double skin)
{
    backing_unitcell_ = &ucell;
    nat_ = ucell.nat;
    lat0_ = ucell.lat0;
    omega_ = ucell.omega;
    latvec_ = ucell.latvec;
    gt_ = ucell.GT;
    type_labels_.resize(static_cast<std::size_t>(ucell.ntype));
    type_masses_.resize(static_cast<std::size_t>(ucell.ntype));
    type_atom_counts_.resize(static_cast<std::size_t>(ucell.ntype));
    for (int it = 0; it < ucell.ntype; ++it)
    {
        type_labels_[static_cast<std::size_t>(it)] = ucell.atoms[it].label;
        type_masses_[static_cast<std::size_t>(it)] = ucell.atoms[it].mass;
        type_atom_counts_[static_cast<std::size_t>(it)] = ucell.atoms[it].na;
    }
    cutoff_ = cutoff;
    skin_ = skin;
    owned_atoms_.clear();
    ghost_atoms_.clear();

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
                                             0));
        }
    }
    exchange_ghost_atoms();
}

void MDCell::initialize_from_owned_atoms_(double cutoff, double skin)
{
    cutoff_ = cutoff;
    skin_ = skin;
    clear_forces_(owned_atoms_);
    exchange_ghost_atoms();
}
#endif


void MDCell::initialize_from_unitcell(UnitCell& ucell,
                                      double cutoff,
                                      double skin,
                                      const ModuleBase::CommunicationDomain& comm_domain)
{
#ifdef __MPI
    initialize_from_ucell_(ucell, comm_domain.communicator(), cutoff, skin);
#else
    static_cast<void>(comm_domain);
    initialize_from_ucell_(ucell, cutoff, skin);
#endif
}

void MDCell::initialize_from_owned_atoms(const ModuleBase::Matrix3& latvec,
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
                                         const ModuleBase::CommunicationDomain& comm_domain)
{
    latvec_ = latvec;
    gt_ = gt;
    lat0_ = lat0;
    omega_ = omega;
    nat_ = nat;
    owned_atoms_ = owned_atoms;
    type_labels_ = type_labels;
    type_masses_ = type_masses;
    type_atom_counts_ = type_atom_counts;
#ifdef __MPI
    initialize_from_owned_atoms_(comm_domain.communicator(), cutoff, skin);
#else
    static_cast<void>(comm_domain);
    initialize_from_owned_atoms_(cutoff, skin);
#endif
}

#ifdef __MPI
int MDCell::mpi_rank() const
{
    return rank_;
}

int MDCell::mpi_size() const
{
    return size_;
}

#endif

void MDCell::exchange_ghost_atoms()
{
#ifdef __MPI
    decomp_.exchange_ghost_atoms(owned_atoms_, ghost_atoms_);
    clear_forces_(ghost_atoms_);
    return;
#endif

    ghost_atoms_.clear();

    if (cutoff_ <= 0.0)
    {
        return;
    }

    const ModuleBase::Vector3<double> a1(latvec_.e11, latvec_.e12, latvec_.e13);
    const ModuleBase::Vector3<double> a2(latvec_.e21, latvec_.e22, latvec_.e23);
    const ModuleBase::Vector3<double> a3(latvec_.e31, latvec_.e32, latvec_.e33);
    const ModuleBase::Vector3<double> a2xa3(a2.y * a3.z - a2.z * a3.y,
                                             a2.z * a3.x - a2.x * a3.z,
                                             a2.x * a3.y - a2.y * a3.x);
    const ModuleBase::Vector3<double> a3xa1(a3.y * a1.z - a3.z * a1.y,
                                             a3.z * a1.x - a3.x * a1.z,
                                             a3.x * a1.y - a3.y * a1.x);
    const ModuleBase::Vector3<double> a1xa2(a1.y * a2.z - a1.z * a2.y,
                                             a1.z * a2.x - a1.x * a2.z,
                                             a1.x * a2.y - a1.y * a2.x);
    const double volume = std::abs(a1.x * a2xa3.x + a1.y * a2xa3.y + a1.z * a2xa3.z);
    if (volume <= 0.0)
    {
        throw std::runtime_error("MDCell requires a nonzero cell volume for periodic ghosts.");
    }

    const double search_radius = (cutoff_ + skin_) / lat0_;
    const int layers[3] = {
        static_cast<int>(std::ceil(a2xa3.norm() * search_radius / volume)),
        static_cast<int>(std::ceil(a3xa1.norm() * search_radius / volume)),
        static_cast<int>(std::ceil(a1xa2.norm() * search_radius / volume))
    };
    for (int ix = -layers[0]; ix <= layers[0]; ++ix)
    {
        for (int iy = -layers[1]; iy <= layers[1]; ++iy)
        {
            for (int iz = -layers[2]; iz <= layers[2]; ++iz)
            {
                if (ix == 0 && iy == 0 && iz == 0)
                {
                    continue;
                }
                for (std::size_t iat = 0; iat < owned_atoms_.size(); ++iat)
                {
                    LocalAtom image = owned_atoms_[iat];
                    const ModuleBase::Vector3<double> shifted_frac(image.frac.x + ix,
                                                                     image.frac.y + iy,
                                                                     image.frac.z + iz);
                    image.cart = shifted_frac * latvec_;
                    image.force.set(0.0, 0.0, 0.0);
                    ghost_atoms_.push_back(image);
                }
            }
        }
    }
}

void MDCell::accumulate_ghost_forces()
{
#ifdef __MPI
    decomp_.accumulate_ghost_forces(owned_atoms_, ghost_atoms_);
#else
    for (std::size_t ighost = 0; ighost < ghost_atoms_.size(); ++ighost)
    {
        const LocalAtom& ghost = ghost_atoms_[ighost];
        for (std::size_t iowned = 0; iowned < owned_atoms_.size(); ++iowned)
        {
            LocalAtom& owned = owned_atoms_[iowned];
            if (owned.type == ghost.type && owned.type_index == ghost.type_index)
            {
                owned.force += ghost.force;
                break;
            }
        }
    }
#endif
}

void MDCell::migrate_owned_atoms()
{
#ifdef __MPI
    decomp_.migrate_owned_atoms(owned_atoms_);
    exchange_ghost_atoms();
    neighbor_layout_valid_ = false;
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
    exchange_ghost_atoms();
    neighbor_layout_valid_ = false;
}

void MDCell::prepare_neighbors()
{
    bool rebuild = !neighbor_layout_valid_ || neighbor_reference_frac_.size() != owned_atoms_.size();
    double local_max_displacement = 0.0;
    if (!rebuild)
    {
        for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
        {
            ModuleBase::Vector3<double> delta = owned_atoms_[i].frac - neighbor_reference_frac_[i];
            delta.x -= std::nearbyint(delta.x);
            delta.y -= std::nearbyint(delta.y);
            delta.z -= std::nearbyint(delta.z);
            local_max_displacement = std::max(local_max_displacement, (delta * latvec_).norm() * lat0_);
        }
#ifdef __MPI
        if (comm_ != MPI_COMM_NULL)
        {
            MPI_Allreduce(MPI_IN_PLACE, &local_max_displacement, 1, MPI_DOUBLE, MPI_MAX, comm_);
        }
#endif
        rebuild = local_max_displacement >= skin_ * 0.5;
    }

    if (rebuild)
    {
        migrate_owned_atoms();
        neighbor_search_.reset(new NeighborSearch);
        neighbor_search_->init(*this, cutoff_ + skin_);
        neighbor_search_->build_neighbors();
        neighbor_search_->refresh_mdcell(*this, cutoff_);
        neighbor_reference_frac_.resize(owned_atoms_.size());
        for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
        {
            neighbor_reference_frac_[i] = owned_atoms_[i].frac;
        }
        neighbor_layout_valid_ = true;
        return;
    }

#ifdef __MPI
    decomp_.update_ghost_atom_positions(owned_atoms_, ghost_atoms_);
#else
    exchange_ghost_atoms();
#endif
    neighbor_search_->refresh_mdcell(*this, cutoff_);
}

const NeighborSearch& MDCell::neighbor_search() const
{
    if (!neighbor_search_)
    {
        throw std::runtime_error("MDCell neighbor list has not been prepared.");
    }
    return *neighbor_search_;
}

bool MDCell::has_neighbor_search() const
{
    return neighbor_layout_valid_ && neighbor_search_ != NULL;
}

void MDCell::set_lattice_vectors(const ModuleBase::Matrix3& latvec)
{
    latvec_ = latvec;
    gt_ = latvec_.Inverse();
    omega_ = std::abs(latvec_.Det()) * lat0_ * lat0_ * lat0_;
    neighbor_layout_valid_ = false;
#ifdef __MPI
    if (comm_ != MPI_COMM_NULL)
    {
        decomp_.init(comm_, latvec_, lat0_, cutoff_, skin_);
    }
#endif
    sync_backing_unitcell_geometry_();
    if (backing_unitcell_ != nullptr)
    {
        backing_unitcell_->cell_parameter_updated = true;
    }
}

void MDCell::refresh_cart_from_frac()
{
    for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
    {
        owned_atoms_[i].frac.x = wrap_fractional_(owned_atoms_[i].frac.x);
        owned_atoms_[i].frac.y = wrap_fractional_(owned_atoms_[i].frac.y);
        owned_atoms_[i].frac.z = wrap_fractional_(owned_atoms_[i].frac.z);
        owned_atoms_[i].cart = owned_atoms_[i].frac * latvec_;
    }
    neighbor_layout_valid_ = false;
}

const std::vector<LocalAtom>& MDCell::ghost_atoms() const
{
    return ghost_atoms_;
}

std::vector<LocalAtom>& MDCell::mutable_owned_atoms()
{
    return owned_atoms_;
}

std::vector<LocalAtom>& MDCell::mutable_ghost_atoms()
{
    return ghost_atoms_;
}

double MDCell::cutoff() const
{
    return cutoff_;
}

bool MDCell::has_backing_unitcell() const
{
    return backing_unitcell_ != nullptr;
}

UnitCell& MDCell::backing_unitcell()
{
    assert(backing_unitcell_ != nullptr);
    return *backing_unitcell_;
}

const UnitCell& MDCell::backing_unitcell() const
{
    assert(backing_unitcell_ != nullptr);
    return *backing_unitcell_;
}

void MDCell::sync_backing_unitcell()
{
    if (backing_unitcell_ == nullptr)
    {
        return;
    }

    sync_backing_unitcell_geometry_();

#ifdef __MPI
    if (size_ > 1)
    {
        std::vector<int> type_offset(backing_unitcell_->ntype + 1, 0);
        for (int it = 0; it < backing_unitcell_->ntype; ++it)
        {
            type_offset[it + 1] = type_offset[it] + backing_unitcell_->atoms[it].na;
        }

        std::vector<double> cart(3 * nat_, 0.0);
        std::vector<double> frac(3 * nat_, 0.0);
        std::vector<double> vel(3 * nat_, 0.0);
        std::vector<int> mbl(3 * nat_, 0);
        std::vector<int> owner(nat_, 0);

        for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
        {
            const LocalAtom& atom = owned_atoms_[i];
            const int iat = type_offset[atom.type] + atom.type_index;
            cart[3 * iat] = atom.cart.x;
            cart[3 * iat + 1] = atom.cart.y;
            cart[3 * iat + 2] = atom.cart.z;
            frac[3 * iat] = atom.frac.x;
            frac[3 * iat + 1] = atom.frac.y;
            frac[3 * iat + 2] = atom.frac.z;
            vel[3 * iat] = atom.vel.x;
            vel[3 * iat + 1] = atom.vel.y;
            vel[3 * iat + 2] = atom.vel.z;
            mbl[3 * iat] = atom.mbl.x;
            mbl[3 * iat + 1] = atom.mbl.y;
            mbl[3 * iat + 2] = atom.mbl.z;
            owner[iat] = 1;
        }

        MPI_Allreduce(MPI_IN_PLACE, cart.data(), 3 * nat_, MPI_DOUBLE, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, frac.data(), 3 * nat_, MPI_DOUBLE, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, vel.data(), 3 * nat_, MPI_DOUBLE, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, mbl.data(), 3 * nat_, MPI_INT, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, owner.data(), nat_, MPI_INT, MPI_SUM, comm_);

        for (int it = 0; it < backing_unitcell_->ntype; ++it)
        {
            for (int ia = 0; ia < backing_unitcell_->atoms[it].na; ++ia)
            {
                const int iat = type_offset[it] + ia;
                if (owner[iat] != 1)
                {
                    throw std::runtime_error("MDCell backing UnitCell atom ownership is invalid.");
                }
                backing_unitcell_->atoms[it].tau[ia].set(cart[3 * iat], cart[3 * iat + 1], cart[3 * iat + 2]);
                ModuleBase::Vector3<double> displacement(frac[3 * iat] - backing_unitcell_->atoms[it].taud[ia].x,
                                                         frac[3 * iat + 1] - backing_unitcell_->atoms[it].taud[ia].y,
                                                         frac[3 * iat + 2] - backing_unitcell_->atoms[it].taud[ia].z);
                for (int k = 0; k < 3; ++k)
                {
                    if (displacement[k] > 0.5)
                    {
                        displacement[k] -= 1.0;
                    }
                    else if (displacement[k] < -0.5)
                    {
                        displacement[k] += 1.0;
                    }
                }
                backing_unitcell_->atoms[it].taud[ia].set(frac[3 * iat], frac[3 * iat + 1], frac[3 * iat + 2]);
                backing_unitcell_->atoms[it].dis[ia] = displacement;
                backing_unitcell_->atoms[it].vel[ia].set(vel[3 * iat], vel[3 * iat + 1], vel[3 * iat + 2]);
                backing_unitcell_->atoms[it].mbl[ia].set(mbl[3 * iat], mbl[3 * iat + 1], mbl[3 * iat + 2]);
            }
        }
        return;
    }
#endif

    sync_backing_unitcell_owned_atoms_();
}

BaseCell::Kind MDCell::get_kind() const
{
    return Kind::mdcell;
}

std::int64_t MDCell::get_nat() const
{
    return nat_;
}

double MDCell::get_lat0() const
{
    return lat0_;
}

double MDCell::get_omega() const
{
    return omega_;
}

const ModuleBase::Matrix3& MDCell::get_latvec() const
{
    return latvec_;
}

const ModuleBase::Matrix3& MDCell::get_GT() const
{
    return gt_;
}
