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

const GhostExchangePlan& MDCell::get_ghost_exchange_plan() const
{
    if (!ghost_exchange_plan_valid_)
    {
        throw std::runtime_error("MDCell ghost exchange plan has not been prepared.");
    }
    return ghost_exchange_plan_;
}

GhostExchangePlan& MDCell::get_ghost_exchange_plan()
{
    if (!ghost_exchange_plan_valid_)
    {
        ModuleBase::WARNING_QUIT("MDCell", "Ghost exchange plan is not available.");
    }
    return ghost_exchange_plan_;
}

void MDCell::clear_ghost_exchange_plan_()
{
    ghost_exchange_plan_ = GhostExchangePlan();
    ghost_exchange_plan_valid_ = false;
}

void MDCell::set_ghost_exchange_plan_(const GhostExchangePlan& plan)
{
    ghost_exchange_plan_ = plan;
    ghost_exchange_plan_valid_ = true;
}

void MDCell::set_backing_unitcell(UnitCell& ucell)
{
    backing_unitcell_ = &ucell;
}

UnitCell& MDCell::backing_unitcell()
{
    return *backing_unitcell_;
}

const UnitCell& MDCell::backing_unitcell() const
{
    return *backing_unitcell_;
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
    backing_unitcell_ = nullptr;
    cutoff_ = 0.0;
    skin_ = skin;
    neighbor_search_.reset();
    neighbor_layout_valid_ = false;
    ghost_atoms_.clear();
    clear_ghost_exchange_plan_();
#ifdef __MPI
    comm_ = comm_domain.communicator();
    rank_ = comm_domain.rank();
    size_ = comm_domain.size();
#else
    static_cast<void>(comm_domain);
#endif
    clear_forces_(owned_atoms_);
}

void MDCell::set_neighbor_cutoff(double cutoff)
{
    if (cutoff <= 0.0)
    {
        throw std::runtime_error("MDCell neighbor cutoff must be positive.");
    }
    cutoff_ = cutoff;
    neighbor_search_.reset();
    neighbor_reference_frac_.clear();
    neighbor_layout_valid_ = false;
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

void MDCell::sync_backing_unitcell(const bool positions_updated)
{
    if (backing_unitcell_ == nullptr)
    {
        return;
    }

    sync_backing_unitcell_geometry_();
    for (int it = 0; it < backing_unitcell_->ntype; ++it)
    {
        backing_unitcell_->atoms[it].boundary_shift.resize(backing_unitcell_->atoms[it].na);
    }

#ifdef __MPI
    if (size_ > 1)
    {
        std::vector<int> type_offset(backing_unitcell_->ntype + 1, 0);
        for (int it = 0; it < backing_unitcell_->ntype; ++it)
        {
            type_offset[it + 1] = type_offset[it] + backing_unitcell_->atoms[it].na;
        }

        const int coordinate_count = 3 * nat_;
        std::vector<double> cart(coordinate_count, 0.0);
        std::vector<double> frac(coordinate_count, 0.0);
        std::vector<double> vel(coordinate_count, 0.0);
        std::vector<double> displacement(coordinate_count, 0.0);
        std::vector<int> boundary_shift(coordinate_count, 0);
        std::vector<int> mbl(coordinate_count, 0);
        std::vector<int> owner(nat_, 0);

        for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
        {
            const LocalAtom& atom = owned_atoms_[i];
            const int iat = type_offset[atom.type] + atom.type_index;
            for (int axis = 0; axis < 3; ++axis)
            {
                const int component = 3 * iat + axis;
                cart[component] = atom.cart[axis];
                frac[component] = atom.frac[axis];
                vel[component] = atom.vel[axis];
                displacement[component] = atom.displacement[axis];
                boundary_shift[component] = atom.boundary_shift[axis];
                mbl[component] = atom.mbl[axis];
            }
            owner[iat] = 1;
        }

        ModuleBase::CommunicationDomain domain;
        domain.initialize(comm_);
        domain.sum(cart.data(), coordinate_count);
        domain.sum(frac.data(), coordinate_count);
        domain.sum(vel.data(), coordinate_count);
        domain.sum(displacement.data(), coordinate_count);
        domain.sum(boundary_shift.data(), coordinate_count);
        domain.sum(mbl.data(), coordinate_count);
        const int atom_count = nat_;
        domain.sum(owner.data(), atom_count);

        for (int it = 0; it < backing_unitcell_->ntype; ++it)
        {
            for (int ia = 0; ia < backing_unitcell_->atoms[it].na; ++ia)
            {
                const int iat = type_offset[it] + ia;
                if (owner[iat] != 1)
                {
                    throw std::runtime_error("MDCell backing UnitCell atom ownership is invalid.");
                }
                for (int axis = 0; axis < 3; ++axis)
                {
                    const int component = 3 * iat + axis;
                    backing_unitcell_->atoms[it].tau[ia][axis] = cart[component];
                    backing_unitcell_->atoms[it].taud[ia][axis] = frac[component];
                    backing_unitcell_->atoms[it].vel[ia][axis] = vel[component];
                    backing_unitcell_->atoms[it].dis[ia][axis] = displacement[component];
                    backing_unitcell_->atoms[it].boundary_shift[ia][axis] = boundary_shift[component];
                    backing_unitcell_->atoms[it].mbl[ia][axis] = mbl[component];
                }
            }
        }
        if (positions_updated)
        {
            backing_unitcell_->ionic_position_updated = true;
        }
        return;
    }
#endif

    for (std::size_t i = 0; i < owned_atoms_.size(); ++i)
    {
        const LocalAtom& atom = owned_atoms_[i];
        backing_unitcell_->atoms[atom.type].tau[atom.type_index] = atom.cart;
        backing_unitcell_->atoms[atom.type].taud[atom.type_index] = atom.frac;
        backing_unitcell_->atoms[atom.type].dis[atom.type_index] = atom.displacement;
        backing_unitcell_->atoms[atom.type].boundary_shift[atom.type_index] = atom.boundary_shift;
        backing_unitcell_->atoms[atom.type].vel[atom.type_index] = atom.vel;
        backing_unitcell_->atoms[atom.type].mbl[atom.type_index] = atom.mbl;
    }
    if (positions_updated)
    {
        backing_unitcell_->ionic_position_updated = true;
    }
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
