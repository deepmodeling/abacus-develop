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

BaseCell::Kind MDCell::kind() const
{
    return Kind::mdcell;
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
