#include "source_md/md_state_view.h"

#include "source_cell/md_cell.h"
#include "source_cell/unitcell.h"

#include <cassert>
#include <stdexcept>

namespace
{
MdAtomRef& checked_atom_ref(std::vector<MdAtomRef>& atoms, const int i)
{
    if (i < 0 || i >= static_cast<int>(atoms.size()))
    {
        throw std::out_of_range("MdStateView atom index out of range.");
    }
    return atoms[static_cast<std::size_t>(i)];
}

const MdAtomRef& checked_atom_ref(const std::vector<MdAtomRef>& atoms, const int i)
{
    if (i < 0 || i >= static_cast<int>(atoms.size()))
    {
        throw std::out_of_range("MdStateView atom index out of range.");
    }
    return atoms[static_cast<std::size_t>(i)];
}
}

MdStateView MdStateView::from_unitcell(UnitCell& ucell,
                                       ModuleBase::Vector3<double>* force,
                                       ModuleBase::Vector3<double>* vel,
                                       ModuleBase::Vector3<int>* mbl,
                                       double* mass)
{
    MdStateView view;
    view.distributed_ = false;
    view.atoms_.reserve(static_cast<std::size_t>(ucell.nat));

    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        Atom& atom = ucell.atoms[it];
        for (int ia = 0; ia < atom.na; ++ia, ++iat)
        {
            MdAtomRef ref;
            ref.cart = &atom.tau[ia];
            ref.frac = &atom.taud[ia];
            ref.vel = vel != 0 ? &vel[iat] : &atom.vel[ia];
            ref.mbl = mbl != 0 ? &mbl[iat] : &atom.mbl[ia];
            ref.force = force != 0 ? &force[iat]
                                   : (static_cast<int>(atom.force.size()) == atom.na ? &atom.force[ia] : 0);
            ref.mass = mass != 0 ? mass[iat] : atom.mass;
            ref.type = it;
            ref.type_index = ia;
            ref.global_id = iat;
            view.atoms_.push_back(ref);
        }
    }

    assert(static_cast<int>(view.atoms_.size()) == ucell.nat);
    return view;
}

MdStateView MdStateView::from_mdcell(MdCell& mdcell)
{
    MdStateView view;
    view.distributed_ = true;
    const std::vector<LocalAtom>& owned_atoms = mdcell.owned_atoms();
    view.atoms_.reserve(owned_atoms.size());

    for (std::size_t i = 0; i < owned_atoms.size(); ++i)
    {
        LocalAtom& atom = mdcell.mutable_owned_atoms()[i];
        MdAtomRef ref;
        ref.cart = &atom.cart;
        ref.frac = &atom.frac;
        ref.vel = &atom.vel;
        ref.mbl = &atom.mbl;
        ref.force = &atom.force;
        ref.mass = atom.mass;
        ref.type = atom.type;
        ref.type_index = atom.type_index;
        ref.global_id = atom.global_id;
        view.atoms_.push_back(ref);
    }

    return view;
}

int MdStateView::size() const
{
    return static_cast<int>(atoms_.size());
}

bool MdStateView::uses_distributed_storage() const
{
    return distributed_;
}

ModuleBase::Vector3<double>& MdStateView::cart(const int i) const
{
    MdAtomRef& atom = checked_atom_ref(const_cast<std::vector<MdAtomRef>&>(atoms_), i);
    assert(atom.cart != 0);
    return *atom.cart;
}

ModuleBase::Vector3<double>& MdStateView::frac(const int i) const
{
    MdAtomRef& atom = checked_atom_ref(const_cast<std::vector<MdAtomRef>&>(atoms_), i);
    assert(atom.frac != 0);
    return *atom.frac;
}

ModuleBase::Vector3<double>& MdStateView::vel(const int i) const
{
    MdAtomRef& atom = checked_atom_ref(const_cast<std::vector<MdAtomRef>&>(atoms_), i);
    assert(atom.vel != 0);
    return *atom.vel;
}

ModuleBase::Vector3<int>& MdStateView::mbl(const int i) const
{
    MdAtomRef& atom = checked_atom_ref(const_cast<std::vector<MdAtomRef>&>(atoms_), i);
    assert(atom.mbl != 0);
    return *atom.mbl;
}

ModuleBase::Vector3<double>& MdStateView::force(const int i) const
{
    MdAtomRef& atom = checked_atom_ref(const_cast<std::vector<MdAtomRef>&>(atoms_), i);
    if (atom.force == 0)
    {
        throw std::runtime_error("MdStateView force() requested but no force storage is attached.");
    }
    return *atom.force;
}

double MdStateView::mass(const int i) const
{
    return checked_atom_ref(atoms_, i).mass;
}

int MdStateView::type(const int i) const
{
    return checked_atom_ref(atoms_, i).type;
}

int MdStateView::type_index(const int i) const
{
    return checked_atom_ref(atoms_, i).type_index;
}

ModuleNeighList::GlobalAtomId MdStateView::global_id(const int i) const
{
    return checked_atom_ref(atoms_, i).global_id;
}
