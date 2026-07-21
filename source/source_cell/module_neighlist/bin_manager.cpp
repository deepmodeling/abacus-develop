#include <limits>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include "bin_manager.h"

// ========== Bin class implementation ==========

const std::vector<ModuleNeighList::LocalAtomIndex>& Bin::get_atom_indices() const {
    return atom_indices_;
}

void Bin::set_id(int ix, int iy, int iz) {
    id_x_ = ix;
    id_y_ = iy;
    id_z_ = iz;
}

void Bin::clear_atoms() {
    atom_indices_.clear();
}

void Bin::add_atom_index(ModuleNeighList::LocalAtomIndex atom_index) {
    atom_indices_.push_back(atom_index);
}

// ========== BinManager getter methods ==========

int BinManager::get_nbinx() const {
    return nbinx_;
}

int BinManager::get_nbiny() const {
    return nbiny_;
}

int BinManager::get_nbinz() const {
    return nbinz_;
}

int BinManager::get_total_bins() const {
    return ModuleNeighList::checked_int_size(bins_.size(), "BinManager total bin count");
}

int BinManager::get_bin_atom_count(int bin_index) const {
    if (bin_index < 0 || static_cast<std::size_t>(bin_index) >= bins_.size()) {
        return 0;
    }
    return ModuleNeighList::checked_int_size(bins_[bin_index].get_atom_indices().size(),
                                             "Bin atom count");
}

// ========== BinManager main methods ==========

void BinManager::init_bins(
    double sr,
    const std::vector<NeighborAtom>& all_atoms
)
{
    sradius_ = sr;
    if (!std::isfinite(sradius_) || sradius_ <= 0.0)
    {
        throw std::invalid_argument("BinManager search radius must be finite and positive.");
    }
    if(all_atoms.empty())
    {
        x_min_ = y_min_ = z_min_ = 0;
        x_max_ = y_max_ = z_max_ = 0;
        nbinx_ = nbiny_ = nbinz_ = 1;
        bins_.clear();
        bins_.resize(1);
        return;
    }

    x_min_ = y_min_ = z_min_ = std::numeric_limits<double>::max();
    x_max_ = y_max_ = z_max_ = std::numeric_limits<double>::lowest();

    auto update_bounds = [&](const std::vector<NeighborAtom>& atoms)
    {
        for (const auto& atom : atoms)
        {
            x_min_ = std::min(x_min_, atom.position_x);
            x_max_ = std::max(x_max_, atom.position_x);

            y_min_ = std::min(y_min_, atom.position_y);
            y_max_ = std::max(y_max_, atom.position_y);

            z_min_ = std::min(z_min_, atom.position_z);
            z_max_ = std::max(z_max_, atom.position_z);
        }
    };

    update_bounds(all_atoms);

    bin_sizex_ = bin_sizey_ = bin_sizez_ = sradius_;

    const auto checked_bin_dimension = [](const double span, const double bin_size, const char* context) {
        const double count = std::ceil(span / bin_size);
        if (!std::isfinite(count) || count > static_cast<double>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error(std::string(context) + " exceeds int range.");
        }
        return static_cast<int>(count);
    };

    nbinx_ = checked_bin_dimension(x_max_ - x_min_, bin_sizex_, "BinManager X bin count");
    nbiny_ = checked_bin_dimension(y_max_ - y_min_, bin_sizey_, "BinManager Y bin count");
    nbinz_ = checked_bin_dimension(z_max_ - z_min_, bin_sizez_, "BinManager Z bin count");

    nbinx_ = std::max(1, nbinx_);
    nbiny_ = std::max(1, nbiny_);
    nbinz_ = std::max(1, nbinz_);

    const std::size_t nbins_xy = ModuleNeighList::checked_size_product(static_cast<std::size_t>(nbinx_),
                                                                        static_cast<std::size_t>(nbiny_),
                                                                        "BinManager bin count");
    const std::size_t nbins_size = ModuleNeighList::checked_size_product(nbins_xy,
                                                                         static_cast<std::size_t>(nbinz_),
                                                                         "BinManager bin count");
    const int nbins = ModuleNeighList::checked_int_size(nbins_size, "BinManager bin count");

    bins_.clear();
    bins_.resize(nbins);

    for (int ix = 0; ix < nbinx_; ++ix)
    {
        for (int iy = 0; iy < nbiny_; ++iy)
        {
            for (int iz = 0; iz < nbinz_; ++iz)
            {
                int idx = bin_index(ix, iy, iz);

                bins_[idx].set_id(ix, iy, iz);
                bins_[idx].clear_atoms();
            }
        }
    }
}

void BinManager::do_binning(
    const std::vector<NeighborAtom>& atoms
)
{
    if (atoms.size() > static_cast<std::size_t>(std::numeric_limits<ModuleNeighList::LocalAtomIndex>::max()))
    {
        throw std::overflow_error("BinManager binned atom count exceeds local atom index range.");
    }

    for (std::size_t iatom = 0; iatom < atoms.size(); ++iatom)
    {
        const NeighborAtom& atom = atoms[iatom];
        int ix = std::min(
            std::max(int((atom.position_x - x_min_) / bin_sizex_), 0),
            nbinx_ - 1
        );

        int iy = std::min(
            std::max(int((atom.position_y - y_min_) / bin_sizey_), 0),
            nbiny_ - 1
        );

        int iz = std::min(
            std::max(int((atom.position_z - z_min_) / bin_sizez_), 0),
            nbinz_ - 1
        );

        int idx = bin_index(ix, iy, iz);

        const ModuleNeighList::LocalAtomIndex atom_index = static_cast<ModuleNeighList::LocalAtomIndex>(iatom);
        bins_[idx].add_atom_index(atom_index);
    }
}

int BinManager::bin_index(int ix, int iy, int iz) const {
    return ix * nbiny_ * nbinz_ + iy * nbinz_ + iz;
}

void BinManager::build_atom_neighbors(
    NeighborList& neighbor_list,
    const std::vector<NeighborAtom>& atoms,
    const std::vector<NeighborAtom>& binned_atoms
)
{
    assert(atoms.size() == static_cast<size_t>(neighbor_list.get_nlocal()));

    double sradius2 = sradius_ * sradius_;

    neighbor_list.reset();

    std::vector<int> neigh_tmp;

    const int nlocal = neighbor_list.get_nlocal();
    for (int i = 0; i < nlocal; i++)
    {
        neigh_tmp.clear();
        const NeighborAtom& atom = atoms[i];

        int ix = std::min(
            std::max(int((atom.position_x - x_min_) / bin_sizex_), 0),
            nbinx_ - 1
        );

        int iy = std::min(
            std::max(int((atom.position_y - y_min_) / bin_sizey_), 0),
            nbiny_ - 1
        );

        int iz = std::min(
            std::max(int((atom.position_z - z_min_) / bin_sizez_), 0),
            nbinz_ - 1
        );

        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                for (int dz = -1; dz <= 1; dz++)
                {
                    int jx = ix + dx;
                    int jy = iy + dy;
                    int jz = iz + dz;

                    if (jx < 0 || jx >= nbinx_ ||
                        jy < 0 || jy >= nbiny_ ||
                        jz < 0 || jz >= nbinz_)
                        continue;

                    int nidx = bin_index(jx, jy, jz);

                    for (const ModuleNeighList::LocalAtomIndex binned_atom_index : bins_[nidx].get_atom_indices())
                    {
                        const NeighborAtom& natom = binned_atoms[static_cast<std::size_t>(binned_atom_index)];
                        double dx = atom.position_x - natom.position_x;
                        double dy = atom.position_y - natom.position_y;
                        double dz = atom.position_z - natom.position_z;

                        double dist2 = dx * dx + dy * dy + dz * dz;

                        if (natom.atom_id == atom.atom_id)
                        {
                            continue;
                        }
                        if (dist2 <= sradius2)
                        {
                            neigh_tmp.push_back(natom.atom_id);
                        }
                    }
                }
            }
        }

        const int n = ModuleNeighList::checked_int_size(neigh_tmp.size(), "BinManager neighbor count");

        int* ptr = neighbor_list.allocator_.allocate(n);

        for (int k = 0; k < n; k++)
        {
            assert(ptr != nullptr);
            ptr[k] = neigh_tmp[k];
        }

        neighbor_list.firstneigh_[i] = ptr;
        neighbor_list.numneigh_[i] = n;
    }
}

void BinManager::clear()
{
    for (auto& bin : bins_)
    {
        bin.clear_atoms();
    }

    bins_.clear();
}
