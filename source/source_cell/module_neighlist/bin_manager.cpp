#include <limits>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "bin_manager.h"

// ========== Bin class implementation ==========

int Bin::get_id_x() const {
    return id_x_;
}

int Bin::get_id_y() const {
    return id_y_;
}

int Bin::get_id_z() const {
    return id_z_;
}

const std::vector<NeighborAtom>& Bin::get_atoms() const {
    return atoms_;
}

void Bin::set_id(int ix, int iy, int iz) {
    id_x_ = ix;
    id_y_ = iy;
    id_z_ = iz;
}

void Bin::clear_atoms() {
    atoms_.clear();
}

void Bin::add_atom(const NeighborAtom& atom) {
    atoms_.push_back(atom);
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
    return static_cast<int>(bins_.size());
}

int BinManager::get_bin_atom_count(int bin_index) const {
    if (bin_index < 0 || bin_index >= static_cast<int>(bins_.size())) {
        return 0;
    }
    return static_cast<int>(bins_[bin_index].get_atoms().size());
}

// ========== BinManager main methods ==========

void BinManager::init_bins(
    double sr,
    const std::vector<NeighborAtom>& inside_atoms,
    const std::vector<NeighborAtom>& ghost_atoms
)
{
    sradius_ = sr;
    if(inside_atoms.empty() && ghost_atoms.empty())
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

    update_bounds(inside_atoms);
    update_bounds(ghost_atoms);

    bin_sizex_ = bin_sizey_ = bin_sizez_ = sradius_;

    nbinx_ = std::ceil((x_max_ - x_min_) / bin_sizex_);
    nbiny_ = std::ceil((y_max_ - y_min_) / bin_sizey_);
    nbinz_ = std::ceil((z_max_ - z_min_) / bin_sizez_);

    nbinx_ = std::max(1, nbinx_);
    nbiny_ = std::max(1, nbiny_);
    nbinz_ = std::max(1, nbinz_);

    int nbins = nbinx_ * nbiny_ * nbinz_;

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
    const std::vector<NeighborAtom>& inside_atoms,
    const std::vector<NeighborAtom>& ghost_atoms
)
{
    auto bin_atom = [&](const NeighborAtom& atom)
    {
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

        bins_[idx].add_atom(atom);
    };

    for (const auto& atom : inside_atoms) bin_atom(atom);
    for (const auto& atom : ghost_atoms) bin_atom(atom);
}

int BinManager::bin_index(int ix, int iy, int iz) const {
    return ix * nbiny_ * nbinz_ + iy * nbinz_ + iz;
}

void BinManager::build_atom_neighbors(
    NeighborList& neighbor_list,
    std::vector<NeighborAtom>& atoms
)
{
    assert(atoms.size() == static_cast<size_t>(neighbor_list.get_nlocal()));

    double sradius2 = sradius_ * sradius_;

    neighbor_list.reset();

    std::vector<int> neigh_tmp;

    for (int i = 0; i < atoms.size(); i++)
    {
        neigh_tmp.clear();

        int ix = std::min(
            std::max(int((atoms[i].position_x - x_min_) / bin_sizex_), 0),
            nbinx_ - 1
        );

        int iy = std::min(
            std::max(int((atoms[i].position_y - y_min_) / bin_sizey_), 0),
            nbiny_ - 1
        );

        int iz = std::min(
            std::max(int((atoms[i].position_z - z_min_) / bin_sizez_), 0),
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

                    for (const NeighborAtom& natom : bins_[nidx].get_atoms())
                    {
                        double dx = atoms[i].position_x - natom.position_x;
                        double dy = atoms[i].position_y - natom.position_y;
                        double dz = atoms[i].position_z - natom.position_z;

                        double dist2 = dx * dx + dy * dy + dz * dz;

                        if (dist2 <= sradius2 && dist2 != 0)
                        {
                            neigh_tmp.push_back(natom.atom_id);
                        }
                    }
                }
            }
        }

        int n = neigh_tmp.size();

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