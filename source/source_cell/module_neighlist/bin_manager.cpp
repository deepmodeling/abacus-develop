#include <limits>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "bin_manager.h"
#include "source_base/timer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
constexpr int neighbor_build_openmp_threshold = 256;
}

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
    ModuleBase::timer::start("BinManager", "init_bins");
    sradius_ = sr;
    if(inside_atoms.empty() && ghost_atoms.empty())
    {
        x_min_ = y_min_ = z_min_ = 0;
        x_max_ = y_max_ = z_max_ = 0;
        nbinx_ = nbiny_ = nbinz_ = 1;
        bins_.clear();
        bins_.resize(1);
        ModuleBase::timer::end("BinManager", "init_bins");
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
    ModuleBase::timer::end("BinManager", "init_bins");
}

void BinManager::do_binning(
    const std::vector<NeighborAtom>& inside_atoms,
    const std::vector<NeighborAtom>& ghost_atoms
)
{
    ModuleBase::timer::start("BinManager", "do_binning");
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
    ModuleBase::timer::end("BinManager", "do_binning");
}

int BinManager::bin_index(int ix, int iy, int iz) const {
    return ix * nbiny_ * nbinz_ + iy * nbinz_ + iz;
}

template <typename Emit>
void BinManager::visit_neighbors(
    int i,
    const std::vector<NeighborAtom>& atoms,
    double sradius2,
    const Emit& emit
) const
{
    const int ix = std::min(
        std::max(int((atoms[i].position_x - x_min_) / bin_sizex_), 0),
        nbinx_ - 1
    );

    const int iy = std::min(
        std::max(int((atoms[i].position_y - y_min_) / bin_sizey_), 0),
        nbiny_ - 1
    );

    const int iz = std::min(
        std::max(int((atoms[i].position_z - z_min_) / bin_sizez_), 0),
        nbinz_ - 1
    );

    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            for (int dz = -1; dz <= 1; dz++)
            {
                const int jx = ix + dx;
                const int jy = iy + dy;
                const int jz = iz + dz;

                if (jx < 0 || jx >= nbinx_ ||
                    jy < 0 || jy >= nbiny_ ||
                    jz < 0 || jz >= nbinz_)
                {
                    continue;
                }

                const int nidx = bin_index(jx, jy, jz);

                for (const NeighborAtom& natom : bins_[nidx].get_atoms())
                {
                    const double dx = atoms[i].position_x - natom.position_x;
                    const double dy = atoms[i].position_y - natom.position_y;
                    const double dz = atoms[i].position_z - natom.position_z;

                    const double dist2 = dx * dx + dy * dy + dz * dz;

                    if (dist2 <= sradius2 && dist2 != 0)
                    {
                        emit(natom.atom_id);
                    }
                }
            }
        }
    }
}

void BinManager::build_atom_neighbors(
    NeighborList& neighbor_list,
    std::vector<NeighborAtom>& atoms
)
{
    ModuleBase::timer::start("BinManager", "build_atom_neighbors");
    assert(atoms.size() == static_cast<size_t>(neighbor_list.get_nlocal()));

    const int nlocal = static_cast<int>(atoms.size());
    const double sradius2 = sradius_ * sradius_;

    neighbor_list.reset();

#ifdef _OPENMP
    const bool use_parallel = nlocal >= neighbor_build_openmp_threshold && omp_get_max_threads() > 1;
    if (use_parallel)
    {
        std::vector<int> neighbor_counts(nlocal, 0);

#pragma omp parallel for schedule(static)
        for (int i = 0; i < nlocal; i++)
        {
            int count = 0;
            visit_neighbors(i, atoms, sradius2, [&](const int) { ++count; });
            neighbor_counts[i] = count;
        }

        for (int i = 0; i < nlocal; i++)
        {
            const int n = neighbor_counts[i];
            neighbor_list.firstneigh_[i] = neighbor_list.allocator_.allocate(n);
            neighbor_list.numneigh_[i] = n;
        }

#pragma omp parallel for schedule(static)
        for (int i = 0; i < nlocal; i++)
        {
            int* ptr = neighbor_list.firstneigh_[i];
            int k = 0;
            visit_neighbors(i, atoms, sradius2, [&](const int atom_id)
            {
                assert(ptr != nullptr);
                ptr[k++] = atom_id;
            });
            assert(k == neighbor_counts[i]);
        }

        ModuleBase::timer::end("BinManager", "build_atom_neighbors");
        return;
    }
#endif

    std::vector<int> neigh_tmp;

    for (int i = 0; i < nlocal; i++)
    {
        neigh_tmp.clear();
        visit_neighbors(i, atoms, sradius2, [&](const int atom_id) { neigh_tmp.push_back(atom_id); });

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
    ModuleBase::timer::end("BinManager", "build_atom_neighbors");
}

void BinManager::clear()
{
    for (auto& bin : bins_)
    {
        bin.clear_atoms();
    }

    bins_.clear();
}
