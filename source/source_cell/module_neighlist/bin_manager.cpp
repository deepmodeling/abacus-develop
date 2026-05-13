#include <limits>
#include <cmath>
#include <algorithm>
#include "bin_manager.h"



void BinManager::init_bins(
    double sr,
    const std::vector<NeighborAtom>& inside_atoms,
    const std::vector<NeighborAtom>& ghost_atoms
)
{
    sradius = sr;
    if(inside_atoms.empty() && ghost_atoms.empty())
    {
        x_min=y_min=z_min=0;
        x_max=y_max=z_max=0;
        nbinx=nbiny=nbinz=1;
        bins.clear();
        bins.resize(1);
        return;
    }

    x_min = y_min = z_min = std::numeric_limits<double>::max();

    x_max = y_max = z_max = std::numeric_limits<double>::lowest();
    

    auto update_bounds = [&](const std::vector<NeighborAtom>& atoms)
    {
        for (const auto& atom : atoms)
        {
            x_min = std::min(x_min, atom.position_x);
            x_max = std::max(x_max, atom.position_x);

            y_min = std::min(y_min, atom.position_y);
            y_max = std::max(y_max, atom.position_y);

            z_min = std::min(z_min, atom.position_z);
            z_max = std::max(z_max, atom.position_z);
        }
    };

    update_bounds(inside_atoms);
    update_bounds(ghost_atoms);

    bin_sizex = bin_sizey = bin_sizez = sradius;

    nbinx = std::ceil((x_max - x_min) / bin_sizex);
    nbiny = std::ceil((y_max - y_min) / bin_sizey);
    nbinz = std::ceil((z_max - z_min) / bin_sizez);

    nbinx = std::max(1, nbinx);
    nbiny = std::max(1, nbiny);
    nbinz = std::max(1, nbinz);

    int nbins = nbinx * nbiny * nbinz;

    bins.clear();

    bins.resize(nbins);

    for (int ix = 0; ix < nbinx; ++ix)
    {
        for (int iy = 0; iy < nbiny; ++iy)
        {
            for (int iz = 0; iz < nbinz; ++iz)
            {
                int idx = ix * nbiny * nbinz + iy * nbinz + iz;

                bins[idx].id_x = ix;
                bins[idx].id_y = iy;
                bins[idx].id_z = iz;

                bins[idx].atoms.clear();
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
            std::max(int((atom.position_x - x_min) / bin_sizex), 0),
            nbinx - 1
        );

        int iy = std::min(
            std::max(int((atom.position_y - y_min) / bin_sizey), 0),
            nbiny - 1
        );

        int iz = std::min(
            std::max(int((atom.position_z - z_min) / bin_sizez), 0),
            nbinz - 1
        );

        int idx = ix * nbiny * nbinz + iy * nbinz + iz;

        bins[idx].atoms.push_back(atom);
    };

    for (const auto& atom : inside_atoms) bin_atom(atom);

    for (const auto& atom : ghost_atoms) bin_atom(atom);
}

void BinManager::build_atom_neighbors(
    NeighborList& neighbor_list,
    std::vector<NeighborAtom>& atoms
)
{
    assert(atoms.size() == neighbor_list.numneigh.size());

    double sradius2 = sradius * sradius;

    neighbor_list.reset();

    for (int i = 0; i < atoms.size(); i++)
    {
        std::vector<int> neigh_tmp;

        int ix = std::min(
            std::max(int((atoms[i].position_x - x_min) / bin_sizex), 0),
            nbinx - 1
        );

        int iy = std::min(
            std::max(int((atoms[i].position_y - y_min) / bin_sizey), 0),
            nbiny - 1
        );

        int iz = std::min(
            std::max(int((atoms[i].position_z - z_min) / bin_sizez), 0),
            nbinz - 1
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

                    if (jx < 0 || jx >= nbinx ||
                        jy < 0 || jy >= nbiny ||
                        jz < 0 || jz >= nbinz)
                        continue;

                    int nidx = jx * nbiny * nbinz + jy * nbinz + jz;

                    for (const NeighborAtom& natom : bins[nidx].atoms)
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

        //std::cout<<n<<std::endl;

        int* ptr = neighbor_list.allocator.allocate(n);
    
        for (int k = 0; k < n; k++)
        {
            assert(ptr != nullptr);
            ptr[k] = neigh_tmp[k];
        }

        neighbor_list.firstneigh[i] = ptr;
        neighbor_list.numneigh[i] = n;
    }
}

void BinManager::clear()
{
    for (auto& bin : bins)
    {
        bin.atoms.clear();
    }

    bins.clear();
}