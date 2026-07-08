#include "source_cell/module_neighlist/neighbor_search.h"
#include "source_cell/module_neighlist/domain_decomposition.h"
#include "source_cell/module_neighlist/neighbor_types.h"
#include "source_cell/module_neighlist/unitcell_lite.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
int read_int_arg(int argc, char** argv, int index, int fallback)
{
    return argc <= index ? fallback : std::atoi(argv[index]);
}

double read_double_arg(int argc, char** argv, int index, double fallback)
{
    return argc <= index ? fallback : std::atof(argv[index]);
}

double cell_volume(const ModuleBase::Matrix3& latvec)
{
    const double cx = latvec.e22 * latvec.e33 - latvec.e23 * latvec.e32;
    const double cy = latvec.e23 * latvec.e31 - latvec.e21 * latvec.e33;
    const double cz = latvec.e21 * latvec.e32 - latvec.e22 * latvec.e31;
    return std::abs(latvec.e11 * cx + latvec.e12 * cy + latvec.e13 * cz);
}

ModuleBase::Matrix3 make_simple_lattice_latvec(int nx, int ny, int nz, double spacing, double skew)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = nx * spacing;
    latvec.e12 = 0.0;
    latvec.e13 = 0.0;
    latvec.e21 = skew * ny * spacing;
    latvec.e22 = ny * spacing;
    latvec.e23 = 0.0;
    latvec.e31 = 0.25 * skew * nz * spacing;
    latvec.e32 = 0.5 * skew * nz * spacing;
    latvec.e33 = nz * spacing;
    return latvec;
}

ModuleBase::Vector3<double> direct_to_cartesian(const ModuleBase::Matrix3& latvec,
                                                double fx,
                                                double fy,
                                                double fz)
{
    return ModuleBase::Vector3<double>(fx * latvec.e11 + fy * latvec.e21 + fz * latvec.e31,
                                       fx * latvec.e12 + fy * latvec.e22 + fz * latvec.e32,
                                       fx * latvec.e13 + fy * latvec.e23 + fz * latvec.e33);
}

UnitCellLite make_simple_lattice_ucell(int nx, int ny, int nz, double spacing, double skew)
{
    const ModuleBase::Matrix3 latvec = make_simple_lattice_latvec(nx, ny, nz, spacing, skew);

    std::vector<ModuleBase::Vector3<double>> tau;
    tau.reserve(static_cast<size_t>(nx) * ny * nz);
    for (int ix = 0; ix < nx; ++ix)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int iz = 0; iz < nz; ++iz)
            {
                tau.push_back(direct_to_cartesian(latvec,
                                                  static_cast<double>(ix) / nx,
                                                  static_cast<double>(iy) / ny,
                                                  static_cast<double>(iz) / nz));
            }
        }
    }

    UnitCellLite ucell;
    const double omega = cell_volume(latvec);
    ucell.set_lattice(1.0, omega, latvec);
    ucell.set_atoms(1, {static_cast<int>(tau.size())}, tau);
    return ucell;
}

long long checked_lattice_atom_count(int nx, int ny, int nz)
{
    const long long lx = nx;
    const long long ly = ny;
    const long long lz = nz;
    if (lx > std::numeric_limits<long long>::max() / ly ||
        lx * ly > std::numeric_limits<long long>::max() / lz)
    {
        throw std::overflow_error("benchmark lattice atom count overflows.");
    }
    return lx * ly * lz;
}

long long owner_begin_index(long long n, int coord, int dims)
{
    return (static_cast<long long>(coord) * n + dims - 1) / dims;
}

long long owner_end_index(long long n, int coord, int dims)
{
    return (static_cast<long long>(coord + 1) * n + dims - 1) / dims;
}

void generate_owned_atoms_from_lattice(const DomainDecomposition& decomp,
                                       const ModuleBase::Matrix3& latvec,
                                       int nx,
                                       int ny,
                                       int nz,
                                       std::vector<LocalAtom>& owned_atoms)
{
    owned_atoms.clear();

    const auto& coords = decomp.coords();
    const auto& dims = decomp.dims();

    const long long ix_begin = owner_begin_index(nx, coords[0], dims[0]);
    const long long ix_end = owner_end_index(nx, coords[0], dims[0]);
    const long long iy_begin = owner_begin_index(ny, coords[1], dims[1]);
    const long long iy_end = owner_end_index(ny, coords[1], dims[1]);
    const long long iz_begin = owner_begin_index(nz, coords[2], dims[2]);
    const long long iz_end = owner_end_index(nz, coords[2], dims[2]);

    const std::size_t local_count
        = ModuleNeighList::checked_size_product(
            static_cast<std::size_t>(ix_end - ix_begin),
            ModuleNeighList::checked_size_product(static_cast<std::size_t>(iy_end - iy_begin),
                                                  static_cast<std::size_t>(iz_end - iz_begin),
                                                  "benchmark local atom count"),
            "benchmark local atom count");
    owned_atoms.reserve(local_count);

    for (long long ix = ix_begin; ix < ix_end; ++ix)
    {
        for (long long iy = iy_begin; iy < iy_end; ++iy)
        {
            for (long long iz = iz_begin; iz < iz_end; ++iz)
            {
                const double fx = static_cast<double>(ix) / nx;
                const double fy = static_cast<double>(iy) / ny;
                const double fz = static_cast<double>(iz) / nz;
                const ModuleBase::Vector3<double> frac(fx, fy, fz);
                const ModuleBase::Vector3<double> cart = direct_to_cartesian(latvec, fx, fy, fz);
                const ModuleNeighList::GlobalAtomId global_id
                    = static_cast<ModuleNeighList::GlobalAtomId>((ix * ny + iy) * nz + iz);

                owned_atoms.push_back(LocalAtom(cart,
                                                frac,
                                                0,
                                                0,
                                                global_id,
                                                decomp.rank(),
                                                false));
            }
        }
    }
}

long long count_neighbor_pairs(const NeighborList& list)
{
    long long pairs = 0;
    for (int local_i = 0; local_i < list.get_nlocal(); ++local_i)
    {
        pairs += list.get_numneigh(local_i);
    }
    return pairs;
}

long long square_sum(long long n)
{
    const __int128 value = static_cast<__int128>(n) * (n - 1) * (2 * n - 1) / 6;
    if (value > std::numeric_limits<long long>::max())
    {
        throw std::overflow_error("benchmark square sum exceeds long long range.");
    }
    return static_cast<long long>(value);
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int mpi_rank = 0;
    int mpi_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (argc > 1 && std::string(argv[1]) == "--help")
    {
        if (mpi_rank == 0)
        {
            std::cout << "Usage: neighbor_search_mpi_benchmark [nx ny nz repeat cutoff spacing skew check_serial]\n"
                      << "Defaults: nx=16 ny=16 nz=16 repeat=5 cutoff=1.75 spacing=1.0 skew=0.0 check_serial=1\n";
        }
        MPI_Finalize();
        return 0;
    }

    const int nx = read_int_arg(argc, argv, 1, 16);
    const int ny = read_int_arg(argc, argv, 2, 16);
    const int nz = read_int_arg(argc, argv, 3, 16);
    const int repeat = read_int_arg(argc, argv, 4, 5);
    const double cutoff = read_double_arg(argc, argv, 5, 1.75);
    const double spacing = read_double_arg(argc, argv, 6, 1.0);
    const double skew = read_double_arg(argc, argv, 7, 0.0);
    const int check_serial = read_int_arg(argc, argv, 8, 1);

    if (nx <= 0 || ny <= 0 || nz <= 0 || repeat <= 0 || cutoff <= 0.0 || spacing <= 0.0)
    {
        if (mpi_rank == 0)
        {
            std::cerr << "All dimensions, repeat, cutoff, and spacing must be positive.\n";
        }
        MPI_Finalize();
        return 2;
    }

    const ModuleBase::Matrix3 latvec = make_simple_lattice_latvec(nx, ny, nz, spacing, skew);
    const double lat0 = 1.0;
    const long long nat = checked_lattice_atom_count(nx, ny, nz);

    long long serial_all_atoms = -1;
    long long serial_neighbor_pairs = -1;
    double serial_init_time = 0.0;
    double serial_build_time = 0.0;
    if (mpi_rank == 0 && check_serial)
    {
        UnitCellLite ucell = make_simple_lattice_ucell(nx, ny, nz, spacing, skew);
        NeighborSearch serial;
        const double t0 = MPI_Wtime();
        serial.init(ucell, cutoff);
        const double t1 = MPI_Wtime();
        serial.build_neighbors();
        const double t2 = MPI_Wtime();
        serial_all_atoms = static_cast<long long>(serial.get_all_atoms().size());
        serial_neighbor_pairs = count_neighbor_pairs(serial.get_neighbor_list());
        serial_init_time = t1 - t0;
        serial_build_time = t2 - t1;
    }
    MPI_Bcast(&serial_all_atoms, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&serial_neighbor_pairs, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    double init_time = 0.0;
    double build_time = 0.0;
    double total_time = 0.0;
    long long last_inside = 0;
    long long last_ghost = 0;
    long long last_all = 0;
    long long last_pairs = 0;
    long long inside_index_sum = 0;
    long long inside_index_square_sum = 0;
    int local_failure = 0;

    for (int i = 0; i < repeat; ++i)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        const double t0 = MPI_Wtime();
        DomainDecomposition decomp;
        std::vector<LocalAtom> owned_atoms;
        std::vector<LocalAtom> ghost_atoms;
        NeighborSearch ns;
        decomp.init(MPI_COMM_WORLD, latvec, lat0, cutoff, 0.0);
        generate_owned_atoms_from_lattice(decomp, latvec, nx, ny, nz, owned_atoms);
        decomp.exchange_ghost_atoms(owned_atoms, ghost_atoms);
        ns.init_distributed(owned_atoms, ghost_atoms, cutoff, lat0);
        const double t1 = MPI_Wtime();
        ns.build_neighbors();
        const double t2 = MPI_Wtime();

        init_time += t1 - t0;
        build_time += t2 - t1;
        total_time += t2 - t0;

        if (i == repeat - 1)
        {
            const auto& inside_atoms = ns.get_inside_atoms();
            const auto& ghost_atoms = ns.get_ghost_atoms();
            const auto& all_atoms = ns.get_all_atoms();
            const auto& list = ns.get_neighbor_list();

            last_inside = static_cast<long long>(inside_atoms.size());
            last_ghost = static_cast<long long>(ghost_atoms.size());
            last_all = static_cast<long long>(all_atoms.size());
            last_pairs = 0;
            inside_index_sum = 0;
            inside_index_square_sum = 0;

            for (size_t atom_id = 0; atom_id < all_atoms.size(); ++atom_id)
            {
                if (all_atoms[atom_id].atom_id !=
                    ModuleNeighList::checked_local_atom_index(atom_id, "benchmark atom id"))
                {
                    local_failure = 1;
                }
            }

            for (const NeighborAtom& atom : inside_atoms)
            {
                inside_index_sum += atom.global_id;
                inside_index_square_sum += static_cast<long long>(atom.global_id) * atom.global_id;
            }

            for (int local_i = 0; local_i < list.get_nlocal(); ++local_i)
            {
                last_pairs += list.get_numneigh(local_i);
                for (int ad = 0; ad < list.get_numneigh(local_i); ++ad)
                {
                    const int neighbor_id = list.get_firstneigh(local_i)[ad];
                    if (neighbor_id < 0 || static_cast<size_t>(neighbor_id) >= all_atoms.size())
                    {
                        local_failure = 1;
                    }
                }
            }
        }
    }

    long long global_inside = 0;
    long long global_ghost = 0;
    long long global_all = 0;
    long long global_pairs = 0;
    long long global_index_sum = 0;
    long long global_index_square_sum = 0;
    long long min_all = 0;
    long long max_all = 0;
    long long min_inside = 0;
    long long max_inside = 0;
    long long min_ghost = 0;
    long long max_ghost = 0;
    long long min_pairs = 0;
    long long max_pairs = 0;
    int global_failure = 0;
    MPI_Allreduce(&last_inside, &global_inside, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&last_ghost, &global_ghost, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&last_all, &global_all, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&last_pairs, &global_pairs, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&inside_index_sum, &global_index_sum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&inside_index_square_sum, &global_index_square_sum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&last_all, &min_all, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&last_all, &max_all, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&last_inside, &min_inside, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&last_inside, &max_inside, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&last_ghost, &min_ghost, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&last_ghost, &max_ghost, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&last_pairs, &min_pairs, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&last_pairs, &max_pairs, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    double max_init_time = 0.0;
    double max_build_time = 0.0;
    double max_total_time = 0.0;
    MPI_Reduce(&init_time, &max_init_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&build_time, &max_build_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    const bool ownership_ok = global_inside == nat &&
                              global_index_sum == nat * (nat - 1) / 2 &&
                              global_index_square_sum == square_sum(nat);
    const bool neighbor_pairs_ok = !check_serial || global_pairs == serial_neighbor_pairs;
    const bool all_ok = ownership_ok && global_failure == 0 && neighbor_pairs_ok;

    if (mpi_rank == 0)
    {
        std::cout << "NeighborSearch MPI halo benchmark\n"
                  << "algorithm fractional_halo_bins\n"
                  << "np " << mpi_size << "\n"
                  << "atoms " << nat << "\n"
                  << "grid " << nx << " " << ny << " " << nz << "\n"
                  << "repeat " << repeat << "\n"
                  << "cutoff " << cutoff << "\n"
                  << "spacing " << spacing << "\n"
                  << "skew " << skew << "\n"
                  << "check_serial " << check_serial << "\n"
                  << "serial_all_atoms " << serial_all_atoms << "\n"
                  << "serial_neighbor_pairs " << serial_neighbor_pairs << "\n"
                  << "inside_sum " << global_inside << "\n"
                  << "inside_min " << min_inside << "\n"
                  << "inside_max " << max_inside << "\n"
                  << "ghost_sum " << global_ghost << "\n"
                  << "ghost_min " << min_ghost << "\n"
                  << "ghost_max " << max_ghost << "\n"
                  << "all_atoms_sum " << global_all << "\n"
                  << "all_atoms_min " << min_all << "\n"
                  << "all_atoms_max " << max_all << "\n"
                  << "neighbor_pairs_sum " << global_pairs << "\n"
                  << "neighbor_pairs_min " << min_pairs << "\n"
                  << "neighbor_pairs_max " << max_pairs << "\n"
                  << "time_serial_ref_init " << serial_init_time << "\n"
                  << "time_serial_ref_build " << serial_build_time << "\n"
                  << "time_serial_ref_total " << serial_init_time + serial_build_time << "\n"
                  << "time_init_max_total " << max_init_time << "\n"
                  << "time_build_max_total " << max_build_time << "\n"
                  << "time_total_max_total " << max_total_time << "\n"
                  << "time_init_max_avg " << max_init_time / repeat << "\n"
                  << "time_build_max_avg " << max_build_time / repeat << "\n"
                  << "time_total_max_avg " << max_total_time / repeat << "\n"
                  << "ownership_ok " << (ownership_ok ? 1 : 0) << "\n"
                  << "neighbor_pairs_ok " << (neighbor_pairs_ok ? 1 : 0) << "\n"
                  << "neighbor_ids_ok " << (global_failure == 0 ? 1 : 0) << "\n";
    }

    MPI_Finalize();
    return all_ok ? 0 : 1;
}
