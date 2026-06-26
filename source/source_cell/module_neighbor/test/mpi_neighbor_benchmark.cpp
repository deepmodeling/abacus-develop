#include "../atom_pack.h"
#include "../mpi_domain.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct Options
{
    int atoms_per_rank = 4096;
    int repeats = 5;
    int shape_x = 0;
    int shape_y = 0;
    int shape_z = 0;
    double spacing = 1.5;
    double radius = 2.2;
    std::string cell = "orthogonal";
};

int parse_positive_int(const char* value, const char* name)
{
    const long parsed = std::strtol(value, nullptr, 10);
    if (parsed <= 0 || parsed > std::numeric_limits<int>::max())
    {
        throw std::invalid_argument(std::string(name) + " must be a positive integer.");
    }
    return static_cast<int>(parsed);
}

double parse_positive_double(const char* value, const char* name)
{
    const double parsed = std::strtod(value, nullptr);
    if (!std::isfinite(parsed) || parsed <= 0.0)
    {
        throw std::invalid_argument(std::string(name) + " must be finite and positive.");
    }
    return parsed;
}

Options parse_options(const int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument = argv[i];
        if (i + 1 >= argc)
        {
            throw std::invalid_argument("Every benchmark option requires a value.");
        }
        if (argument == "--atoms-per-rank")
        {
            options.atoms_per_rank = parse_positive_int(argv[++i], "--atoms-per-rank");
        }
        else if (argument == "--repeats")
        {
            options.repeats = parse_positive_int(argv[++i], "--repeats");
        }
        else if (argument == "--spacing")
        {
            options.spacing = parse_positive_double(argv[++i], "--spacing");
        }
        else if (argument == "--radius")
        {
            options.radius = parse_positive_double(argv[++i], "--radius");
        }
        else if (argument == "--shape-x")
        {
            options.shape_x = parse_positive_int(argv[++i], "--shape-x");
        }
        else if (argument == "--shape-y")
        {
            options.shape_y = parse_positive_int(argv[++i], "--shape-y");
        }
        else if (argument == "--shape-z")
        {
            options.shape_z = parse_positive_int(argv[++i], "--shape-z");
        }
        else if (argument == "--cell")
        {
            options.cell = argv[++i];
        }
        else
        {
            throw std::invalid_argument("Unknown benchmark option: " + argument);
        }
    }
    const bool any_shape = options.shape_x > 0 || options.shape_y > 0 || options.shape_z > 0;
    const bool complete_shape = options.shape_x > 0 && options.shape_y > 0 && options.shape_z > 0;
    if (any_shape != complete_shape)
    {
        throw std::invalid_argument("--shape-x, --shape-y and --shape-z must be supplied together.");
    }
    if (options.cell != "orthogonal" && options.cell != "triclinic")
    {
        throw std::invalid_argument("--cell must be orthogonal or triclinic.");
    }
    return options;
}

std::array<double, 3> fractional_to_cartesian(
    const std::array<std::array<double, 3>, 3>& lattice,
    const double fx,
    const double fy,
    const double fz)
{
    return std::array<double, 3>{{fx * lattice[0][0] + fy * lattice[1][0] + fz * lattice[2][0],
                                  fx * lattice[0][1] + fy * lattice[1][1] + fz * lattice[2][1],
                                  fx * lattice[0][2] + fy * lattice[1][2] + fz * lattice[2][2]}};
}

std::vector<ModuleNeighbor::MpiAtomRecord> make_global_atoms(
    const int nx,
    const int ny,
    const int nz,
    const std::array<std::array<double, 3>, 3>& lattice)
{
    std::vector<ModuleNeighbor::MpiAtomRecord> records;
    records.reserve(nx * ny * nz);
    int index = 0;
    for (int ix = 0; ix < nx; ++ix)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int iz = 0; iz < nz; ++iz)
            {
                const std::array<double, 3> position
                    = fractional_to_cartesian(lattice,
                                              (ix + 0.5) / nx,
                                              (iy + 0.5) / ny,
                                              (iz + 0.5) / nz);
                records.emplace_back(position[0],
                                     position[1],
                                     position[2],
                                     index,
                                     std::array<int, 3>{{0, 0, 0}},
                                     false,
                                     0,
                                     index);
                ++index;
            }
        }
    }
    return records;
}

long long elapsed_microseconds(const std::function<void()>& operation)
{
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    operation();
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
}

} // namespace

int main(int argc, char** argv)
{
#ifndef __MPI
    std::cerr << "MPI benchmark requires ENABLE_MPI=ON.\n";
    return 1;
#else
    MPI_Init(&argc, &argv);
    int world_rank = 0;
    int world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    try
    {
        const Options options = parse_options(argc, argv);
        int rank_dims[3] = {0, 0, 0};
        MPI_Dims_create(world_size, 3, rank_dims);

        int nx = options.shape_x;
        int ny = options.shape_y;
        int nz = options.shape_z;
        const bool strong_scaling = nx > 0;
        if (!strong_scaling)
        {
            const int local_side = static_cast<int>(
                std::llround(std::cbrt(static_cast<double>(options.atoms_per_rank))));
            if (local_side * local_side * local_side != options.atoms_per_rank)
            {
                throw std::invalid_argument("Weak-scaling atoms per rank must be a perfect cube.");
            }
            nx = local_side * rank_dims[0];
            ny = local_side * rank_dims[1];
            nz = local_side * rank_dims[2];
        }
        if (nx <= 0 || ny <= 0 || nz <= 0
            || nx > std::numeric_limits<int>::max() / ny
            || nx * ny > std::numeric_limits<int>::max() / nz)
        {
            throw std::overflow_error("Global benchmark atom count exceeds the int range.");
        }
        const int global_atoms = nx * ny * nz;
        const double lx = nx * options.spacing;
        const double ly = ny * options.spacing;
        const double lz = nz * options.spacing;
        std::array<std::array<double, 3>, 3> lattice{{
            std::array<double, 3>{{lx, 0.0, 0.0}},
            std::array<double, 3>{{0.0, ly, 0.0}},
            std::array<double, 3>{{0.0, 0.0, lz}},
        }};
        if (options.cell == "triclinic")
        {
            lattice[1][0] = 0.20 * ly;
            lattice[2][0] = 0.10 * lz;
            lattice[2][1] = 0.15 * lz;
        }
        const std::vector<ModuleNeighbor::MpiAtomRecord> all_atoms
            = make_global_atoms(nx, ny, nz, lattice);

        ModuleNeighbor::MpiDomain domain;
        domain.initialize_lattice(MPI_COMM_WORLD,
                                  std::array<double, 3>{{0.0, 0.0, 0.0}},
                                  lattice,
                                  options.radius,
                                  true);
        const std::vector<int> local_indices = domain.select_local_atoms(all_atoms);
        std::vector<ModuleNeighbor::MpiAtomRecord> local_atoms;
        local_atoms.reserve(local_indices.size());
        for (const int index: local_indices)
        {
            local_atoms.push_back(all_atoms[index]);
        }

        for (int repeat = 0; repeat < options.repeats; ++repeat)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            ModuleNeighbor::MpiGhostExchangeStats stats;
            std::vector<ModuleNeighbor::MpiAtomRecord> ghosts;
            const long long exchange_us
                = elapsed_microseconds([&]() { ghosts = domain.exchange_ghost_atoms(local_atoms, &stats); });

            ModuleNeighbor::AtomPack pack;
            pack.reserve(static_cast<int>(local_atoms.size() + ghosts.size()));
            for (const ModuleNeighbor::MpiAtomRecord& record: local_atoms)
            {
                pack.append_mpi_record(record, record.type, record.natom);
            }
            for (const ModuleNeighbor::MpiAtomRecord& record: ghosts)
            {
                pack.append_mpi_record(record, record.type, record.natom);
            }

            ModuleNeighbor::GridStorage storage;
            std::vector<ModuleNeighbor::NeighborPair> pairs;
            const long long search_us = elapsed_microseconds([&]() {
                storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, options.radius + 0.1);
                pairs = ModuleNeighbor::build_neighbor_pairs_14(pack, storage, options.radius);
            });

            long long local_values[] = {exchange_us,
                                        search_us,
                                        static_cast<long long>(local_atoms.size()),
                                        static_cast<long long>(ghosts.size()),
                                        stats.neighbor_rank_count,
                                        stats.nonempty_send_rank_count,
                                        stats.sent_atom_count,
                                        stats.received_ghost_count,
                                        stats.sent_payload_count,
                                        stats.received_payload_count,
                                        static_cast<long long>(pairs.size())};
            long long maximum[11] = {};
            long long sum[11] = {};
            MPI_Reduce(local_values, maximum, 11, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(local_values, sum, 11, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            if (world_rank == 0)
            {
                std::cout << "neighbor_mpi," << (strong_scaling ? "strong" : "weak") << ","
                          << options.cell << "," << world_size << ","
                          << static_cast<double>(global_atoms) / world_size << ","
                          << global_atoms << "," << repeat << "," << options.radius;
                for (const long long value: maximum)
                {
                    std::cout << "," << value;
                }
                for (const long long value: sum)
                {
                    std::cout << "," << value;
                }
                std::cout << "\n";
            }
        }
        MPI_Finalize();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "MPI neighbor benchmark failed on rank " << world_rank << ": " << error.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
#endif
}
