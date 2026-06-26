#include "../sltk_grid.h"
#include "prepare_unitcell.h"

#ifndef NEIGHBOR_BENCHMARK_BETA1
#include "../atom_pack.h"
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
#include <sys/resource.h>
#include <tuple>
#include <valarray>
#include <vector>

#ifdef __LCAO
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
LCAO_Orbitals::LCAO_Orbitals() {}
LCAO_Orbitals::~LCAO_Orbitals() {}
#endif
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
    this->start_mag = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_mag;
}

namespace
{

using PairKey = std::tuple<int, int, int, int, int, int, int>;

struct Options
{
    int atoms = 4096;
    int repeats = 10;
    double spacing = 1.5;
    double radius = 2.2;
    std::string cell = "orthogonal";
#ifndef NEIGHBOR_BENCHMARK_BETA1
    std::string mode = "half14";
#endif
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
    for (int index = 1; index < argc; ++index)
    {
        const std::string argument = argv[index];
        if (index + 1 >= argc)
        {
            throw std::invalid_argument("Every benchmark option requires a value.");
        }
        if (argument == "--atoms")
        {
            options.atoms = parse_positive_int(argv[++index], "--atoms");
        }
        else if (argument == "--repeats")
        {
            options.repeats = parse_positive_int(argv[++index], "--repeats");
        }
        else if (argument == "--spacing")
        {
            options.spacing = parse_positive_double(argv[++index], "--spacing");
        }
        else if (argument == "--radius")
        {
            options.radius = parse_positive_double(argv[++index], "--radius");
        }
        else if (argument == "--cell")
        {
            options.cell = argv[++index];
        }
#ifndef NEIGHBOR_BENCHMARK_BETA1
        else if (argument == "--mode")
        {
            options.mode = argv[++index];
        }
#endif
        else
        {
            throw std::invalid_argument("Unknown benchmark option: " + argument);
        }
    }
    if (options.cell != "orthogonal" && options.cell != "triclinic")
    {
        throw std::invalid_argument("--cell must be orthogonal or triclinic.");
    }
#ifndef NEIGHBOR_BENCHMARK_BETA1
    if (options.mode != "half14" && options.mode != "full27")
    {
        throw std::invalid_argument("--mode must be half14 or full27.");
    }
#endif
    return options;
}

std::array<int, 3> shape_for_atoms(const int atoms)
{
    if (atoms == 4096)
    {
        return std::array<int, 3>{{16, 16, 16}};
    }
    if (atoms == 32768)
    {
        return std::array<int, 3>{{32, 32, 32}};
    }
    if (atoms == 131072)
    {
        return std::array<int, 3>{{64, 64, 32}};
    }
    const int side = static_cast<int>(std::llround(std::cbrt(static_cast<double>(atoms))));
    if (side <= 0 || side * side * side != atoms)
    {
        throw std::invalid_argument("--atoms must be a perfect cube or one of 4096, 32768 and 131072.");
    }
    return std::array<int, 3>{{side, side, side}};
}

UcellTestPrepare make_unitcell_input(const Options& options)
{
    const std::array<int, 3> shape = shape_for_atoms(options.atoms);
    const double lx = shape[0] * options.spacing;
    const double ly = shape[1] * options.spacing;
    const double lz = shape[2] * options.spacing;
    std::valarray<double> lattice{
        lx, 0.0, 0.0,
        0.0, ly, 0.0,
        0.0, 0.0, lz};
    if (options.cell == "triclinic")
    {
        lattice[3] = 0.20 * ly;
        lattice[6] = 0.10 * lz;
        lattice[7] = 0.15 * lz;
    }

    std::valarray<double> coordinates(3 * options.atoms);
    int atom = 0;
    for (int ix = 0; ix < shape[0]; ++ix)
    {
        for (int iy = 0; iy < shape[1]; ++iy)
        {
            for (int iz = 0; iz < shape[2]; ++iz)
            {
                const double fx = (ix + 0.5) / shape[0];
                const double fy = (iy + 0.5) / shape[1];
                const double fz = (iz + 0.5) / shape[2];
                coordinates[3 * atom] = fx * lattice[0] + fy * lattice[3] + fz * lattice[6];
                coordinates[3 * atom + 1] = fx * lattice[1] + fy * lattice[4] + fz * lattice[7];
                coordinates[3 * atom + 2] = fx * lattice[2] + fy * lattice[5] + fz * lattice[8];
                ++atom;
            }
        }
    }

    return UcellTestPrepare("benchmark",
                            0,
                            false,
                            false,
                            false,
                            "None",
                            1.0,
                            lattice,
                            std::vector<std::string>{"X"},
                            std::vector<std::string>{"X.upf"},
                            std::vector<std::string>{"auto"},
                            std::vector<std::string>{"X.orb"},
                            std::valarray<int>{options.atoms},
                            std::vector<double>{1.0},
                            "Cartesian",
                            coordinates);
}

std::vector<PairKey> collect_legacy_pairs(const Grid& grid)
{
    std::vector<PairKey> keys;
    for (int type = 0; type < static_cast<int>(grid.all_adj_info.size()); ++type)
    {
        for (int atom = 0; atom < static_cast<int>(grid.all_adj_info[type].size()); ++atom)
        {
            for (const FAtom* neighbor: grid.all_adj_info[type][atom])
            {
                keys.push_back(PairKey(type,
                                       atom,
                                       neighbor->type,
                                       neighbor->natom,
                                       neighbor->cell_x,
                                       neighbor->cell_y,
                                       neighbor->cell_z));
            }
        }
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

#ifndef NEIGHBOR_BENCHMARK_BETA1
std::vector<PairKey> collect_current_pairs(const std::vector<ModuleNeighbor::NeighborPair>& pairs)
{
    std::vector<PairKey> keys;
    keys.reserve(pairs.size());
    for (const ModuleNeighbor::NeighborPair& pair: pairs)
    {
        keys.push_back(PairKey(pair.center_type,
                               pair.center_natom,
                               pair.neighbor_type,
                               pair.neighbor_natom,
                               pair.cell_x,
                               pair.cell_y,
                               pair.cell_z));
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}
#endif

unsigned long long pair_hash(const std::vector<PairKey>& pairs)
{
    unsigned long long hash = 1469598103934665603ULL;
    for (const PairKey& pair: pairs)
    {
        const int values[] = {std::get<0>(pair),
                              std::get<1>(pair),
                              std::get<2>(pair),
                              std::get<3>(pair),
                              std::get<4>(pair),
                              std::get<5>(pair),
                              std::get<6>(pair)};
        for (const int value: values)
        {
            hash ^= static_cast<unsigned int>(value);
            hash *= 1099511628211ULL;
        }
    }
    return hash;
}

long long elapsed_microseconds(const std::function<void()>& operation)
{
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    operation();
    return std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - start)
        .count();
}

long long peak_rss_kib()
{
    struct rusage usage
    {
    };
    return getrusage(RUSAGE_SELF, &usage) == 0 ? usage.ru_maxrss : -1;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = parse_options(argc, argv);
        UcellTestPrepare input = make_unitcell_input(options);
        UnitCell* unitcell = input.SetUcellInfo();
        const std::string implementation
#ifdef NEIGHBOR_BENCHMARK_BETA1
            = "beta1_legacy";
#else
            = options.mode == "full27" ? "current_full27" : "current_half14";
#endif

        for (int repeat = -1; repeat < options.repeats; ++repeat)
        {
            long long pack_us = -1;
            long long storage_us = -1;
            long long search_us = -1;
            long long total_us = 0;
            std::vector<PairKey> pairs;

#ifdef NEIGHBOR_BENCHMARK_BETA1
            std::ofstream output("/tmp/beta1_grid_benchmark.out");
            Grid grid(0);
            total_us = elapsed_microseconds(
                [&]() { grid.init(output, *unitcell, options.radius, true); });
            output.close();
            pairs = collect_legacy_pairs(grid);
#else
            ModuleNeighbor::AtomPack pack;
            ModuleNeighbor::GridStorage storage;
            std::vector<ModuleNeighbor::NeighborPair> neighbor_pairs;
            total_us = elapsed_microseconds([&]() {
                pack_us = elapsed_microseconds([&]() {
                    pack = ModuleNeighbor::build_atom_pack_from_unitcell(*unitcell,
                                                                          true,
                                                                          options.radius);
                });
                storage_us = elapsed_microseconds([&]() {
                    storage = ModuleNeighbor::build_grid_storage_from_atom_pack(
                        pack,
                        options.radius + 0.1);
                });
                search_us = elapsed_microseconds([&]() {
                    neighbor_pairs
                        = options.mode == "full27"
                              ? ModuleNeighbor::build_neighbor_pairs_27(pack,
                                                                        storage,
                                                                        options.radius)
                              : ModuleNeighbor::build_neighbor_pairs_14(pack,
                                                                        storage,
                                                                        options.radius);
                });
            });
            pairs = collect_current_pairs(neighbor_pairs);
#endif
            if (repeat >= 0)
            {
                std::cout << "neighbor_grid_path," << implementation << ","
                          << options.cell << "," << options.atoms << "," << options.radius
                          << "," << repeat << "," << pack_us << "," << storage_us << ","
                          << search_us << "," << total_us << "," << pairs.size() << ","
                          << peak_rss_kib() << "," << pair_hash(pairs) << "\n";
            }
        }
        delete unitcell;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "grid path benchmark failed: " << error.what() << "\n";
        return 1;
    }
}
