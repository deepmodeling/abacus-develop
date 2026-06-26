#include "../atom_pack.h"

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
#include <vector>

namespace
{

struct Options
{
    int atoms = 4096;
    int repeats = 10;
    double spacing = 1.5;
    double radius = 2.2;
    std::string mode = "half14";
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
        if (argument == "--atoms")
        {
            options.atoms = parse_positive_int(argv[++i], "--atoms");
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
        else if (argument == "--mode")
        {
            options.mode = argv[++i];
        }
        else
        {
            throw std::invalid_argument("Unknown benchmark option: " + argument);
        }
    }
    if (options.mode != "half14" && options.mode != "full27" && options.mode != "half14_ref")
    {
        throw std::invalid_argument("--mode must be half14, full27 or half14_ref.");
    }
    return options;
}

ModuleNeighbor::AtomPack make_pack(const Options& options)
{
    ModuleNeighbor::AtomPack pack;
    pack.reserve(options.atoms);
    const int side = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(options.atoms))));
    for (int index = 0; index < options.atoms; ++index)
    {
        const int ix = index % side;
        const int iy = (index / side) % side;
        const int iz = index / (side * side);
        pack.append_atom(ix * options.spacing,
                         iy * options.spacing,
                         iz * options.spacing,
                         0,
                         index,
                         0,
                         0,
                         0,
                         index,
                         false);
    }
    return pack;
}

unsigned long long pair_hash(const std::vector<ModuleNeighbor::NeighborPair>& pairs)
{
    unsigned long long hash = 1469598103934665603ULL;
    for (const ModuleNeighbor::NeighborPair& pair: pairs)
    {
        const int values[] = {pair.center_type,
                              pair.center_natom,
                              pair.neighbor_type,
                              pair.neighbor_natom,
                              pair.cell_x,
                              pair.cell_y,
                              pair.cell_z};
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
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
}

long long peak_rss_kib()
{
    struct rusage usage
    {
    };
    if (getrusage(RUSAGE_SELF, &usage) != 0)
    {
        return -1;
    }
    return usage.ru_maxrss;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = parse_options(argc, argv);
        const ModuleNeighbor::AtomPack pack = make_pack(options);
        const double box_edge = options.radius + 0.1;

        // One warm-up removes one-time allocator and page-fault noise from the
        // reported repetitions.
        ModuleNeighbor::GridStorage warm_storage
            = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, box_edge);
        if (options.mode == "full27")
        {
            ModuleNeighbor::build_neighbor_pairs_27(pack, warm_storage, options.radius);
        }
        else
        {
            ModuleNeighbor::build_neighbor_pairs_14(pack, warm_storage, options.radius);
        }

        for (int repeat = 0; repeat < options.repeats; ++repeat)
        {
            ModuleNeighbor::GridStorage storage;
            std::vector<ModuleNeighbor::NeighborPair> pairs;
            std::vector<ModuleNeighbor::NeighborPair> reference;
            ModuleNeighbor::PagedNeighborList pages;

            const long long storage_us = elapsed_microseconds(
                [&]() { storage = ModuleNeighbor::build_grid_storage_from_atom_pack(pack, box_edge); });
            const long long search_us = elapsed_microseconds([&]() {
                pairs = options.mode == "full27"
                            ? ModuleNeighbor::build_neighbor_pairs_27(pack, storage, options.radius)
                            : ModuleNeighbor::build_neighbor_pairs_14(pack, storage, options.radius);
                if (options.mode == "half14_ref")
                {
                    reference = ModuleNeighbor::build_neighbor_pairs_27(pack, storage, options.radius);
                    if (pairs != reference)
                    {
                        throw std::runtime_error("Half14 result differs from the Full27 reference.");
                    }
                }
            });
            const long long page_us = elapsed_microseconds([&]() {
                pages.build(pairs, std::vector<int>{options.atoms});
            });

            std::cout << "neighbor_kernel," << options.mode << "," << options.atoms << "," << repeat << ","
                      << storage_us << "," << search_us << "," << page_us << "," << pairs.size() << ","
                      << storage.box_size() << "," << pages.page_count() << "," << pages.utilization() << ","
                      << pages.memory_usage_bytes() << "," << peak_rss_kib() << "," << pair_hash(pairs) << "\n";
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "neighbor benchmark failed: " << error.what() << "\n";
        return 1;
    }
}
