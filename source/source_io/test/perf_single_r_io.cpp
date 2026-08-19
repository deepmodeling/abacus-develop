#include "benchmark/benchmark.h"

#include "source_io/module_hs/single_r_io.h"

#include <cstdint>
#include <cstdio>
#include <fstream>

namespace
{
void benchmark_single_R(benchmark::State& state)
{
    const int nlocal = static_cast<int>(state.range(0));
    const int entries_per_row = static_cast<int>(state.range(1));
    Parallel_Orbitals pv;
    pv.set_serial(nlocal, nlocal);

    ModuleIO::SparseRBlock<double> matrix;
    for (int row = 0; row < nlocal; ++row)
    {
        for (int offset = 0; offset < entries_per_row; ++offset)
        {
            const int col = (row + offset * 17) % nlocal;
            matrix[row][col] = 1.0 + 0.01 * offset;
        }
    }

    ModuleIO::SparseWriteOptions options;
    options.binary = true;
    options.reduce = false;
    options.threshold = 1e-10;
    const char* output_filename = "/tmp/abacus_perf_single_r_io.dat";

    for (auto _: state)
    {
        std::ofstream output(output_filename, std::ios::binary | std::ios::trunc);
        ModuleIO::output_single_R(output, matrix, pv, options);
        output.close();
        benchmark::ClobberMemory();
    }
    std::remove(output_filename);
    state.SetItemsProcessed(state.iterations() * nlocal * entries_per_row);
    state.SetBytesProcessed(state.iterations() * nlocal * entries_per_row
                            * static_cast<int64_t>(sizeof(double) + sizeof(int)));
}
} // namespace

BENCHMARK(benchmark_single_R)
    ->Args({2000, 8})
    ->Args({10000, 8})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
