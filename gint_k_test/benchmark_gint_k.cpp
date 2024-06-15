#include <benchmark/benchmark.h>
#include "gint_k.h"

static void BM_GintKIntegration(benchmark::State& state) {
    GintK gint;
    gint.initializeParameters();
    for (auto _ : state) {
        gint.computeIntegration();
    }
}
BENCHMARK(BM_GintKIntegration);

BENCHMARK_MAIN();

