#include"../math_sphbes.h"
#include<fstream>
#include <benchmark/benchmark.h>
#include <iostream>
#include <cstring>
#include <cmath>

/************************************************
*  performace test of class Sphbes
***********************************************/

/**
 * Tested function: 
 *      - sphbesj
 *      - Spherical_Bessel
 */

class PerfSphbes : public benchmark::Fixture {
public:
    const double q = 1;
    const int n = 1000;
    double stop = 1000.0;
    double dr = 0.0;
    double* rc, *rinf, *jc, *jinf;
    void SetUp(const benchmark::State& state){
        const double rcut = state.range(0) + 0.5;
        rc = new double[n + 10]; 
        rinf = new double[n + 10];
        jc = new double[n + 10];
        jinf = new double[n + 10];

        // generate data points in (0, rcut] in log scale
        double rmin = 0.0001;
        double log_rmin = std::log(rmin);
        double log_rcut = std::log(rcut);
        dr = (log_rcut - log_rmin) / (n-1);
        memset(rc, 0, (n+10) * sizeof(double));
        for (int i = 0; i < n; i++)
            rc[i] = std::exp(log_rmin + i * dr);
        
        // generate data points in [rcut, stop] in linear scale
        memset(rinf, 0, (n+10) * sizeof(double));
        rinf[0] = rcut;
        dr = (stop - rcut) / (n-1);
        for (int i = 1; i < n; i++)
            rinf[i] += rinf[i-1] + dr;
    }
    void TearDown(const benchmark::State& state){
        delete[] rc;
        delete[] rinf;
        delete[] jc;
        delete[]  jinf;
    }
};    

BENCHMARK_DEFINE_F(PerfSphbes, BM_Spherical_Bessel)(benchmark::State& state) {
    for (auto _ : state) {
        ModuleBase::Sphbes::Spherical_Bessel(n, rc, q, state.range(0), jc);
        ModuleBase::Sphbes::Spherical_Bessel(n, rinf, q, state.range(0), jinf);
    }
}

BENCHMARK_DEFINE_F(PerfSphbes, BM_sphbesj)(benchmark::State& state) {
    for (auto _ : state) {
        ModuleBase::Sphbes::sphbesj(n, rc, q, state.range(0), jc);
        ModuleBase::Sphbes::sphbesj(n, rinf, q, state.range(0), jinf);
    }
}

BENCHMARK_REGISTER_F(PerfSphbes, BM_sphbesj)->DenseRange(0, 11, 1)->Unit(benchmark::kMicrosecond);
BENCHMARK_REGISTER_F(PerfSphbes, BM_Spherical_Bessel)->DenseRange(0, 11, 1)->Unit(benchmark::kMicrosecond);
BENCHMARK_MAIN(); 